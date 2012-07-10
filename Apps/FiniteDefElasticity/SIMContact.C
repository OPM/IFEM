// $Id$
//==============================================================================
//!
//! \file SIMContact.C
//!
//! \date Jun 05 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Administration of finite deformation contact analysis.
//!
//==============================================================================

#include "SIMContact.h"
#include "MortarContact.h"
#include "RigidBody.h"
#include "ASMbase.h"
#include "Utilities.h"
#include "Functions.h"
#include "Vec3Oper.h"
#include "MPC.h"
#include "tinyxml.h"


SIMContact::~SIMContact ()
{
  size_t i;
  for (i = 0; i < myBodies.size(); i++)
    delete myBodies[i];
  for (i = 0; i < myFuncs.size(); i++)
    delete myFuncs[i];
}


bool SIMContact::addContactElms (std::vector<ContactPair>& contacts,
				 RigidBody* master,
				 const std::string& slave,
				 const TopologySet& entitys,
				 const std::vector<ASMbase*>& model)
{
  TopologySet::const_iterator tit = entitys.find(slave);
  if (tit == entitys.end())
  {
    std::cerr <<" *** SIMContact::addContactElms: Undefined topology set \""
              << slave <<"\""<< std::endl;
    return false;
  }

  // Create contact elements connecting the slave boundary
  // to the master nodes of the rigid body coming in contact
  size_t nMast = master->getNoNodes();
  std::vector<int> nodes;
  nodes.reserve(nMast);

  TopEntity::const_iterator top;
  for (top = tit->second.begin(); top != tit->second.end(); top++)
    if (top->patch > 0 && top->patch <= model.size())
    {
      if (!model[top->patch-1]->addXElms(top->idim,top->item,nMast,nodes))
      {
        std::cerr <<" *** SIMContact::addContactElms: Invalid slave definition."
                  << std::endl;
        return false;
      }
      contacts.push_back(std::make_pair(master,model[top->patch-1]));
    }

  // Assign global master node numbers
  master->initNodes(nodes);
  return true;
}


/*!
  When we use the Augmented Lagrange formulation, a set of Lagrange multipliers
  is added as DOFs in the FE model. One multiplier (representing the contact
  force) is added for each FE node on the slave boundary. The boundary nodes
  are identified through a positive node number in the element connectivity of
  the "extra-ordinary" elements of the patch (the non-boundary nodes have a
  negative number for these elements).
*/

bool SIMContact::addLagrangeMultipliers (ASMbase* patch, int& ngnod)
{
  if (contactMethod != AUGMENTED_LAGRANGE) return true;

  size_t nXelm = patch->getNoXelms();
  if (nXelm == 0) return true;

  // Loop over the extra-ordinary (contact) elements in the patch
  size_t nno = patch->getNoNodes();
  size_t iel = patch->getNoElms(true) - nXelm;
  IntMat::const_iterator eit = patch->begin_elm() + iel;
  for (++iel; eit != patch->end_elm(); eit++, iel++)
  {
    IntVec mlagel;
    for (IntVec::const_iterator nit = eit->begin(); nit != eit->end(); nit++)
      if (*nit >= 0 && patch->getNodeType(*nit+1) == 'D')
      {
        // This is a boundary node, assign a Lagrange multiplier to it
        int node = patch->getNodeID(*nit+1);
        std::map<int,int>::const_iterator it = myALp.find(node);
        if (it == myALp.end())
          myALp[node] = ++ngnod;
        mlagel.push_back(myALp[node]);
      }

    if (!mlagel.empty())
      if (!patch->addLagrangeMultipliers(iel,mlagel))
        return false;
  }

  int nLag = patch->getNoNodes() - nno;
  if (nLag > 0)
    std::cout <<"\nAdded "<< nLag <<" Lagrange multipliers to patch "
              << patch->idx+1 << std::endl;

  return true;
}


bool SIMContact::parseContactTag (const TiXmlElement* elem,
				  const std::vector<ASMbase*>& model,
				  const TopologySet& entitys)
{
  if (contactMethod != NONE)
  {
    std::cerr <<" *** SIMContact::parseContactTag: Duplicated contact tags."
              << std::endl;
    return false;
  }

  // Determine the contact formulation to use
  std::string formulation;
  utl::getAttribute(elem,"formulation",formulation,true);
  if (formulation == "penalty")
    contactMethod = PENALTY;
  else
    contactMethod = AUGMENTED_LAGRANGE;

  double R = 1.0;
  RigidBody* body = NULL;
  std::vector<ContactPair> cset;

  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
  {
    // Define a new rigid master body
    utl::getAttribute(child,"R",R);
    if (!strcasecmp(child->Value(),"sphere"))
      body = new RigidSphere(R);
    else if (!strcasecmp(child->Value(),"cylinder"))
      body = new RigidCylinder(R,model.front()->getNoSpaceDim());
    else if (!strcasecmp(child->Value(),"plane"))
      body = new RigidPlane(model.front()->getNoSpaceDim());
    else
      continue;

    size_t p = 0, nMast = body->getNoNodes();
    std::vector<Vec3> points(nMast);
    std::string slaveSet;

    const TiXmlElement* c = child->FirstChildElement();
    for (; c; c = c->NextSiblingElement())
    {
      const TiXmlNode* cval = c->FirstChild();
      if (!strcasecmp(c->Value(),"point") && cval && p < nMast)
      {
	// Coordinates of the master nodes
	std::istringstream value(cval->Value());
	value >> points[p++];
      }
      else if (!strcasecmp(c->Value(),"slave"))
      {
	// Identification of the slave boundary on the flexible body
	utl::getAttribute(c,"set",slaveSet);
	if (!this->createContactSet(slaveSet,body->code))
	  return false;
	std::cout <<"\tContact code "<< body->code << std::endl;
	if (!this->addContactElms(cset,body,slaveSet,entitys,model))
	  return false;
      }
      else if (!strcasecmp(c->Value(),"dirichlet"))
      {
	// Define the Dirichlet properties of the rigid body
	int comp = 0;
	std::string type;
	utl::getAttribute(c,"comp",comp);
	utl::getAttribute(c,"type",type,true);
	if (!cval || (atof(cval->Value()) == 0.0 && type != "expression"))
	{
	  std::cout <<"\tFixed contact body"<< std::endl;
	  type = "fixed";
	}
	else if (cval && comp > 0 && comp <= 3)
	{
	  std::cout <<"\tPrescribed contact body";
	  myFuncs.push_back(utl::parseTimeFunc(cval->Value(),type));
	  type = "prescribed";
	}
	else
	  continue;

	// Map the Dirichlet properties onto the master nodes
	std::vector<ContactPair>::iterator sit = cset.begin();
	while ((sit = std::find_if(sit,cset.end(),hasBody(body))) != cset.end())
	{
	  int iMast = sit->second->getNoNodes() - nMast;
	  for (size_t i = 1; i <= nMast; i++)
	    if (type == "fixed")
	      sit->second->fix(++iMast,comp);
	    else
	    {
	      sit->second->prescribe(++iMast,comp,-1);
	      MPC* spc = sit->second->findMPC(body->getNodeID(i),comp);
	      dMap[spc] = myFuncs.back();
	    }
	  ++sit;
	}
      }
      else
      {
	// Get the penalty parameter
	const char* value = NULL;
	if ((value = utl::getValue(c,"eps")))
	  body->eps = atof(value);
      }
    }

    body->initPoints(points);
    body->print(std::cout);
    myBodies.push_back(body);
  }

  return true;
}


bool SIMContact::preprocessContact (IntegrandMap& itgs,
				    const SAM& sam, size_t nsd)
{
  // Count the total number of master nodes
  size_t nMaster = 0;
  std::vector<RigidBody*>::const_iterator it;
  for (it = myBodies.begin(); it != myBodies.end(); it++)
    nMaster += (*it)->getNoNodes();

  // Establish the integrands of the Mortar-based contact analysis
  IntegrandBase* contactInt = 0;
  for (it = myBodies.begin(); it != myBodies.end(); it++)
    if ((*it)->eps > 0.0)
    {
      MortarMats* mats = new MortarMats(sam,nMaster,nsd);
      if (contactMethod == PENALTY)
	contactInt = new MortarPenalty(*it,*mats,nsd);
      else if (contactMethod == AUGMENTED_LAGRANGE)
	contactInt = new MortarAugmentedLag(*it,*mats,myALp,nsd);
      else
      {
	delete mats;
	return false;
      }

      // Integrand for the Mortar matrices
      itgs.insert(std::make_pair((*it)->code,new MortarContact(*it,mats,nsd)));
      // Integrand for the residual and tangent contributions
      itgs.insert(std::make_pair((*it)->code,contactInt));
    }

  return true;
}


void SIMContact::updateDirichlet (double time)
{
  std::map<MPC*,ScalarFunc*>::iterator cit;
  for (cit = dMap.begin(); cit != dMap.end(); cit++)
    if (cit->first && cit->second)
      cit->first->setSlaveCoeff((*cit->second)(time));
}


bool SIMContact::updateContactBodies (const std::vector<double>& displ)
{
  for (size_t i = 0; i < myBodies.size(); i++)
    if (!myBodies[i]->update(displ))
      return false;

  return true;
}


bool SIMContact::assembleMortarTangent (const IntegrandBase* problem,
					SystemMatrix* Ktan, SystemVector* Res)
{
  const MortarContact* contp = dynamic_cast<const MortarContact*>(problem);
  if (!contp || !Ktan || !Res) return true; // silently ignore

#if INT_DEBUG > 2
  // Must assemble into a temporary matrix to see the explicit contributions
  SparseMatrix A(Ktan->dim());
  StdVector B(Res->size());
  if (contp->assemble(A,B) && A.size() > 0)
    std::cout <<"SIMContact::assembleMortarTangent: Kmain "<< A;
#endif

  return contp->assemble(*Ktan,*Res);
}
