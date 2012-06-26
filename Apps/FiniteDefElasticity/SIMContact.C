// $Id$
//==============================================================================
//!
//! \file SIMContact.C
//!
//! \date Jun 05 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution drivers for finite deformation contact analysis.
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
  for (size_t i = 0; i < myBodies.size(); i++)
    delete myBodies[i];
}


bool SIMContact::addContactElms (RigidBody* master,
				 const std::vector<ASMbase*>& model,
				 const TopologySet& entitys,
				 const std::string& slave)
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
      myContacts.push_back(std::make_pair(master,model[top->patch-1]));
    }

  // Assign global master node numbers
  master->initNodes(nodes);
  return true;
}


bool SIMContact::parseContactTag (const TiXmlElement* elem,
				  const std::vector<ASMbase*>& model,
				  const TopologySet& entitys)
{
  double R = 1.0;
  RigidBody* body = NULL;
  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
  {
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

    const char* value = 0;
    const TiXmlElement* c = child->FirstChildElement();
    for (; c; c = c->NextSiblingElement())
    {
      const TiXmlNode* cval = c->FirstChild();
      if (!strcasecmp(c->Value(),"point") && cval && p < nMast)
      {
	std::istringstream value(cval->Value());
	value >> points[p++];
      }
      else if (!strcasecmp(c->Value(),"slave"))
      {
	utl::getAttribute(c,"set",slaveSet);
	if (!this->createContactSet(slaveSet,body->code))
	  return false;
	std::cout <<"\tContact code "<< body->code << std::endl;
	if (!this->addContactElms(body,model,entitys,slaveSet))
	  return false;
      }
      else if ((value = utl::getValue(c,"penalty")))
	body->eps = atof(value);

      else if (!strcasecmp(c->Value(),"dirichlet"))
      {
	ASMbase* slave = NULL;
	for (size_t i = 0; i < myContacts.size() && !slave; i++)
	  if (myContacts[i].first == body) slave = myContacts[i].second;

	if (!slave)
	{
	  std::cerr <<" *** SIMContact::parseContactTag: No slave definition."
		    <<"\n                                  <slave> must be "
		    <<"specified before <dirichlet>."<< std::endl;
	  return false;
	}

	int comp = 0;
	std::string type;
	utl::getAttribute(c,"comp",comp);
	utl::getAttribute(c,"type",type,true);
	if (!cval || (atof(cval->Value()) == 0.0 && type != "expression"))
	{
	  std::cout <<"\tFixed contact body"<< std::endl;
	  for (size_t i = 1; i <= nMast; i++)
	    slave->fix(slave->getNoNodes()+i-nMast,comp);
	}
	else if (cval && comp > 0 && comp <= 3)
	{
	  std::cout <<"\tPrescribed contact body";
	  if (type == "expression") std::cout <<" (expression)";
	  char* cstr = strdup(cval->Value());
	  const ScalarFunc* sf = utl::parseTimeFunc(type.c_str(),cstr);
	  std::cout << std::endl;
	  free(cstr);
	  myFuncs.push_back(const_cast<ScalarFunc*>(sf));
	  for (size_t i = 1; i <= nMast; i++)
	  {
	    slave->prescribe(slave->getNoNodes()+i-nMast,comp,-1);
	    dMap[slave->findMPC(body->getNodeID(i),comp)] = myFuncs.back();
	  }
	}
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
  for (it = myBodies.begin(); it != myBodies.end(); it++)
    if ((*it)->eps > 0.0)
    {
      int code = (*it)->code;
      MortarMats* mats = new MortarMats(sam,nMaster,nsd);
      itgs.insert(std::make_pair(code,new MortarContact(*it,mats,nsd)));
      itgs.insert(std::make_pair(code,new MortarPenalty(*it,*mats,nsd)));
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
