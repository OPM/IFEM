// $Id$
//==============================================================================
//!
//! \file SIM1D.C
//!
//! \date Feb 04 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for 1D NURBS-based FEM analysis.
//!
//==============================================================================

#include "SIM1D.h"
#include "ModelGenerator.h"
#include "ASMs1D.h"
#include "Functions.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "IFEM.h"
#include "tinyxml.h"
#include <sstream>


SIM1D::SIM1D (unsigned char n1, bool)
{
  nsd = 1;
  nf = n1;
  twist = nullptr;
}


SIM1D::SIM1D (const CharVec& fields, bool)
{
  nsd = 1;
  nf = fields.empty()?1:fields.front();
  twist = nullptr;
}


SIM1D::SIM1D (IntegrandBase* itg, unsigned char n) : SIMgeneric(itg)
{
  nsd = 1;
  nf = n;
  twist = nullptr;
}


bool SIM1D::addConnection (int master, int slave, int mIdx, int sIdx,
                           int, int basis, bool, int dim, int thick)
{
  if (basis > 0)
  {
    std::cerr <<" *** SIM1D::addConnection: Mixed not implemented."<< std::endl;
    return false;
  }

  int lmaster = this->getLocalPatchIndex(master);
  int lslave = this->getLocalPatchIndex(slave);

  if (lmaster > 0 && lslave > 0)
  {
    if (dim != 0) return false;

    IFEM::cout <<"\tConnecting P"<< slave <<" V"<< sIdx
               <<" to P"<< master <<" V"<< mIdx
               <<" basis="<< basis <<" thick="<< thick << std::endl;

    ASMs1D* spch = static_cast<ASMs1D*>(myModel[lslave-1]);
    ASMs1D* mpch = static_cast<ASMs1D*>(myModel[lmaster-1]);

    if (!spch->connectPatch(sIdx,*mpch,mIdx,thick))
      return false;
  }
  else
    adm.dd.ghostConnections.insert(DomainDecomposition::Interface{master, slave,
                                                                  mIdx, sIdx, 0,
                                                                  dim, basis, thick});

  return true;
}


bool SIM1D::parseGeometryTag (const TiXmlElement* elem)
{
  IFEM::cout <<"  Parsing <"<< elem->Value() <<">"<< std::endl;

  if (!strcasecmp(elem->Value(),"refine"))
  {
    IntVec patches;
    if (!this->parseTopologySet(elem,patches))
      return false;

    ASM1D* pch = nullptr;
    RealArray xi;
    if (!utl::parseKnots(elem,xi))
    {
      int addu = 0;
      utl::getAttribute(elem,"u",addu);
      for (int j : patches)
        if ((pch = dynamic_cast<ASM1D*>(this->getPatch(j,true))))
        {
          IFEM::cout <<"\tRefining P"<< j <<" "<< addu << std::endl;
          pch->uniformRefine(addu);
        }
    }
    else
      for (int j : patches)
        if ((pch = dynamic_cast<ASM1D*>(this->getPatch(j,true))))
        {
          IFEM::cout <<"\tRefining P"<< j
                     <<" with grading "<< elem->FirstChild()->Value() <<":";
          for (size_t i = 0; i < xi.size(); i++)
            IFEM::cout << (i%10 || xi.size() < 11 ? " " : "\n\t") << xi[i];
          IFEM::cout << std::endl;
          pch->refine(xi);
        }
  }

  else if (!strcasecmp(elem->Value(),"raiseorder"))
  {
    IntVec patches;
    if (!this->parseTopologySet(elem,patches))
      return false;

    ASM1D* pch = nullptr;
    int addu = 0;
    utl::getAttribute(elem,"u",addu);
    for (int j : patches)
      if ((pch = dynamic_cast<ASM1D*>(this->getPatch(j,true))))
      {
        IFEM::cout <<"\tRaising order of P"<< j <<" "<< addu << std::endl;
        static_cast<ASM1D*>(pch)->raiseOrder(addu);
      }
  }

  else if (!strcasecmp(elem->Value(),"topology"))
  {
    if (!this->createFEMmodel()) return false;

    const TiXmlElement* child = elem->FirstChildElement("connection");
    for (; child; child = child->NextSiblingElement())
    {
      int master = 0, slave = 0, mVert = 0, sVert = 0;
      utl::getAttribute(child,"master",master);
      utl::getAttribute(child,"mvert",mVert);
      utl::getAttribute(child,"slave",slave);
      utl::getAttribute(child,"svert",sVert);

      if (master == slave ||
          master < 1 || master > (int)myModel.size() ||
          slave  < 1 || slave  > (int)myModel.size())
      {
        std::cerr <<" *** SIM1D::parse: Invalid patch indices "
                  << master <<" "<< slave << std::endl;
        return false;
      }
      IFEM::cout <<"\tConnecting P"<< slave <<" V"<< sVert
                 <<" to P"<< master <<" E"<< mVert << std::endl;
      ASMs1D* spch = static_cast<ASMs1D*>(myModel[slave-1]);
      ASMs1D* mpch = static_cast<ASMs1D*>(myModel[master-1]);
      if (!spch->connectPatch(sVert,*mpch,mVert))
        return false;
    }
  }

  else if (!strcasecmp(elem->Value(),"periodic"))
  {
    if (!this->createFEMmodel()) return false;

    int patch = 0;
    utl::getAttribute(elem,"patch",patch);

    if (patch < 1 || patch > (int)myModel.size())
    {
      std::cerr <<" *** SIM1D::parse: Invalid patch index "
                << patch << std::endl;
      return false;
    }
    IFEM::cout <<"\tPeriodic P"<< patch << std::endl;
    static_cast<ASMs1D*>(myModel[patch-1])->closeEnds();
  }

  else if (!strcasecmp(elem->Value(),"Zdirection"))
  {
    utl::getAttribute(elem,"x",XZp.x);
    utl::getAttribute(elem,"y",XZp.y);
    utl::getAttribute(elem,"z",XZp.z);
    IFEM::cout <<"\tZ-direction vector: "<< XZp << std::endl;
  }

  return true;
}


bool SIM1D::parseBCTag (const TiXmlElement* elem)
{
  if (!strcasecmp(elem->Value(),"fixpoint") && !ignoreDirichlet)
  {
    if (!this->createFEMmodel()) return false;

    int patch = 0, code = 123;
    double rx = 0.0;
    utl::getAttribute(elem,"patch",patch);
    utl::getAttribute(elem,"code",code);
    utl::getAttribute(elem,"rx",rx);
    int pid = this->getLocalPatchIndex(patch);
    if (pid < 1) return pid == 0;

    ASM1D* pch = dynamic_cast<ASM1D*>(myModel[pid-1]);
    if (!pch) return false;

    IFEM::cout <<"\tConstraining P"<< patch
               <<" point at "<< rx <<" with code "<< code << std::endl;
    pch->constrainNode(rx,code);
  }

  return true;
}


/*!
  The twist angle is used to define the local element axes of beam elements
  along the spline curves. The angle described the rotation of the local
  Y-axis, relative to the globalized Y-axis of the beam.
*/

bool SIM1D::parseTwist (const TiXmlElement* elem)
{
  if (!elem->FirstChild())
    return false;

  std::string type;
  utl::getAttribute(elem,"type",type);
  IFEM::cout <<"    Continuous twist angle:";
  if (!type.empty()) IFEM::cout <<" ("<< type <<")";
  twist = utl::parseRealFunc(elem->FirstChild()->Value(),type);
  IFEM::cout << std::endl;

  return true;
}


bool SIM1D::parse (const TiXmlElement* elem)
{
  bool result = this->SIMgeneric::parse(elem);

  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(elem->Value(),"geometry"))
      result &= this->parseGeometryTag(child);
    else if (!strcasecmp(elem->Value(),"boundaryconditions"))
      result &= this->parseBCTag(child);

  if (myGen && result)
    result = myGen->createTopology(*this);

  delete myGen;
  myGen = nullptr;

  return result;
}


bool SIM1D::parse (char* keyWord, std::istream& is)
{
  char* cline = nullptr;
  if (!strncasecmp(keyWord,"REFINE",6))
  {
    int nref = atoi(keyWord+6);
    IFEM::cout <<"\nNumber of patch refinements: "<< nref << std::endl;
    for (int i = 0; i < nref && (cline = utl::readLine(is)); i++)
    {
      bool uniform = !strchr(cline,'.');
      int patch = atoi(strtok(cline," "));
      if (patch == 0 || abs(patch) > (int)myModel.size())
      {
	std::cerr <<" *** SIM1D::parse: Invalid patch index "
		  << patch << std::endl;
	return false;
      }
      int ipatch = patch-1;
      if (patch < 0)
      {
	ipatch = 0;
	patch = -patch;
      }
      if (uniform)
      {
	int addu = atoi(strtok(nullptr," "));
	for (int j = ipatch; j < patch; j++)
	{
	  IFEM::cout <<"\tRefining P"<< j+1 <<" "<< addu << std::endl;
	  static_cast<ASMs1D*>(myModel[j])->uniformRefine(addu);
	}
      }
      else
      {
	RealArray xi;
	if (utl::parseKnots(xi))
	  for (int j = ipatch; j < patch; j++)
	  {
	    IFEM::cout <<"\tRefining P"<< j+1;
	    for (size_t i = 0; i < xi.size(); i++)
	      IFEM::cout <<" "<< xi[i];
	    IFEM::cout << std::endl;
	    static_cast<ASMs1D*>(myModel[j])->refine(xi);
	  }
      }
    }
  }

  else if (!strncasecmp(keyWord,"RAISEORDER",10))
  {
    int nref = atoi(keyWord+10);
    IFEM::cout <<"\nNumber of order raise: "<< nref << std::endl;
    for (int i = 0; i < nref && (cline = utl::readLine(is)); i++)
    {
      int patch = atoi(strtok(cline," "));
      int addu  = atoi(strtok(nullptr," "));
      if (patch == 0 || abs(patch) > (int)myModel.size())
      {
	std::cerr <<" *** SIM1D::parse: Invalid patch index "
		  << patch << std::endl;
	return false;
      }
      int ipatch = patch-1;
      if (patch < 0)
      {
	ipatch = 0;
	patch = -patch;
      }
      for (int j = ipatch; j < patch; j++)
      {
	IFEM::cout <<"\tRaising order of P"<< j+1
		   <<" "<< addu << std::endl;
	static_cast<ASMs1D*>(myModel[j])->raiseOrder(addu);
      }
    }
  }

  else if (!strncasecmp(keyWord,"TOPOLOGY",8))
  {
    if (!this->createFEMmodel()) return false;

    int ntop = atoi(keyWord+8);
    IFEM::cout <<"\nNumber of patch connections: "<< ntop << std::endl;
    for (int i = 0; i < ntop && (cline = utl::readLine(is)); i++)
    {
      int master = atoi(strtok(cline," "));
      int mVert  = atoi(strtok(nullptr," "));
      int slave  = atoi(strtok(nullptr," "));
      int sVert  = atoi(strtok(nullptr," "));
      if (master == slave ||
	  master < 1 || master > (int)myModel.size() ||
	  slave  < 1 || slave  > (int)myModel.size())
      {
	std::cerr <<" *** SIM1D::parse: Invalid patch indices "
		  << master <<" "<< slave << std::endl;
	return false;
      }
      IFEM::cout <<"\tConnecting P"<< slave <<" V"<< sVert
                 <<" to P"<< master <<" E"<< mVert << std::endl;
      ASMs1D* spch = static_cast<ASMs1D*>(myModel[slave-1]);
      ASMs1D* mpch = static_cast<ASMs1D*>(myModel[master-1]);
      if (!spch->connectPatch(sVert,*mpch,mVert))
	return false;
    }
  }

  else if (!strncasecmp(keyWord,"PERIODIC",8))
  {
    if (!this->createFEMmodel()) return false;

    int nper = atoi(keyWord+8);
    IFEM::cout <<"\nNumber of periodicities: "<< nper << std::endl;
    for (int i = 0; i < nper && (cline = utl::readLine(is)); i++)
    {
      int patch = atoi(strtok(cline," "));
      if (patch < 1 || patch > (int)myModel.size())
      {
	std::cerr <<" *** SIM1D::parse: Invalid patch index "
		  << patch << std::endl;
	return false;
      }
      IFEM::cout <<"\tPeriodic P"<< patch << std::endl;
      static_cast<ASMs1D*>(myModel[patch-1])->closeEnds();
    }
  }

  else if (!strncasecmp(keyWord,"CONSTRAINTS",11))
  {
    if (ignoreDirichlet) return true; // Ignore all boundary conditions
    if (!this->createFEMmodel()) return false;

    int ngno = 0;
    int ncon = atoi(keyWord+11);
    IFEM::cout <<"\nNumber of constraints: "<< ncon << std::endl;
    for (int i = 0; i < ncon && (cline = utl::readLine(is)); i++)
    {
      int patch = atoi(strtok(cline," "));
      int pvert = atoi(strtok(nullptr," "));
      int bcode = atoi(strtok(nullptr," "));
      double pd = (cline = strtok(nullptr," ")) ? atof(cline) : 0.0;
      if (pd == 0.0)
      {
	if (!this->addConstraint(patch,pvert,0,bcode%1000000,0,ngno))
	  return false;
      }
      else
      {
	int code = 1000000 + bcode;
	while (myScalars.find(code) != myScalars.end())
	  code += 1000000;

	if (!this->addConstraint(patch,pvert,0,bcode%1000000,code,ngno))
	  return false;

	IFEM::cout <<" ";
	cline = strtok(nullptr," ");
	myScalars[code] = const_cast<RealFunc*>(utl::parseRealFunc(cline,pd));
      }
      IFEM::cout << std::endl;
    }
  }

  else if (!strncasecmp(keyWord,"FIXPOINTS",9))
  {
    if (ignoreDirichlet) return true; // Ignore all boundary conditions
    if (!this->createFEMmodel()) return false;

    int nfix = atoi(keyWord+9);
    IFEM::cout <<"\nNumber of fixed points: "<< nfix << std::endl;
    for (int i = 0; i < nfix && (cline = utl::readLine(is)); i++)
    {
      int patch = atoi(strtok(cline," "));
      double rx = atof(strtok(nullptr," "));
      int bcode = (cline = strtok(nullptr," ")) ? atoi(cline) : 123;
      if (patch < 1 || patch > (int)myModel.size())
      {
	std::cerr <<" *** SIM1D::parse: Invalid patch index "
		  << patch << std::endl;
	return false;
      }
      IFEM::cout <<"\tConstraining P"<< patch
                 <<" point at "<< rx <<" with code "<< bcode << std::endl;
      static_cast<ASMs1D*>(myModel[patch-1])->constrainNode(rx,bcode);
    }
  }

  else
    return this->SIMgeneric::parse(keyWord,is);

  return true;
}


/*!
  \brief Local-scope convenience function for error message generation.
*/

static bool constrError (const char* lab, int idx)
{
  std::cerr <<" *** SIM1D::addConstraint: Invalid "<< lab << idx << std::endl;
  return false;
}


bool SIM1D::addConstraint (int patch, int lndx, int ldim, int dirs, int code,
                           int&, char basis)
{
  if (patch < 1 || patch > (int)myModel.size())
    return constrError("patch index ",patch);

  IFEM::cout <<"\tConstraining P"<< patch;
  if (ldim == 0)
    IFEM::cout << " V"<< lndx;
  IFEM::cout <<" in direction(s) "<< dirs;
  if (code != 0) IFEM::cout <<" code = "<< code;
  if (basis > 1) IFEM::cout <<" basis = "<< (int)basis;
#if SP_DEBUG > 1
  std::cout << std::endl;
#endif

  ASMs1D* pch = static_cast<ASMs1D*>(myModel[patch-1]);
  if (abs(ldim) > 0)
    pch->constrainPatch(dirs,code);
  else switch (lndx) // Vertex constraints
    {
    case 1: pch->constrainNode(0.0,dirs,code,basis); break;
    case 2: pch->constrainNode(1.0,dirs,code,basis); break;
    default:
      IFEM::cout << std::endl;
      return constrError("vertex index ",lndx);
    }

  return true;
}


ASMbase* SIM1D::readPatch (std::istream& isp, int pchInd,
                           const CharVec& unf) const
{
  ASMbase* pch = ASM1D::create(opt.discretization,nsd,
                               unf.empty() ? nf : unf.front());
  if (pch)
  {
    if (!pch->read(isp))
      delete pch, pch = nullptr;
    else if (pch->empty() || this->getLocalPatchIndex(pchInd+1) < 1)
      delete pch, pch = nullptr;
    else
      pch->idx = myModel.size();
  }

  return pch;
}


bool SIM1D::readPatches (std::istream& isp, PatchVec& patches,
                         const char* whiteSpace) const
{
  ASMbase* pch = nullptr;
  for (int pchInd = 1; isp.good(); pchInd++)
    if ((pch = ASM1D::create(opt.discretization,nsd,nf)))
    {
      if (whiteSpace)
        IFEM::cout << whiteSpace <<"Reading patch "<< pchInd << std::endl;
      if (!pch->read(isp))
      {
        delete pch;
        return false;
      }
      else if (pch->empty() || this->getLocalPatchIndex(pchInd) < 1)
        delete pch;
      else
      {
        pch->idx = patches.size();
        patches.push_back(pch);
      }
    }

  return true;
}


bool SIM1D::createFEMmodel (char)
{
  bool ok = true;
  ASMstruct::resetNumbering();
  for (size_t i = 0; i < myModel.size() && ok; i++)
    if (twist)
      ok = static_cast<ASMs1D*>(myModel[i])->generateTwistedFEModel(*twist,XZp);
    else if (!XZp.isZero())
      ok = static_cast<ASMs1D*>(myModel[i])->generateOrientedFEModel(XZp);
    else
      ok = myModel[i]->generateFEMTopology();

  return ok;
}


ModelGenerator* SIM1D::getModelGenerator (const TiXmlElement* geo) const
{
  return new DefaultGeometry1D(geo);
}


Vector SIM1D::getSolution (const Vector& psol, double u,
                           int deriv, int patch) const
{
  return this->SIMgeneric::getSolution(psol,&u,deriv,patch);
}


bool SIM1D::updateRotations (const Vector& incSol, double alpha)
{
  if (nf != 6) return true;

  bool ok = true;
  for (size_t i = 0; i < myModel.size() && ok; i++)
    if (incSol.empty())
      // Update the rotations of last converged load/time step
      static_cast<ASMs1D*>(myModel[i])->updateRotations();
    else
    {
      Vector locSol;
      myModel[i]->extractNodeVec(incSol,locSol);
      if (alpha != 0.0 && alpha != 1.0) locSol *= alpha;
      ok = static_cast<ASMs1D*>(myModel[i])->updateRotations(locSol,alpha==0.0);
    }

  return ok;
}
