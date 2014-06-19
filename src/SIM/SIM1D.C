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
#include "ASMs1D.h"
#include "Functions.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "tinyxml.h"
#include <sstream>


SIM1D::SIM1D (unsigned char n1, unsigned char, bool)
{
  nd = 1;
  nf = n1;
  twist = NULL;
}


SIM1D::SIM1D (IntegrandBase* itg, unsigned char n) : SIMgeneric(itg)
{
  nd = 1;
  nf = n;
  twist = NULL;
}


bool SIM1D::parseGeometryTag (const TiXmlElement* elem)
{
  std::cout <<"  Parsing <"<< elem->Value() <<">"<< std::endl;

  if (!strcasecmp(elem->Value(),"refine"))
  {
    int lowpatch = 1, uppatch = 2;
    if (utl::getAttribute(elem,"patch",lowpatch))
      uppatch = lowpatch;
    if (utl::getAttribute(elem,"lowerpatch",lowpatch))
      uppatch = myModel.size();
    utl::getAttribute(elem,"upperpatch",uppatch);

    if (lowpatch < 1 || uppatch > (int)myModel.size())
    {
      std::cerr <<" *** SIM1D::parse: Invalid patch indices, lower="
                << lowpatch <<" upper="<< uppatch << std::endl;
      return false;
    }

    RealArray xi;
    if (!utl::parseKnots(elem,xi))
    {
      int addu = 0;
      utl::getAttribute(elem,"u",addu);
      for (int j = lowpatch-1; j < uppatch; j++)
      {
        std::cout <<"\tRefining P"<< j+1 <<" "<< addu << std::endl;
        static_cast<ASMs1D*>(myModel[j])->uniformRefine(addu);
      }
    }
    else
      for (int j = lowpatch-1; j < uppatch; j++)
      {
        std::cout <<"\tRefining P"<< j+1;
        for (size_t i = 0; i < xi.size(); i++)
          std::cout <<" "<< xi[i];
        std::cout << std::endl;
        static_cast<ASMs1D*>(myModel[j])->refine(xi);
      }
  }

  else if (!strcasecmp(elem->Value(),"raiseorder"))
  {
    int lowpatch = 1, uppatch = 2;
    if (utl::getAttribute(elem,"patch",lowpatch))
      uppatch = lowpatch;
    if (utl::getAttribute(elem,"lowerpatch",lowpatch))
      uppatch = myModel.size();
    utl::getAttribute(elem,"upperpatch",uppatch);

    if (lowpatch < 1 || uppatch > (int)myModel.size())
    {
      std::cerr <<" *** SIM1D::parse: Invalid patch indices, lower="
                << lowpatch <<" upper="<< uppatch << std::endl;
      return false;
    }

    int addu = 0;
    utl::getAttribute(elem,"u",addu);
    for (int j = lowpatch-1; j < uppatch; j++)
    {
      std::cout <<"\tRaising order of P"<< j+1 <<" "<< addu << std::endl;
      static_cast<ASMs1D*>(myModel[j])->raiseOrder(addu);
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
      std::cout <<"\tConnecting P"<< slave <<" V"<< sVert
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
    std::cout <<"\tPeriodic P"<< patch << std::endl;
    static_cast<ASMs1D*>(myModel[patch-1])->closeEnds();
  }

  else if (!strcasecmp(elem->Value(),"Zdirection"))
  {
    utl::getAttribute(elem,"x",XZp.x);
    utl::getAttribute(elem,"y",XZp.y);
    utl::getAttribute(elem,"z",XZp.z);
    std::cout <<"\tZ-direction vector: "<< XZp << std::endl;
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

    std::cout <<"\tConstraining P"<< patch
              <<" point at "<< rx <<" with code "<< code << std::endl;
    static_cast<ASMs1D*>(myModel[patch-1])->constrainNode(rx,code);
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
  std::cout <<"    Continuous twist angle:";
  if (!type.empty()) std::cout <<" ("<< type <<")";
  twist = utl::parseRealFunc(elem->FirstChild()->Value(),type);
  std::cout << std::endl;

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

  return result;
}


bool SIM1D::parse (char* keyWord, std::istream& is)
{
  char* cline = NULL;
  if (!strncasecmp(keyWord,"REFINE",6))
  {
    int nref = atoi(keyWord+6);
    std::cout <<"\nNumber of patch refinements: "<< nref << std::endl;
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
	int addu = atoi(strtok(NULL," "));
	for (int j = ipatch; j < patch; j++)
	{
	  std::cout <<"\tRefining P"<< j+1 <<" "<< addu << std::endl;
	  static_cast<ASMs1D*>(myModel[j])->uniformRefine(addu);
	}
      }
      else
      {
	RealArray xi;
	if (utl::parseKnots(xi))
	  for (int j = ipatch; j < patch; j++)
	  {
	    std::cout <<"\tRefining P"<< j+1;
	    for (size_t i = 0; i < xi.size(); i++)
	      std::cout <<" "<< xi[i];
	    std::cout << std::endl;
	    static_cast<ASMs1D*>(myModel[j])->refine(xi);
	  }
      }
    }
  }

  else if (!strncasecmp(keyWord,"RAISEORDER",10))
  {
    int nref = atoi(keyWord+10);
    std::cout <<"\nNumber of order raise: "<< nref << std::endl;
    for (int i = 0; i < nref && (cline = utl::readLine(is)); i++)
    {
      int patch = atoi(strtok(cline," "));
      int addu  = atoi(strtok(NULL," "));
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
	std::cout <<"\tRaising order of P"<< j+1
		  <<" "<< addu << std::endl;
	static_cast<ASMs1D*>(myModel[j])->raiseOrder(addu);
      }
    }
  }

  else if (!strncasecmp(keyWord,"TOPOLOGY",8))
  {
    if (!this->createFEMmodel()) return false;

    int ntop = atoi(keyWord+8);
    std::cout <<"\nNumber of patch connections: "<< ntop << std::endl;
    for (int i = 0; i < ntop && (cline = utl::readLine(is)); i++)
    {
      int master = atoi(strtok(cline," "));
      int mVert  = atoi(strtok(NULL," "));
      int slave  = atoi(strtok(NULL," "));
      int sVert  = atoi(strtok(NULL," "));
      if (master == slave ||
	  master < 1 || master > (int)myModel.size() ||
	  slave  < 1 || slave  > (int)myModel.size())
      {
	std::cerr <<" *** SIM1D::parse: Invalid patch indices "
		  << master <<" "<< slave << std::endl;
	return false;
      }
      std::cout <<"\tConnecting P"<< slave <<" V"<< sVert
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
    std::cout <<"\nNumber of periodicities: "<< nper << std::endl;
    for (int i = 0; i < nper && (cline = utl::readLine(is)); i++)
    {
      int patch = atoi(strtok(cline," "));
      if (patch < 1 || patch > (int)myModel.size())
      {
	std::cerr <<" *** SIM1D::parse: Invalid patch index "
		  << patch << std::endl;
	return false;
      }
      std::cout <<"\tPeriodic P"<< patch << std::endl;
      static_cast<ASMs1D*>(myModel[patch-1])->closeEnds();
    }
  }

  else if (!strncasecmp(keyWord,"CONSTRAINTS",11))
  {
    if (ignoreDirichlet) return true; // Ignore all boundary conditions
    if (!this->createFEMmodel()) return false;

    int ngno = 0;
    int ncon = atoi(keyWord+11);
    std::cout <<"\nNumber of constraints: "<< ncon << std::endl;
    for (int i = 0; i < ncon && (cline = utl::readLine(is)); i++)
    {
      int patch = atoi(strtok(cline," "));
      int pvert = atoi(strtok(NULL," "));
      int bcode = atoi(strtok(NULL," "));
      double pd = (cline = strtok(NULL," ")) ? atof(cline) : 0.0;
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

	cline = strtok(NULL," ");
	myScalars[code] = const_cast<RealFunc*>(utl::parseRealFunc(cline,pd));
      }
      std::cout << std::endl;
    }
  }

  else if (!strncasecmp(keyWord,"FIXPOINTS",9))
  {
    if (ignoreDirichlet) return true; // Ignore all boundary conditions
    if (!this->createFEMmodel()) return false;

    int nfix = atoi(keyWord+9);
    std::cout <<"\nNumber of fixed points: "<< nfix << std::endl;
    for (int i = 0; i < nfix && (cline = utl::readLine(is)); i++)
    {
      int patch = atoi(strtok(cline," "));
      double rx = atof(strtok(NULL," "));
      int bcode = (cline = strtok(NULL," ")) ? atoi(cline) : 123;
      if (patch < 1 || patch > (int)myModel.size())
      {
	std::cerr <<" *** SIM1D::parse: Invalid patch index "
		  << patch << std::endl;
	return false;
      }
      std::cout <<"\tConstraining P"<< patch
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
                           int&)
{
  if (patch < 1 || patch > (int)myModel.size())
    return constrError("patch index ",patch);

  std::cout <<"\tConstraining P"<< patch;
  if (ldim == 0)
    std::cout << " V"<< lndx;
  std::cout <<" in direction(s) "<< dirs;
  if (code) std::cout <<" code = "<< code <<" ";
#if SP_DEBUG > 1
  std::cout << std::endl;
#endif

  ASMs1D* pch = static_cast<ASMs1D*>(myModel[patch-1]);
  if (abs(ldim) > 0)
    pch->constrainPatch(dirs,code);
  else switch (lndx) // Vertex constraints
    {
    case 1: pch->constrainNode(0.0,dirs,code); break;
    case 2: pch->constrainNode(1.0,dirs,code); break;
    default: std::cout << std::endl;
      return constrError("vertex index ",lndx);
    }

  return true;
}


ASMbase* SIM1D::readPatch (std::istream& isp, int pchInd) const
{
  ASMbase* pch = ASM1D::create(opt.discretization,nd,nf);
  if (pch)
  {
    if (!pch->read(isp))
      delete pch, pch = NULL;
    else if (pch->empty() || this->getLocalPatchIndex(pchInd+1) < 1)
      delete pch, pch = NULL;
    else
      pch->idx = myModel.size();
  }

  return pch;
}


bool SIM1D::readPatches (std::istream& isp, PatchVec& patches,
                         const char* whiteSpace)
{
  ASMbase* pch = NULL;
  for (int pchInd = 1; isp.good(); pchInd++)
    if ((pch = ASM1D::create(opt.discretization,nd,nf)))
    {
      std::cout << whiteSpace <<"Reading patch "<< pchInd << std::endl;
      if (!pch->read(isp))
      {
        delete pch;
        return false;
      }
      else if (pch->empty() || this->getLocalPatchIndex(pchInd+1) < 1)
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
    else
      ok = static_cast<ASMs1D*>(myModel[i])->generateOrientedFEModel(XZp);

  if (nGlPatches == 0 && !adm.isParallel())
    nGlPatches = myModel.size();

  return ok;
}


ASMbase* SIM1D::createDefaultGeometry () const
{
  std::istringstream unitLine("100 1 0 0\n1 0\n2 2\n0 0 1 1\n0\n1\n");

  return this->readPatch(unitLine,1);
}


Vector SIM1D::getSolution (const Vector& psol, double u,
                           int deriv, int patch) const
{
  return this->SIMgeneric::getSolution(psol,&u,deriv,patch);
}


bool SIM1D::updateRotations (const Vector& incSol, double alpha)
{
  if (nf != 6) return true;
  if (incSol.empty()) return false;

  bool ok = true;
  Vector locSol;
  for (size_t i = 0; i < myModel.size() && ok; i++)
  {
    myModel[i]->extractNodeVec(incSol,locSol);
    if (alpha != 0.0 && alpha != 1.0) locSol *= alpha;
    ok = static_cast<ASMs1D*>(myModel[i])->updateRotations(locSol,alpha==0.0);
  }

  return ok;
}
