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
#include "ASMs1DLag.h"
#include "ASMs1DSpec.h"
#include "Functions.h"
#include "Utilities.h"
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "tinyxml.h"


SIM1D::SIM1D (unsigned char n_f)
{
  nf = n_f;
}


bool SIM1D::parse (char* keyWord, std::istream& is)
{
  char* cline = 0;
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
      std::cout <<"\tPeriodic P" << patch << std::endl;
      static_cast<ASMs1D*>(myModel[patch-1])->closeEnds();
    }
  }

  else if (!strncasecmp(keyWord,"CONSTRAINTS",11))
  {
    if (ignoreDirichlet) return true; // Ignore all boundary conditions
    if (!this->createFEMmodel()) return false;

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
	if (!this->addConstraint(patch,pvert,0,bcode%1000))
	  return false;
      }
      else
      {
	int code = 1000 + bcode;
	while (myScalars.find(code) != myScalars.end())
	  code += 1000;

	if (!this->addConstraint(patch,pvert,0,bcode%1000,code))
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
    return this->SIMbase::parse(keyWord,is);

  return true;
}


bool SIM1D::parseGeometryTag(const TiXmlElement* elem)
{
  // The remaining keywords are retained for backward compatibility with the
  // prototype version. They enable direct specification of topology and
  // properties as well as grid refinement without using the GPM module.

  if (!strcasecmp(elem->Value(),"refine")) {
    ASMs1D* pch = 0;
    bool uniform=true;
    if (elem->Attribute("type") && !strcasecmp(elem->Attribute("type"),"nonuniform"))
      uniform = false;
    size_t lowpatch = 1, uppatch = 2;

    if (elem->Attribute("patch"))
      lowpatch = uppatch = atoi(elem->Attribute("patch"));
    if (elem->Attribute("lowerpatch")) {
      lowpatch = atoi(elem->Attribute("lowerpatch"));
      uppatch = myModel.size();
    }
    if (elem->Attribute("upperpatch"))
      uppatch = atoi(elem->Attribute("upperpatch"));

    if (lowpatch < 1 || uppatch > myModel.size())
    {
      std::cerr <<" *** SIM1D::parse: Invalid patch index "
        << "lower: " << lowpatch << " upper: " << uppatch << std::endl;
      return false;
    }

    if (uniform)
    {
      int addu=0;
      if (elem->Attribute("u"))
        addu = atoi(elem->Attribute("u"));
      for (size_t j = lowpatch-1; j < uppatch; j++) {
        if ((pch = dynamic_cast<ASMs1D*>(myModel[j]))) {
          std::cout <<"\tRefining P"<< j+1
            <<" "<< addu << std::endl;
          pch->uniformRefine(addu);
        }
      }
    } else {
      char* cline2 = strdup(elem->FirstChildElement()->Value());
      char* cline = cline2;
      strtok(cline," ");
      RealArray xi;
      while (cline) {
        xi.push_back(atof(cline));
        cline = strtok(NULL," ");
      }
      free(cline2);
      for (size_t j = lowpatch-1; j < uppatch; j++) {
        if ((pch = dynamic_cast<ASMs1D*>(myModel[j])))
        {
          std::cout <<"\tRefining P"<< j+1;
          for (size_t i = 0; i < xi.size(); i++)
            std::cout <<" "<< xi[i];
          std::cout << std::endl;
          pch->refine(xi);
        }
      }
    }
  } else if (!strcasecmp(elem->Value(),"raiseorder")) {
    size_t lowpatch = 1, uppatch = 2;
    if (elem->Attribute("patch"))
      lowpatch = uppatch = atoi(elem->Attribute("patch"));
    if (elem->Attribute("lowerpatch")) {
      lowpatch = atoi(elem->Attribute("lowerpatch"));
      uppatch = myModel.size();
    }
    if (elem->Attribute("upperpatch"))
      uppatch = atoi(elem->Attribute("upperpatch"));

    if (lowpatch < 1 || uppatch > myModel.size())
    {
      std::cerr <<" *** SIM2D::parse: Invalid patch index "
        << "lower: " << lowpatch << " upper: " << uppatch << std::endl;
      return false;
    }
    int addu=0;
    if (elem->Attribute("u"))
      addu = atoi(elem->Attribute("u"));
    ASMs1D* pch;
    for (size_t j = lowpatch-1; j < uppatch; j++) {
      if ((pch = dynamic_cast<ASMs1D*>(myModel[j]))) {
        std::cout <<"\tRaising order of P"<< j+1 <<" "<< addu << std::endl;
        pch->raiseOrder(addu);
      }
    }
  }
  else if (!strcasecmp(elem->Value(),"topology")) {
    if (createFEMmodel()) return false;

    const TiXmlElement* child = elem->FirstChildElement("connection");
    while (child) {
      int master=0, slave=0, mVert=0, sVert=0;
      if (child->Attribute("master"))
        master = atoi(child->Attribute("master"));
      if (child->Attribute("mvert"))
        mVert = atoi(child->Attribute("mvert"));
      if (child->Attribute("slave"))
        slave = atoi(child->Attribute("slave"));
      if (child->Attribute("svert"))
        sVert = atoi(child->Attribute("svert"));

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

      child = child->NextSiblingElement();
    }
  } else if (!strcasecmp(elem->Value(),"periodic")) {
    if (createFEMmodel()) return false;
    int patch=0;
    if (elem->Attribute("patch"))
      patch = atoi(elem->Attribute("patch"));

    if (patch < 1 || patch > (int)myModel.size()) {
      std::cerr <<" *** SIM1D::parse: Invalid patch index "
                << patch << std::endl;
      return false;
    }
    std::cout <<"\tPeriodic P" << patch << std::endl;
    static_cast<ASMs1D*>(myModel[patch-1])->closeEnds();
  }
  // TODO: constraints? fixpoints?

  return true;
}


bool SIM1D::parse(const TiXmlElement* elem)
{
  bool result=SIMbase::parse(elem);

  const TiXmlElement* child = elem->FirstChildElement();
  while (child) {
    if (!strcasecmp(elem->Value(),"geometry"))
      result &= parseGeometryTag(child);
    child = child->NextSiblingElement();
  }
  return result;
}


/*!
  \brief Local-scope convenience function for error message generation.
*/

static bool constrError (const char* lab, int idx)
{
  std::cerr <<" *** SIM1D::addConstraint: Invalid "<< lab << idx << std::endl;
  return false;
}


bool SIM1D::addConstraint (int patch, int lndx, int, int dirs, int code)
{
  if (patch < 1 || patch > (int)myModel.size())
    return constrError("patch index ",patch);

  std::cout <<"\tConstraining P"<< patch
	    << " V" << lndx <<" in direction(s) "<< dirs;
  if (code) std::cout <<" code = "<< code <<" ";

  ASMs1D* pch = static_cast<ASMs1D*>(myModel[patch-1]);
  switch (lndx) // Vertex constraints
    {
    case 1: pch->constrainNode(0.0,code); break;
    case 2: pch->constrainNode(1.0,code); break;
    default: std::cout << std::endl;
      return constrError("vertex index ",lndx);
    }

  return true;
}


void SIM1D::setQuadratureRule (size_t ng)
{
  for (size_t i = 0; i < myModel.size(); i++)
    if (!myModel.empty())
      static_cast<ASMs1D*>(myModel[i])->setGauss(ng);
}


bool SIM1D::readPatch (std::istream& isp, int pchInd)
{
  ASMs1D* pch = 0;
  switch (discretization) {
  case ASM::Lagrange:
    pch = new ASMs1DLag(1,nf);
    break;
  case ASM::Spectral:
    pch = new ASMs1DSpec(1,nf);
    break;
  default:
    pch = new ASMs1D(1,nf);
  }

  if (!pch->read(isp))
  {
    delete pch;
    return false;
  }
  else if (pch->empty())
    delete pch;
  else
    myModel.push_back(pch);

  return true;
}


bool SIM1D::readPatches (std::istream& isp)
{
  ASMbase* pch = 0;
  for (int patchNo = 1; isp.good(); patchNo++)
  {
    std::cout <<"Reading patch "<< patchNo << std::endl;
    switch (discretization)
      {
      case ASM::Lagrange:
        pch = new ASMs1DLag(1,nf);
        break;
      case ASM::Spectral:
        pch = new ASMs1DSpec(1,nf);
        break;
      default:
        pch = new ASMs1D(1,nf);
      }
    if (!pch->read(isp))
    {
      delete pch;
      return false;
    }
    else if (pch->empty())
      delete pch;
    else
      myModel.push_back(pch);
  }

  return true;
}
