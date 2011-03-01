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
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>


SIM1D::SIM1D (unsigned char n_f)
{
  nf = n_f;
}


bool SIM1D::parse (char* keyWord, std::istream& is)
{
  char* cline = 0;
  if (!strncasecmp(keyWord,"PATCHES",7))
  {
    ASMbase* pch = 0;
    int npatch = atoi(keyWord+7);
    std::cout <<"\nNumber of patches: "<< npatch << std::endl;
    for (int i = 0; i < npatch && (cline = utl::readLine(is)); i++)
    {
      switch (discretization) {
      case Lagrange:
	pch = new ASMs1DLag(strtok(cline," "),1,nf);
	break;
      case Spectral:
	pch = new ASMs1DSpec(strtok(cline," "),1,nf);
	break;
      default:
	pch = new ASMs1D(strtok(cline," "),1,nf);
      }
      if (pch->empty())
	delete pch;
      else
	myModel.push_back(pch);
    }

    if ((int)myModel.size() < npatch)
    {
      std::cerr <<" *** SIM1D::parse: Expected "<< npatch
		<<" patches but could read only "<< myModel.size()
		<< std::endl;
      return false;
    }
  }

  else if (!strncasecmp(keyWord,"PATCHFILE",9))
  {
    size_t i = 9; while (i < strlen(keyWord) && isspace(keyWord[i])) i++;
    std::cout <<"\nReading data file "<< keyWord+i << std::endl;
    std::ifstream isp(keyWord+i);

    ASMbase* pch = 0;
    for (int patchNo = 1; isp.good(); patchNo++)
    {
      std::cout <<"Reading patch "<< patchNo << std::endl;
      switch (discretization) {
      case Lagrange:
	pch = new ASMs1DLag(isp,1,nf);
	break;
      case Spectral:
	pch = new ASMs1DSpec(isp,1,nf);
	break;
      default:
	pch = new ASMs1D(isp,1,nf);
      }
      if (pch->empty())
	delete pch;
      else
	myModel.push_back(pch);
    }

    if (myModel.empty())
    {
      std::cerr <<" *** SIM1D::parse: No patches read"<< std::endl;
      return false;
    }
  }

  else if (!strncasecmp(keyWord,"PROPERTYFILE",12))
  {
    bool oneBasedIdx = keyWord[12] == '1';
    size_t i = (oneBasedIdx || keyWord[12] == '0') ? 13 : 12;
    while (i < strlen(keyWord) && isspace(keyWord[i])) i++;
    std::cout <<"\nReading data file "<< keyWord+i << std::endl;
    std::ifstream isp(keyWord+i);
    while (isp.good())
    {
      Property p;
      int ldim, lindx = 0;
      isp >> p.pindx >> p.patch >> ldim;
      if (ldim < 1) isp >> lindx;

      if (!oneBasedIdx)
      {
	// We always require the item indices to be 1-based
 	++p.patch;
 	++lindx;
      }

      p.ldim = ldim;
      p.lindx = lindx;
      if (isp.good())
	myProps.push_back(p);
    }
  }

  else if (!strncasecmp(keyWord,"DIRICHLET",9))
  {
    if (ignoreDirichlet) return true; // Ignore all boundary conditions

    int ndir = atoi(keyWord+9);
    std::cout <<"\nNumber of Dirichlet properties: "<< ndir << std::endl;
    for (int i = 0; i < ndir && (cline = utl::readLine(is)); i++)
    {
      int code = atoi(strtok(cline," "));
      double d = (cline = strtok(NULL," ")) ? atof(cline) : 0.0;
      std::cout <<"\tDirichlet code "<< code <<": ";
      if (d == 0.0)
      {
	this->setPropertyType(code,Property::DIRICHLET);
	std::cout <<"(fixed)";
      }
      else
      {
	this->setPropertyType(code,Property::DIRICHLET_INHOM);
	if ((cline = strtok(NULL," ")))
	  myScalars[code] = const_cast<RealFunc*>(utl::parseRealFunc(cline,d));
	else
	{
	  std::cout << d;
	  myScalars[code] = new ConstFunc(d);
	}
      }
      std::cout << std::endl;
    }
  }

  // The remaining keywords are retained for backward compatibility with the
  // prototype version. They enable direct specification of topology and
  // properties as well as uniform refinement without using the GPM module.

  else if (!strncasecmp(keyWord,"REFINE",6))
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
	while ((cline = strtok(NULL," ")))
	  xi.push_back(atof(cline));
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

	if ((cline = strtok(NULL," ")))
	  myScalars[code] = const_cast<RealFunc*>(utl::parseRealFunc(cline,pd));
	else
	  myScalars[code] = new ConstFunc(pd);
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
