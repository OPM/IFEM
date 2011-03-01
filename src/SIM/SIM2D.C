// $Id$
//==============================================================================
//!
//! \file SIM2D.C
//!
//! \date Feb 04 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for 2D NURBS-based FEM analysis.
//!
//==============================================================================

#include "SIM2D.h"
#include "ASMs2Dmx.h"
#include "ASMs2DmxLag.h"
#include "ASMs2DSpec.h"
#include "Functions.h"
#include "Utilities.h"
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>


SIM2D::SIM2D (unsigned char n1, unsigned char n2)
{
  nf[0] = n1;
  nf[1] = n2;
}


bool SIM2D::parse (char* keyWord, std::istream& is)
{
  char* cline = 0;
  if (!strncasecmp(keyWord,"PATCHES",7))
  {
    ASMbase* pch = 0;
    int npatch = atoi(keyWord+7);
    std::cout <<"\nNumber of patches: "<< npatch << std::endl;
    for (int i = 0; i < npatch && (cline = utl::readLine(is)); i++)
    {
      cline = strtok(cline," ");
      switch (discretization) {
      case Lagrange:
	if (nf[1] > 0)
	  pch = new ASMs2DmxLag(cline,2,nf[0],nf[1]);
	else
	  pch = new ASMs2DLag(cline,2,nf[0]);
	break;
      case Spectral:
	pch = new ASMs2DSpec(cline,2,nf[0]);
	break;
      default:
	if (nf[1] > 0)
	  pch = new ASMs2Dmx(cline,2,nf[0],nf[1]);
	else
	  pch = new ASMs2D(cline,2,nf[0]);
      }
      if (pch->empty() || this->getLocalPatchIndex(i+1) < 1)
	delete pch;
      else
	myModel.push_back(pch);
    }

    if ((int)myModel.size() < npatch)
    {
      std::cerr <<" *** SIM2D::parse: Expected "<< npatch
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
	if (nf[1] > 0)
	  pch = new ASMs2DmxLag(isp,2,nf[0],nf[1]);
	else
	  pch = new ASMs2DLag(isp,2,nf[0]);
	break;
      case Spectral:
	pch = new ASMs2DSpec(isp,2,nf[0]);
	break;
      default:
	if (nf[1] > 0)
	  pch = new ASMs2Dmx(isp,2,nf[0],nf[1]);
	else
	  pch = new ASMs2D(isp,2,nf[0]);
      }
      if (pch->empty() || this->getLocalPatchIndex(patchNo) < 1)
	delete pch;
      else
	myModel.push_back(pch);
    }

    if (myModel.empty())
    {
      std::cerr <<" *** SIM2D::parse: No patches read"<< std::endl;
      return false;
    }
  }

  else if (!strncasecmp(keyWord,"NODEFILE",8))
  {
    if (!this->createFEMmodel()) return false;

    bool oneBasedIdx = keyWord[8] == '1';
    size_t i = (oneBasedIdx || keyWord[8] == '0') ? 9 : 8;
    while (i < strlen(keyWord) && isspace(keyWord[i])) i++;
    std::cout <<"\nReading data file "<< keyWord+i << std::endl;
    std::ifstream isn(keyWord+i);
    while (isn.good())
    {
      int patch = 0;
      isn >> patch;
      if (!oneBasedIdx) ++patch;
      int pid = this->getLocalPatchIndex(patch);
      if (pid < 0) return false;

      ASMs2D::BlockNodes n;
      for (i = 0; i < 4 && isn.good(); i++)
	isn >> n.ibnod[i];
      for (i = 0; i < 4 && isn.good(); i++)
	isn >> n.edges[i].icnod >> n.edges[i].incr;
      isn >> n.iinod;

      if (!oneBasedIdx)
      {
	// We always require the node numbers to be 1-based
	for (i = 0; i < 4; i++) ++n.ibnod[i];
	for (i = 0; i < 4; i++) ++n.edges[i].icnod;
	++n.iinod;
      }

      if (isn.good() && pid > 0)
	if (!static_cast<ASMs2D*>(myModel[pid-1])->assignNodeNumbers(n))
	{
	  std::cerr <<" *** SIM2D::parse: Failed to assign node numbers"
		    <<" for patch "<< patch << std::endl;
	  return false;
	}
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
      if (ldim < 2) isp >> lindx;

      if (!oneBasedIdx)
      {
	// We always require the item indices to be 1-based
	++p.patch;
	++lindx;
      }

      p.ldim = ldim;
      p.lindx = lindx;
      p.patch = this->getLocalPatchIndex(p.patch);
      if (p.patch > 0 && isp.good())
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
  // properties as well as grid refinement without using the GPM module.

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
	std::cerr <<" *** SIM2D::parse: Invalid patch index "
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
	int addv = atoi(strtok(NULL," "));
	for (int j = ipatch; j < patch; j++)
	{
	  std::cout <<"\tRefining P"<< j+1
		    <<" "<< addu <<" "<< addv << std::endl;
	  ASMs2D* pch = static_cast<ASMs2D*>(myModel[j]);
	  pch->uniformRefine(0,addu);
	  pch->uniformRefine(1,addv);
	}
      }
      else
      {
	int dir = atoi(strtok(NULL," "));
	RealArray xi;
	while ((cline = strtok(NULL," ")))
	  xi.push_back(atof(cline));
	for (int j = ipatch; j < patch; j++)
	{
	  std::cout <<"\tRefining P"<< j+1 <<" dir="<< dir;
	  for (size_t i = 0; i < xi.size(); i++)
	    std::cout <<" "<< xi[i];
	  std::cout << std::endl;
	  static_cast<ASMs2D*>(myModel[j])->refine(dir-1,xi);
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
      int addv  = atoi(strtok(NULL," "));
      if (patch == 0 || abs(patch) > (int)myModel.size())
      {
	std::cerr <<" *** SIM2D::parse: Invalid patch index "
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
		  <<" "<< addu <<" "<< addv << std::endl;
	static_cast<ASMs2D*>(myModel[j])->raiseOrder(addu,addv);
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
      int mEdge  = atoi(strtok(NULL," "));
      int slave  = atoi(strtok(NULL," "));
      int sEdge  = atoi(strtok(NULL," "));
      bool rever = (cline = strtok(NULL," ")) ? cline[0] == 'R' : false;
      if (master == slave ||
	  master < 1 || master > (int)myModel.size() ||
	  slave  < 1 || slave  > (int)myModel.size())
      {
	std::cerr <<" *** SIM2D::parse: Invalid patch indices "
		  << master <<" "<< slave << std::endl;
	return false;
      }
      std::cout <<"\tConnecting P"<< slave <<" E"<< sEdge
		<<" to P"<< master <<" E"<< mEdge
		<<" reversed? "<< rever << std::endl;
      ASMs2D* spch = static_cast<ASMs2D*>(myModel[slave-1]);
      ASMs2D* mpch = static_cast<ASMs2D*>(myModel[master-1]);
      if (!spch->connectPatch(sEdge,*mpch,mEdge,rever))
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
      int pedir = atoi(strtok(NULL," "));
      if (patch < 1 || patch > (int)myModel.size())
      {
	std::cerr <<" *** SIM2D::parse: Invalid patch index "
		  << patch << std::endl;
	return false;
      }
      std::cout <<"\tPeriodic "<< char('H'+pedir) <<"-direction P"<< patch
		<< std::endl;
      static_cast<ASMs2D*>(myModel[patch-1])->closeEdges(pedir);
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
      int pedge = atoi(strtok(NULL," "));
      int bcode = atoi(strtok(NULL," "));
      double pd = (cline = strtok(NULL," ")) ? atof(cline) : 0.0;

      patch = this->getLocalPatchIndex(patch);
      if (patch < 1) continue;

      int ldim  = pedge >= 0 ? 1 : 0;
      if (pd == 0.0)
      {
	if (!this->addConstraint(patch,abs(pedge),ldim,bcode%1000))
	  return false;
      }
      else
      {
	int code = 1000 + bcode;
	while (myScalars.find(code) != myScalars.end())
	  code += 1000;

	if (!this->addConstraint(patch,abs(pedge),ldim,bcode%1000,code))
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
      double ry = atof(strtok(NULL," "));
      int bcode = (cline = strtok(NULL," ")) ? atoi(cline) : 123;

      int pid = this->getLocalPatchIndex(patch);
      if (pid < 1) continue;

      std::cout <<"\tConstraining P"<< patch
		<<" point at "<< rx <<" "<< ry
		<<" with code "<< bcode << std::endl;
      static_cast<ASMs2D*>(myModel[pid-1])->constrainNode(rx,ry,bcode);
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
  std::cerr <<" *** SIM2D::addConstraint: Invalid "<< lab << idx << std::endl;
  return false;
}


bool SIM2D::addConstraint (int patch, int lndx, int ldim, int dirs, int code)
{
  if (patch < 1 || patch > (int)myModel.size())
    return constrError("patch index ",patch);

  std::cout <<"\tConstraining P"<< patch
	    << (ldim == 0 ? " V" : " E") << lndx <<" in direction(s) "<< dirs;
  if (code) std::cout <<" code = "<< code <<" ";

  ASMs2D* pch = static_cast<ASMs2D*>(myModel[patch-1]);
  switch (ldim)
    {
    case 0: // Vertex constraints
      switch (lndx)
	{
	case 1: pch->constrainCorner(-1,-1,dirs,code); break;
	case 2: pch->constrainCorner( 1,-1,dirs,code); break;
	case 3: pch->constrainCorner(-1, 1,dirs,code); break;
	case 4: pch->constrainCorner( 1, 1,dirs,code); break;
	default: std::cout << std::endl;
	  return constrError("vertex index ",lndx);
	}
      break;

    case 1: // Edge constraints
      switch (lndx)
	{
	case 1: pch->constrainEdge(-1,dirs,code); break;
	case 2: pch->constrainEdge( 1,dirs,code); break;
	case 3: pch->constrainEdge(-2,dirs,code); break;
	case 4: pch->constrainEdge( 2,dirs,code); break;
	default: std::cout << std::endl;
	  return constrError("edge index ",lndx);
	}
      break;

    default:
      std::cout << std::endl;
      return constrError("local dimension switch ",ldim);
    }

  return true;
}


void SIM2D::setQuadratureRule (size_t ng)
{
  for (size_t i = 0; i < myModel.size(); i++)
    if (!myModel.empty())
      static_cast<ASMs2D*>(myModel[i])->setGauss(ng);
}
