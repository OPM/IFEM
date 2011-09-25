// $Id$
//==============================================================================
//!
//! \file SIM3D.C
//!
//! \date Dec 08 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for 3D NURBS-based FEM analysis.
//!
//==============================================================================

#include "SIM3D.h"
#include "ASMs3Dmx.h"
#include "ASMs3DmxLag.h"
#include "ASMs3DSpec.h"
#include "Functions.h"
#include "Utilities.h"
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>


SIM3D::SIM3D (bool checkRHS, unsigned char n1, unsigned char n2)
{
  nf[0] = n1;
  nf[1] = n2;
  checkRHSys = checkRHS;
}


bool SIM3D::parse (char* keyWord, std::istream& is)
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
	  pch = new ASMs3DmxLag(cline,checkRHSys,nf[0],nf[1]);
	else
	  pch = new ASMs3DLag(cline,checkRHSys,nf[0]);
	break;
      case Spectral:
	pch = new ASMs3DSpec(cline,checkRHSys,nf[0]);
	break;
      default:
	if (nf[1] > 0)
	  pch = new ASMs3Dmx(cline,checkRHSys,nf[0],nf[1]);
	else
	  pch = new ASMs3D(cline,checkRHSys,nf[0]);
      }
      if (pch->empty() || this->getLocalPatchIndex(i+1) < 1)
	delete pch;
      else
	myModel.push_back(pch);
    }

    if ((int)myModel.size() < npatch)
    {
      std::cerr <<" *** SIM3D::parse: Expected "<< npatch
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
    this->readPatches(isp);

    if (myModel.empty())
    {
      std::cerr <<" *** SIM3D::parse: No patches read"<< std::endl;
      return false;
    }
  }

  else if (!strncasecmp(keyWord,"NODEFILE",8))
  {
    if (!this->createFEMmodel()) return false;

    bool oneBasedIdx = keyWord[8] == '1';
    size_t i = (oneBasedIdx || keyWord[8] == '0') ? 9 : 8;
    while (i < strlen(keyWord) && isspace(keyWord[i])) i++;
    std::ifstream isn(keyWord+i);
    if (isn)
      std::cout <<"\nReading data file "<< keyWord+i << std::endl;
    else
    {
      std::cerr <<" *** SIM3D::read: Failure opening input file "
		<< std::string(keyWord+i) << std::endl;
      return false;
    }

    while (isn.good())
    {
      int patch = 0;
      isn >> patch;
      if (!oneBasedIdx) ++patch;
      int pid = this->getLocalPatchIndex(patch);
      if (pid < 0) return false;

      ASMs3D::BlockNodes n;
      for (i = 0; i <  8 && isn.good(); i++)
	isn >> n.ibnod[i];
      for (i = 0; i < 12 && isn.good(); i++)
	isn >> n.edges[i].icnod >> n.edges[i].incr;
      for (i = 0; i <  6 && isn.good(); i++)
	isn >> n.faces[i].isnod >> n.faces[i].incrI >> n.faces[i].incrJ;
      isn >> n.iinod;

      if (!oneBasedIdx)
      {
	// We always require the node numbers to be 1-based
	for (i = 0; i <  8; i++) ++n.ibnod[i];
	for (i = 0; i < 12; i++) ++n.edges[i].icnod;
	for (i = 0; i <  6; i++) ++n.faces[i].isnod;
	++n.iinod;
      }

      if (isn.good() && pid > 0)
	if (!static_cast<ASMs3D*>(myModel[pid-1])->assignNodeNumbers(n))
	{
	  std::cerr <<" *** SIM3D::parse: Failed to assign node numbers"
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
    std::ifstream isp(keyWord+i);
    if (isp)
      std::cout <<"\nReading data file "<< keyWord+i << std::endl;
    else
    {
      std::cerr <<" *** SIM3D::read: Failure opening input file "
		<< std::string(keyWord+i) << std::endl;
      return false;
    }

    while (isp.good())
    {
      Property p;
      int ldim, lindx = 0;
      isp >> p.pindx >> p.patch >> ldim;
      if (ldim < 3) isp >> lindx;

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

	cline = strtok(NULL," ");
	myScalars[code] = const_cast<RealFunc*>(utl::parseRealFunc(cline,d));
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
	std::cerr <<" *** SIM3D::parse: Invalid patch index "
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
	int addw = atoi(strtok(NULL," "));
	for (int j = ipatch; j < patch; j++)
	{
	  std::cout <<"\tRefining P"<< j+1
		    <<" "<< addu <<" "<< addv <<" "<< addw << std::endl;
	  ASMs3D* pch = static_cast<ASMs3D*>(myModel[j]);
	  pch->uniformRefine(0,addu);
	  pch->uniformRefine(1,addv);
	  pch->uniformRefine(2,addw);
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
	  static_cast<ASMs3D*>(myModel[j])->refine(dir-1,xi);
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
      int addw  = atoi(strtok(NULL," "));
      if (patch == 0 || abs(patch) > (int)myModel.size())
      {
	std::cerr <<" *** SIM3D::parse: Invalid patch index "
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
		  <<" "<< addu <<" "<< addv <<" "<< addw << std::endl;
	static_cast<ASMs3D*>(myModel[j])->raiseOrder(addu,addv,addw);
      }
    }
  }

  else if (!strncasecmp(keyWord,"TOPOLOGYFILE",12))
  {
    if (!this->createFEMmodel()) return false;

    size_t i = 12; while (i < strlen(keyWord) && isspace(keyWord[i])) i++;
    std::ifstream ist(keyWord+i);
    if (ist)
      std::cout <<"\nReading data file "<< keyWord+i << std::endl;
    else
    {
      std::cerr <<" *** SIM3D::read: Failure opening input file "
		<< std::string(keyWord+i) << std::endl;
      return false;
    }

    while ((cline = utl::readLine(ist)))
    {
      int master = atoi(strtok(cline," "))+1;
      int mFace  = atoi(strtok(NULL," "))+1;
      int slave  = atoi(strtok(NULL," "))+1;
      int sFace  = atoi(strtok(NULL," "))+1;
      int swapd  = atoi(strtok(NULL," "));
      int rev_u  = atoi(strtok(NULL," "));
      int rev_v  = atoi(strtok(NULL," "));
      int orient = 4*swapd+2*rev_u+rev_v;
      if (master == slave ||
	  master < 1 || master > (int)myModel.size() ||
	  slave  < 1 || slave  > (int)myModel.size())
      {
	std::cerr <<" *** SIM3D::parse: Invalid patch indices "
		  << master <<" "<< slave << std::endl;
	return false;
      }
      std::cout <<"\tConnecting P"<< slave <<" F"<< sFace
		<<" to P"<< master <<" F"<< mFace
		<<" orient "<< orient << std::endl;
      ASMs3D* spch = static_cast<ASMs3D*>(myModel[slave-1]);
      ASMs3D* mpch = static_cast<ASMs3D*>(myModel[master-1]);
      if (!spch->connectPatch(sFace,*mpch,mFace,orient))
	return false;
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
      int mFace  = atoi(strtok(NULL," "));
      int slave  = atoi(strtok(NULL," "));
      int sFace  = atoi(strtok(NULL," "));
      int orient = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
      if (master == slave ||
	  master < 1 || master > (int)myModel.size() ||
	  slave  < 1 || slave  > (int)myModel.size())
      {
	std::cerr <<" *** SIM3D::parse: Invalid patch indices "
		  << master <<" "<< slave << std::endl;
	return false;
      }
      std::cout <<"\tConnecting P"<< slave <<" F"<< sFace
		<<" to P"<< master <<" F"<< mFace
		<<" orient "<< orient << std::endl;
      ASMs3D* spch = static_cast<ASMs3D*>(myModel[slave-1]);
      ASMs3D* mpch = static_cast<ASMs3D*>(myModel[master-1]);
      if (!spch->connectPatch(sFace,*mpch,mFace,orient))
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
      int pfdir = atoi(strtok(NULL," "));
      if (patch < 1 || patch > (int)myModel.size())
      {
	std::cerr <<" *** SIM3D::parse: Invalid patch index "
		  << patch << std::endl;
	return false;
      }
      std::cout <<"\tPeriodic "<< char('H'+pfdir) <<"-direction P"<< patch
		<< std::endl;
      static_cast<ASMs3D*>(myModel[patch-1])->closeFaces(pfdir);
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
      int pface = atoi(strtok(NULL," "));
      int bcode = atoi(strtok(NULL," "));
      double pd = (cline = strtok(NULL," ")) ? atof(cline) : 0.0;

      patch = this->getLocalPatchIndex(patch);
      if (patch < 1) continue;

      if (pface > 10)
      {
	if (!this->addConstraint(patch,pface%10,pface/10,pd,bcode))
	  return false;
      }
      else if (pd == 0.0)
      {
	int ldim = pface < 0 ? 0 : 2;
	if (!this->addConstraint(patch,abs(pface),ldim,bcode%1000))
	  return false;
      }
      else
      {
	int ldim = pface < 0 ? 0 : 2;
	int code = 1000 + bcode;
	while (myScalars.find(code) != myScalars.end())
	  code += 1000;

	if (!this->addConstraint(patch,abs(pface),ldim,bcode%1000,code))
	  return false;

	cline = strtok(NULL," ");
	myScalars[code] = const_cast<RealFunc*>(utl::parseRealFunc(cline,pd));
      }
      if (pface < 10) std::cout << std::endl;
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
      double rz = atof(strtok(NULL," "));
      int bcode = (cline = strtok(NULL," ")) ? atoi(cline) : 123;

      int pid = this->getLocalPatchIndex(patch);
      if (pid < 1) continue;

      std::cout <<"\tConstraining P"<< patch
		<<" point at "<< rx <<" "<< ry <<" "<< rz
		<<" with code "<< bcode << std::endl;
      static_cast<ASMs3D*>(myModel[pid-1])->constrainNode(rx,ry,rz,bcode);
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
  std::cerr <<" *** SIM3D::addConstraint: Invalid "<< lab << idx << std::endl;
  return false;
}


bool SIM3D::addConstraint (int patch, int lndx, int ldim, int dirs, int code)
{
  if (patch < 1 || patch > (int)myModel.size())
    return constrError("patch index ",patch);

  std::cout <<"\tConstraining P"<< patch
	    << (ldim == 0 ? " V" : (ldim == 1 ? " E": " F")) << lndx
	    <<" in direction(s) "<< dirs;
  if (code) std::cout <<" code = "<< code <<" ";

  ASMs3D* pch = static_cast<ASMs3D*>(myModel[patch-1]);
  switch (ldim)
    {
    case 0: // Vertex constraints
      switch (lndx)
	{
	case 1: pch->constrainCorner(-1,-1,-1,dirs,code); break;
	case 2: pch->constrainCorner( 1,-1,-1,dirs,code); break;
	case 3: pch->constrainCorner(-1, 1,-1,dirs,code); break;
	case 4: pch->constrainCorner( 1, 1,-1,dirs,code); break;
	case 5: pch->constrainCorner(-1,-1, 1,dirs,code); break;
	case 6: pch->constrainCorner( 1,-1, 1,dirs,code); break;
	case 7: pch->constrainCorner(-1, 1, 1,dirs,code); break;
	case 8: pch->constrainCorner( 1, 1, 1,dirs,code); break;
	default: std::cout << std::endl;
	  return constrError("vertex index ",lndx);
	}
      break;

    case 1: // Edge constraints
      if (lndx > 0 && lndx <= 12)
	pch->constrainEdge(lndx,dirs,code);
      else
      {
	std::cout << std::endl;
	return constrError("edge index ",lndx);
      }
      break;

    case 2:
      switch (lndx)
	{
	case 1: pch->constrainFace(-1,dirs,code); break;
	case 2: pch->constrainFace( 1,dirs,code); break;
	case 3: pch->constrainFace(-2,dirs,code); break;
	case 4: pch->constrainFace( 2,dirs,code); break;
	case 5: pch->constrainFace(-3,dirs,code); break;
	case 6: pch->constrainFace( 3,dirs,code); break;
	default: std::cout << std::endl;
	  return constrError("face index ",lndx);
	}
      break;

    default:
      std::cout << std::endl;
      return constrError("local dimension switch ",ldim);
    }

  return true;
}


bool SIM3D::addConstraint (int patch, int lndx, int line, double xi, int dirs)
{
  if (patch < 1 || patch > (int)myModel.size())
    return constrError("patch index ",patch);

  std::cout <<"\tConstraining P"<< patch
	    <<" F"<< lndx <<" L"<< line <<" at xi="<< xi
	    <<" in direction(s) "<< dirs << std::endl;

  ASMs3D* pch = static_cast<ASMs3D*>(myModel[patch-1]);
  switch (line)
    {
    case 1: // Face line constraints in local I-direction
      switch (lndx)
	{
	case 1: pch->constrainLine(-1,2,xi,dirs); break;
	case 2: pch->constrainLine( 1,2,xi,dirs); break;
	case 3: pch->constrainLine(-2,3,xi,dirs); break;
	case 4: pch->constrainLine( 2,3,xi,dirs); break;
	case 5: pch->constrainLine(-3,1,xi,dirs); break;
	case 6: pch->constrainLine( 3,1,xi,dirs); break;
	default: return constrError("face index ",lndx);
	}
      break;

    case 2: // Face line constraints in local J-direction
      switch (lndx)
	{
	case 1: pch->constrainLine(-1,3,xi,dirs); break;
	case 2: pch->constrainLine( 1,3,xi,dirs); break;
	case 3: pch->constrainLine(-2,1,xi,dirs); break;
	case 4: pch->constrainLine( 2,1,xi,dirs); break;
	case 5: pch->constrainLine(-3,2,xi,dirs); break;
	case 6: pch->constrainLine( 3,2,xi,dirs); break;
	default: return constrError("face index ",lndx);
	}
      break;

    default:
      return constrError("face line index ",line);
    }

  return true;
}


void SIM3D::setQuadratureRule (size_t ng)
{
  for (size_t i = 0; i < myModel.size(); i++)
    if (!myModel.empty())
      static_cast<ASMs3D*>(myModel[i])->setGauss(ng);
}


void SIM3D::readPatches (std::istream& isp)
{
  ASMbase* pch = 0;
  for (int patchNo = 1; isp.good(); patchNo++)
  {
    std::cout <<"Reading patch "<< patchNo << std::endl;
    switch (discretization) {
      case Lagrange:
        if (nf[1] > 0)
          pch = new ASMs3DmxLag(isp,checkRHSys,nf[0],nf[1]);
        else
          pch = new ASMs3DLag(isp,checkRHSys,nf[0]);
        break;
      case Spectral:
        pch = new ASMs3DSpec(isp,checkRHSys,nf[0]);
        break;
      default:
        if (nf[1] > 0)
          pch = new ASMs3Dmx(isp,checkRHSys,nf[0],nf[1]);
        else
          pch = new ASMs3D(isp,checkRHSys,nf[0]);
    }
    if (pch->empty() || this->getLocalPatchIndex(patchNo) < 1)
      delete pch;
    else
      myModel.push_back(pch);
  }
}
