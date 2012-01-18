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
#include "ASMs2DC1.h"
#include "Functions.h"
#include "Utilities.h"

#include "tinyxml.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif


/*!
  A struct defining a patch interface for C1-continuous models.
*/

struct Interface
{
  std::pair<ASMs2DC1*,int> master; //!< Patch and edge index of the master
  std::pair<ASMs2DC1*,int> slave;  //!< Patch and edge index of the slave
  bool reversed;                   //!< Relative orientation toggle
  //! \brief Constructor initializing an Interface instance.
  Interface(ASMs2DC1* m, int me, ASMs2DC1* s, int se, bool r = false)
  {
    master = std::make_pair(m,me);
    slave = std::make_pair(s,se);
    reversed = r;
  }
};


SIM2D::SIM2D (unsigned char n1, unsigned char n2) : isRefined(false)
{
  nf[0] = n1;
  nf[1] = n2;

  if (nf[1] > 0) mixedFEM = true;
}


bool SIM2D::parseGeometryTag(const TiXmlElement* elem)
{
  // The remaining keywords are retained for backward compatibility with the
  // prototype version. They enable direct specification of topology and
  // properties as well as grid refinement without using the GPM module.

  if (!strcasecmp(elem->Value(),"refine") && !isRefined) {
    ASM2D* pch = 0;
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
      std::cerr <<" *** SIM2D::parse: Invalid patch index "
        << "lower: " << lowpatch << " upper: " << uppatch << std::endl;
      return false;
    }

    if (uniform)
    {
      int addu=0, addv = 0;
      if (elem->Attribute("u"))
        addu = atoi(elem->Attribute("u"));
      if (elem->Attribute("v"))
        addv = atoi(elem->Attribute("v"));
      for (size_t j = lowpatch-1; j < uppatch; j++) {
        if ((pch = dynamic_cast<ASM2D*>(myModel[j]))) {
          std::cout <<"\tRefining P"<< j+1
            <<" "<< addu <<" "<< addv << std::endl;
          pch->uniformRefine(0,addu);
          pch->uniformRefine(1,addv);
        }
      }
    } else {
      int dir=1;
      if (elem->Attribute("dir"))
        dir = atoi(elem->Attribute("dir"));
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
        if ((pch = dynamic_cast<ASM2D*>(myModel[j])))
        {
          std::cout <<"\tRefining P"<< j+1 <<" dir="<< dir;
          for (size_t i = 0; i < xi.size(); i++)
            std::cout <<" "<< xi[i];
          std::cout << std::endl;
          pch->refine(dir-1,xi);
        }
      }
    }
  } else if (!strcasecmp(elem->Value(),"raiseorder") && !isRefined) {
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
    int addu=0, addv = 0;
    if (elem->Attribute("u"))
      addu = atoi(elem->Attribute("u"));
    if (elem->Attribute("v"))
      addv = atoi(elem->Attribute("v"));
    ASM2D* pch;
    for (size_t j = lowpatch-1; j < uppatch; j++) {
      if ((pch = dynamic_cast<ASM2D*>(myModel[j]))) {
        std::cout <<"\tRaising order of P"<< j+1
          <<" "<< addu <<" "<< addv << std::endl;
        pch->raiseOrder(addu,addv);
      }
    }
  }
  else if (!strcasecmp(elem->Value(),"topology")) {
    if (createFEMmodel()) return false;

    const TiXmlElement* child = elem->FirstChildElement("connection");
    while (child) {
      int master=0, slave=0, mEdge=0, sEdge=0;
      if (child->Attribute("master"))
        master = atoi(child->Attribute("master"));
      if (child->Attribute("medge"))
        mEdge = atoi(child->Attribute("medge"));
      if (child->Attribute("slave"))
        slave = atoi(child->Attribute("slave"));
      if (child->Attribute("sedge"))
        sEdge = atoi(child->Attribute("sedge"));

      bool rever=false;
      if (child->Attribute("reverse") && !strcasecmp(child->Attribute("reverse"),"true"))
        rever = true;

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

      child = child->NextSiblingElement();
    }
  } else if (!strcasecmp(elem->Value(),"periodic")) {
    if (createFEMmodel()) return false;
    int patch=0, pedir=1;
    if (elem->Attribute("patch"))
      patch = atoi(elem->Attribute("patch"));
    if (elem->Attribute("dir"))
      pedir = atoi(elem->Attribute("dir"));

    if (patch < 1 || patch > (int)myModel.size()) {
      std::cerr <<" *** SIM2D::parse: Invalid patch index "
                << patch << std::endl;
      return false;
    }
    std::cout <<"\tPeriodic "<< char('H'+pedir) <<"-direction P"<< patch
              << std::endl;
    static_cast<ASMs2D*>(myModel[patch-1])->closeEdges(pedir);
    // cannot do multi-threaded assembly with periodicities
#ifdef USE_OPENMP
    omp_set_num_threads(1);
#endif
  }
  // TODO: constraints? fixpoints?

  return true;
}


bool SIM2D::parse(const TiXmlElement* elem)
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


bool SIM2D::parse (char* keyWord, std::istream& is)
{
  char* cline = 0;
  if (!strncasecmp(keyWord,"REFINE",6))
  {
    int nref = atoi(keyWord+6);
    if (isRefined) // just read through the next lines without doing anything
      for (int i = 0; i < nref && utl::readLine(is); i++);
    else
    {
      ASM2D* pch = 0;
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
	    if ((pch = dynamic_cast<ASM2D*>(myModel[j])))
	    {
	      std::cout <<"\tRefining P"<< j+1
			<<" "<< addu <<" "<< addv << std::endl;
	      pch->uniformRefine(0,addu);
	      pch->uniformRefine(1,addv);
	    }
	}
	else
	{
	  RealArray xi;
	  int dir = atoi(strtok(NULL," "));
	  if (utl::parseKnots(xi))
	    for (int j = ipatch; j < patch; j++)
	      if ((pch = dynamic_cast<ASM2D*>(myModel[j])))
	      {
		std::cout <<"\tRefining P"<< j+1 <<" dir="<< dir;
		for (size_t i = 0; i < xi.size(); i++)
		  std::cout <<" "<< xi[i];
		std::cout << std::endl;
		pch->refine(dir-1,xi);
	      }
	}
      }
    }
  }

  else if (!strncasecmp(keyWord,"RAISEORDER",10))
  {
    int nref = atoi(keyWord+10);
    if (isRefined) // just read through the next lines without doing anything
      for (int i = 0; i < nref && utl::readLine(is); i++);
    else
    {
      ASM2D* pch = 0;
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
	  if ((pch = dynamic_cast<ASM2D*>(myModel[j])))
	  {
	    std::cout <<"\tRaising order of P"<< j+1
		      <<" "<< addu <<" "<< addv << std::endl;
	    pch->raiseOrder(addu,addv);
	  }
      }
    }
  }

  else if (!strncasecmp(keyWord,"TOPOLOGY",8))
  {
    if (!this->createFEMmodel()) return false;

    int ntop = atoi(keyWord+8);
    std::cout <<"\nNumber of patch connections: "<< ntop << std::endl;
    std::vector<Interface> top;
    if (discretization == ASM::SplineC1)
      top.reserve(ntop);

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
      else if (discretization == ASM::SplineC1)
	top.push_back(Interface(static_cast<ASMs2DC1*>(mpch),mEdge,
				static_cast<ASMs2DC1*>(spch),sEdge,rever));
    }

    // Second pass for C1-continuous patches, to set up additional constraints
    std::vector<Interface>::const_iterator it;
    for (it = top.begin(); it != top.end(); it++)
      if (!it->slave.first->connectC1(it->slave.second,
				      it->master.first,
				      it->master.second,
				      it->reversed)) return false;
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
    // cannot do multi-threaded assembly with periodicities
#ifdef USE_OPENMP
    omp_set_num_threads(1);
#endif
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
      cline = strtok(NULL," ");
      bool localAxes = cline && toupper(cline[0]) == 'L';
      if (localAxes) cline = strtok(NULL," ");
      int pedge = cline ? atoi(cline) : 0;
      int bcode = (cline = strtok(NULL," ")) ? atoi(cline) : 12;
      double pd = (cline = strtok(NULL," ")) ? atof(cline) : 0.0;

      patch = this->getLocalPatchIndex(patch);
      if (patch < 1) continue;

      int ldim = pedge >= 0 ? 1 : 0;
      if (pedge < 0 || (pedge > 0 && localAxes))
	pedge = -pedge;

      if (pd == 0.0)
      {
	if (!this->addConstraint(patch,pedge,ldim,bcode%1000))
	  return false;
      }
      else
      {
	int code = 1000 + bcode;
	while (myScalars.find(code) != myScalars.end())
	  code += 1000;

	if (!this->addConstraint(patch,pedge,ldim,bcode%1000,code))
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

    ASM2D* pch = 0;
    int nfix = atoi(keyWord+9);
    std::cout <<"\nNumber of fixed points: "<< nfix << std::endl;
    for (int i = 0; i < nfix && (cline = utl::readLine(is)); i++)
    {
      int patch = atoi(strtok(cline," "));
      double rx = atof(strtok(NULL," "));
      double ry = atof(strtok(NULL," "));
      int bcode = (cline = strtok(NULL," ")) ? atoi(cline) : 12;

      int pid = this->getLocalPatchIndex(patch);
      if (pid > 0 && (pch = dynamic_cast<ASM2D*>(myModel[pid-1])))
      {
	std::cout <<"\tConstraining P"<< patch
		  <<" point at "<< rx <<" "<< ry
		  <<" with code "<< bcode << std::endl;
	pch->constrainNode(rx,ry,bcode);
      }
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
            << (ldim == 0 ? " V" : " E") << abs(lndx)
	    <<" in direction(s) "<< dirs;
  if (lndx < 0) // Signals edge constraints in local coordinates
  {
    std::cout <<" (local)";
    preserveNOrder = true; // Because some extra nodes might be added
  }
  if (code) std::cout <<" code = "<< code <<" ";
#if SP_DEBUG > 1
  std::cout << std::endl;
#endif

  // Must dynamic cast here, since ASM2D is not derived from ASMbase
  ASM2D* pch = dynamic_cast<ASM2D*>(myModel[patch-1]);
  if (!pch) return constrError("2D patch ",patch);

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
	case  1: pch->constrainEdge(-1,dirs,code); break;
	case  2: pch->constrainEdge( 1,dirs,code); break;
	case  3: pch->constrainEdge(-2,dirs,code); break;
	case  4: pch->constrainEdge( 2,dirs,code); break;
	case -1: pch->constrainEdgeLocal(-1,dirs,code); break;
	case -2: pch->constrainEdgeLocal( 1,dirs,code); break;
	case -3: pch->constrainEdgeLocal(-2,dirs,code); break;
	case -4: pch->constrainEdgeLocal( 2,dirs,code); break;
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


bool SIM2D::readPatch (std::istream& isp, int pchInd)
{
  ASMbase* pch = ASM2D::create(discretization,nf,mixedFEM);
  if (pch)
  {
    if (!pch->read(isp))
    {
      delete pch;
      return false;
    }
    else if (pch->empty() || this->getLocalPatchIndex(pchInd+1) < 1)
      delete pch;
    else
      myModel.push_back(pch);
  }

  return true;
}


bool SIM2D::readPatches (std::istream& isp)
{
  ASMbase* pch = 0;
  for (int pchInd = 1; isp.good(); pchInd++)
    if ((pch = ASM2D::create(discretization,nf,mixedFEM)))
    {
      std::cout <<"Reading patch "<< pchInd << std::endl;
      if (!pch->read(isp))
      {
	delete pch;
	return false;
      }
      else if (pch->empty() || this->getLocalPatchIndex(pchInd) < 1)
        delete pch;
      else
        myModel.push_back(pch);
    }

  return true;
}


bool SIM2D::readNodes (std::istream& isn, int pchInd, int basis, bool oneBased)
{
  int i;
  ASMs2D::BlockNodes n;

  for (i = 0; i < 4 && isn.good(); i++)
    isn >> n.ibnod[i];
  for (i = 0; i < 4 && isn.good(); i++)
    isn >> n.edges[i].icnod >> n.edges[i].incr;
  isn >> n.iinod;

  if (!isn.good() || pchInd < 0) return true;

  if (!oneBased)
  {
    // We always require the node numbers to be 1-based
    for (i = 0; i < 4; i++) ++n.ibnod[i];
    for (i = 0; i < 4; i++) ++n.edges[i].icnod;
    ++n.iinod;
  }

  return static_cast<ASMs2D*>(myModel[pchInd])->assignNodeNumbers(n,basis);
}


void SIM2D::clonePatches (const FEModelVec& patches,
			  const std::map<int,int>& glb2locN)
{
  ASM2D* pch = 0;
  for (size_t i = 0; i < patches.size(); i++)
    if ((pch = dynamic_cast<ASM2D*>(patches[i])))
      myModel.push_back(pch->clone(nf));

  g2l = &glb2locN;
}


void SIM2D::readNodes (std::istream& isn)
{
  while (isn.good()) {
    int patch = 0;
    isn >> patch;
    ++patch;
    int pid = getLocalPatchIndex(patch);
    if (pid < 0)
      return;

    ASMs2D::BlockNodes n;
    for (size_t i = 0; i < 4 && isn.good(); i++)
      isn >> n.ibnod[i];
    for (size_t i = 0; i < 4 && isn.good(); i++)
      isn >> n.edges[i].icnod >> n.edges[i].incr;
    isn >> n.iinod;

    // We always require the node numbers to be 1-based
    for (size_t i = 0; i < 4; i++)
      ++n.ibnod[i];
    for (size_t i = 0; i < 4; i++)
      ++n.edges[i].icnod;
    ++n.iinod;

    if (isn.good() && pid > 0) {
      if (!static_cast<ASMs2D*>(myModel[pid-1])->assignNodeNumbers(n)) {
        std::cerr <<" *** SIM2D::readNodes: Failed to assign node numbers"
          <<" for patch "<< patch << std::endl;
        return;
      }
    }
  }
}


bool SIM2D::refine (const std::vector<int>& elements,
		    const std::vector<int>& options, const char* fName)
{
  ASM2D* pch = 0;
  for (size_t i = 0; i < myModel.size(); i++)
    if (!myModel[i]->empty() && (pch = dynamic_cast<ASM2D*>(myModel[i])))
      if (!pch->refine(elements,options,fName))
	return false;

  isRefined = true;
  return true;
}
