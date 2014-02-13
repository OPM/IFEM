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
#include "GoTools/geometry/SplineSurface.h"
#ifdef USE_OPENMP
#include <omp.h>
#endif
#include <fstream>
#include <sstream>


/*!
  \brief A struct defining a patch interface for C1-continuous models.
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


SIM2D::SIM2D (unsigned char n1, unsigned char n2, bool)
{
  nf[0] = n1;
  nf[1] = n2;
  nf[2] = 0;
}


SIM2D::SIM2D (IntegrandBase* itg, unsigned char n) : SIMgeneric(itg)
{
  nf[0] = n;
  nf[1] = nf[2] = 0;
}


bool SIM2D::parseGeometryTag (const TiXmlElement* elem)
{
  std::cout <<"  Parsing <"<< elem->Value() <<">"<< std::endl;

  if (!strcasecmp(elem->Value(),"refine") && !isRefined)
  {
    int lowpatch = 1, uppatch = 2;
    if (utl::getAttribute(elem,"patch",lowpatch))
      uppatch = lowpatch;
    if (utl::getAttribute(elem,"lowerpatch",lowpatch))
      uppatch = myModel.size();
    utl::getAttribute(elem,"upperpatch",uppatch);

    if (lowpatch < 1 || uppatch > (int)myModel.size())
    {
      std::cerr <<" *** SIM2D::parse: Invalid patch indices, lower="
                << lowpatch <<" upper="<< uppatch << std::endl;
      return false;
    }

    ASM2D* pch = NULL;
    RealArray xi;
    if (!utl::parseKnots(elem,xi))
    {
      int addu = 0, addv = 0;
      utl::getAttribute(elem,"u",addu);
      utl::getAttribute(elem,"v",addv);
      for (int j = lowpatch-1; j < uppatch; j++)
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
      int dir = 1;
      utl::getAttribute(elem,"dir",dir);
      for (int j = lowpatch-1; j < uppatch; j++)
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

  else if (!strcasecmp(elem->Value(),"raiseorder") && !isRefined)
  {
    int lowpatch = 1, uppatch = 2;
    if (utl::getAttribute(elem,"patch",lowpatch))
      uppatch = lowpatch;
    if (utl::getAttribute(elem,"lowerpatch",lowpatch))
      uppatch = myModel.size();
    utl::getAttribute(elem,"upperpatch",uppatch);

    if (lowpatch < 1 || uppatch > (int)myModel.size())
    {
      std::cerr <<" *** SIM2D::parse: Invalid patch indices, lower="
                << lowpatch <<" upper="<< uppatch << std::endl;
      return false;
    }

    ASM2D* pch = NULL;
    int addu = 0, addv = 0;
    utl::getAttribute(elem,"u",addu);
    utl::getAttribute(elem,"v",addv);
    for (int j = lowpatch-1; j < uppatch; j++)
      if ((pch = dynamic_cast<ASM2D*>(myModel[j])))
      {
        std::cout <<"\tRaising order of P"<< j+1
                  <<" "<< addu <<" "<< addv << std::endl;
        pch->raiseOrder(addu,addv);
      }
  }

  else if (!strcasecmp(elem->Value(),"topology"))
  {
    if (!this->createFEMmodel()) return false;

    std::vector<Interface> top;
    const TiXmlElement* child = elem->FirstChildElement("connection");
    for (; child; child = child->NextSiblingElement())
    {
      int master = 0, slave = 0, mEdge = 0, sEdge = 0;
      bool rever = false;
      utl::getAttribute(child,"master",master);
      utl::getAttribute(child,"medge",mEdge);
      utl::getAttribute(child,"slave",slave);
      utl::getAttribute(child,"sedge",sEdge);
      utl::getAttribute(child,"reverse",rever);

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
      else if (opt.discretization == ASM::SplineC1)
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

  else if (!strcasecmp(elem->Value(),"periodic"))
  {
    if (!this->createFEMmodel()) return false;

    int patch = 0, pedir = 1;
    utl::getAttribute(elem,"patch",patch);
    utl::getAttribute(elem,"dir",pedir);

    if (patch < 1 || patch > (int)myModel.size())
    {
      std::cerr <<" *** SIM2D::parse: Invalid patch index "
                << patch << std::endl;
      return false;
    }
    std::cout <<"\tPeriodic "<< char('H'+pedir) <<"-direction P"<< patch
              << std::endl;
    static_cast<ASMs2D*>(myModel[patch-1])->closeEdges(pedir);
#ifdef USE_OPENMP
    // Cannot do multi-threaded assembly with periodicities
    omp_set_num_threads(1);
#endif
  }

  else if (!strcasecmp(elem->Value(),"immersedboundary"))
  {
    int patch = 0;
    utl::getAttribute(elem,"patch",patch);

    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"Hole") ||
          !strcasecmp(child->Value(),"Circle"))
      {
        double R = 1.0, Xc = 0.0, Yc = 0.0;
        utl::getAttribute(child,"R",R);
        utl::getAttribute(child,"Xc",Xc);
        utl::getAttribute(child,"Yc",Yc);
        if (patch > 0 && patch <= (int)myModel.size())
          myModel[patch-1]->addHole(R,Xc,Yc);
        else
          myModel.front()->addHole(R,Xc,Yc);
      }
      else if (!strcasecmp(child->Value(),"Oval"))
      {
        double R = 1.0, X1 = 0.0, Y1 = 0.0, X2 = 0.0, Y2 = 0.0;
        utl::getAttribute(child,"R",R);
        utl::getAttribute(child,"X1",X1);
        utl::getAttribute(child,"Y1",Y1);
        utl::getAttribute(child,"X2",X2);
        utl::getAttribute(child,"Y2",Y2);
        if (patch > 0 && patch <= (int)myModel.size())
          myModel[patch-1]->addHole(R,X1,Y1,X2,Y2);
        else
          myModel.front()->addHole(R,X1,Y1,X2,Y2);
      }
  }

  return true;
}


bool SIM2D::parseBCTag (const TiXmlElement* elem)
{
  if (!strcasecmp(elem->Value(),"fixpoint") && !ignoreDirichlet)
  {
    if (!this->createFEMmodel()) return false;

    int patch = 0, code = 12;
    double rx = 0.0, ry = 0.0;
    utl::getAttribute(elem,"patch",patch);
    utl::getAttribute(elem,"code",code);
    utl::getAttribute(elem,"rx",rx);
    utl::getAttribute(elem,"ry",ry);
    int pid = this->getLocalPatchIndex(patch);
    if (pid < 1) return pid == 0;

    ASM2D* pch = dynamic_cast<ASM2D*>(myModel[pid-1]);
    if (!pch) return false;

    std::cout <<"\tConstraining P"<< patch
              <<" point at "<< rx <<" "<< ry
              <<" with code "<< code << std::endl;
    pch->constrainNode(rx,ry,code);
  }

  return true;
}


bool SIM2D::parse (const TiXmlElement* elem)
{
  if (!strcasecmp(elem->Value(),"geometry"))
  {
    // Check for immersed boundary calculation.
    // This code must be placed here (and not in parseGeometryTag)
    // due to instanciation of the ASMs2DIB class.
    int maxDepth = 0;
    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"immersedboundary"))
        if (utl::getAttribute(child,"max_depth",maxDepth))
        {
          nf[1] = 'I';
          nf[2] = maxDepth;
          std::cout <<"  Parsing <immersedboundary>\n"
                    <<"\tMax refinement depth : "<< maxDepth << std::endl;
          // Immersed boundary cannot be used with C1-continuous multi-patches
          if (opt.discretization == ASM::SplineC1)
            opt.discretization = ASM::Spline;
      }
  }

  bool result = this->SIMgeneric::parse(elem);

  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(elem->Value(),"geometry"))
      result &= this->parseGeometryTag(child);
    else if (!strcasecmp(elem->Value(),"boundaryconditions"))
      result &= this->parseBCTag(child);

  return result;
}


bool SIM2D::parse (char* keyWord, std::istream& is)
{
  char* cline = NULL;
  if (!strncasecmp(keyWord,"REFINE",6))
  {
    int nref = atoi(keyWord+6);
    if (isRefined) // just read through the next lines without doing anything
      for (int i = 0; i < nref && utl::readLine(is); i++);
    else
    {
      ASM2D* pch = NULL;
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
      ASM2D* pch = NULL;
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
    if (opt.discretization == ASM::SplineC1)
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
      else if (opt.discretization == ASM::SplineC1)
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
    }
#ifdef USE_OPENMP
    // Cannot do multi-threaded assembly with periodicities
    omp_set_num_threads(1);
#endif
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
      int pedge = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
      int bcode = (cline = strtok(NULL," ")) ? atoi(cline) : 12;
      double pd = (cline = strtok(NULL," ")) ? atof(cline) : 0.0;

      patch = this->getLocalPatchIndex(patch);
      if (patch < 1) continue;

      int ldim = pedge >= 0 ? 1 : 0;
      if (pedge < 0) pedge = -pedge;

      if (pd == 0.0)
      {
	if (!this->addConstraint(patch,pedge,ldim,bcode%1000000,0,ngno))
	  return false;
      }
      else
      {
	int code = 1000000 + bcode;
	while (myScalars.find(code) != myScalars.end())
	  code += 1000000;

	if (!this->addConstraint(patch,pedge,ldim,bcode%1000000,-code,ngno))
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

    ASM2D* pch = NULL;
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
    return this->SIMgeneric::parse(keyWord,is);

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


bool SIM2D::addConstraint (int patch, int lndx, int ldim, int dirs, int code,
			   int& ngnod)
{
  if (patch < 1 || patch > (int)myModel.size())
    return constrError("patch index ",patch);

  bool open = ldim < 0; // open means without its end points
  bool project = lndx < -10;
  if (project) lndx += 10;

  std::cout <<"\tConstraining P"<< patch;
  if (abs(ldim) < 2)
    std::cout << (ldim == 0 ? " V" : " E") << abs(lndx);
  std::cout <<" in direction(s) "<< dirs;
  if (lndx < 0) std::cout << (project ? " (local projected)" : " (local)");
  if (code != 0) std::cout <<" code = "<< abs(code) <<" ";
#if SP_DEBUG > 1
  std::cout << std::endl;
#endif

  // Must dynamic cast here, since ASM2D is not derived from ASMbase
  ASM2D* pch = dynamic_cast<ASM2D*>(myModel[patch-1]);
  if (!pch) return constrError("2D patch ",patch);

  switch (abs(ldim))
    {
    case 0: // Vertex constraints
      switch (lndx)
	{
	case 1: pch->constrainCorner(-1,-1,dirs,abs(code)); break;
	case 2: pch->constrainCorner( 1,-1,dirs,abs(code)); break;
	case 3: pch->constrainCorner(-1, 1,dirs,abs(code)); break;
	case 4: pch->constrainCorner( 1, 1,dirs,abs(code)); break;
	default:
	  std::cout << std::endl;
	  return constrError("vertex index ",lndx);
	}
      break;

    case 1: // Edge constraints
      switch (lndx)
	{
	case  1: pch->constrainEdge(-1,open,dirs,code); break;
	case  2: pch->constrainEdge( 1,open,dirs,code); break;
	case  3: pch->constrainEdge(-2,open,dirs,code); break;
	case  4: pch->constrainEdge( 2,open,dirs,code); break;
	case -1:
	  ngnod += pch->constrainEdgeLocal(-1,open,dirs,code,project);
	  break;
	case -2:
	  ngnod += pch->constrainEdgeLocal( 1,open,dirs,code,project);
	  break;
	case -3:
	  ngnod += pch->constrainEdgeLocal(-2,open,dirs,code,project);
	  break;
	case -4:
	  ngnod += pch->constrainEdgeLocal( 2,open,dirs,code,project);
	  break;
	default:
	  std::cout << std::endl;
	  return constrError("edge index ",lndx);
	}
      break;

    case 2: // Face constraints
      myModel[patch-1]->constrainPatch(dirs,code);
      break;

    default:
      std::cout << std::endl;
      return constrError("local dimension switch ",ldim);
    }

  return true;
}


ASMbase* SIM2D::readPatch (std::istream& isp, int pchInd) const
{
  ASMbase* pch = ASM2D::create(opt.discretization,nf,nf[1] > 0);
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


bool SIM2D::readPatches (std::istream& isp, PatchVec& patches,
                         const char* whiteSpace)
{
  ASMbase* pch = NULL;
  for (int pchInd = 1; isp.good(); pchInd++)
    if ((pch = ASM2D::create(opt.discretization,nf,nf[1] > 0)))
    {
      std::cout << whiteSpace <<"Reading patch "<< pchInd << std::endl;
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


void SIM2D::readNodes (std::istream& isn)
{
  while (isn.good())
  {
    int patch = 0;
    isn >> patch;
    int pid = getLocalPatchIndex(patch+1);
    if (pid < 0) return;

    if (!this->readNodes(isn,pid-1))
    {
      std::cerr <<" *** SIM2D::readNodes: Failed to assign node numbers"
		<<" for patch "<< patch+1 << std::endl;
      return;
    }
  }
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


void SIM2D::clonePatches (const PatchVec& patches,
			  const std::map<int,int>& glb2locN)
{
  ASM2D* pch = NULL;
  for (size_t i = 0; i < patches.size(); i++)
    if ((pch = dynamic_cast<ASM2D*>(patches[i])))
      myModel.push_back(pch->clone(nf));

  g2l = &glb2locN;
}


ASMbase* SIM2D::createDefaultGeometry () const
{
  std::istringstream unitSquare("200 1 0 0\n2 0\n"
				"2 2\n0 0 1 1\n"
				"2 2\n0 0 1 1\n"
				"0 0\n"
				"1 0\n"
				"0 1\n"
				"1 1\n");

  return this->readPatch(unitSquare,1);
}


Vector SIM2D::getSolution (const Vector& psol, double u, double v,
                           int deriv, int patch) const
{
  double par[2] = { u, v };
  return this->SIMgeneric::getSolution(psol,par,deriv,patch);
}
