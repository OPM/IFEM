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
#include "ModelGenerator.h"
#include "ASMs2DC1.h"
#include "ImmersedBoundaries.h"
#include "Functions.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "IFEM.h"
#include "tinyxml.h"
#ifdef USE_OPENMP
#include <omp.h>
#endif


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
    : master(std::make_pair(m,me)), slave(std::make_pair(s,se))
  {
    reversed = r;
  }
};


SIM2D::SIM2D (unsigned char n1, bool check)
{
  nsd = 2;
  nf.push_back(n1);
  checkRHSys = check;
}


SIM2D::SIM2D (const CharVec& fields, bool check) : nf(fields)
{
  nsd = 2;
  checkRHSys = check;
}


SIM2D::SIM2D (IntegrandBase* itg, unsigned char n, bool check) : SIMgeneric(itg)
{
  nsd = 2;
  nf.push_back(n);
  checkRHSys = check;
}


bool SIM2D::addConnection (int master, int slave, int mIdx,
                           int sIdx, int orient, int basis,
                           bool coordCheck, int dim, int thick)
{
  if (orient < 0 || orient > 1)
  {
    std::cerr <<" *** SIM2D::addConnection: Invalid orientation "<< orient <<"."
              << std::endl;
    return false;
  }

  int lmaster = this->getLocalPatchIndex(master);
  int lslave = this->getLocalPatchIndex(slave);
  if (lmaster > 0 && lslave > 0)
  {
    if (dim < 1) return true; // ignored in serial

    IFEM::cout <<"\tConnecting P"<< slave <<" E"<< sIdx
               <<" to P"<< master <<" E"<< mIdx
               <<" reversed? "<< orient << std::endl;

    ASM2D* spch = dynamic_cast<ASM2D*>(myModel[lslave-1]);
    ASM2D* mpch = dynamic_cast<ASM2D*>(myModel[lmaster-1]);

    std::set<int> bases;
    if (basis == 0)
      for (size_t b = 1; b <= myModel[lslave-1]->getNoBasis(); b++)
        bases.insert(b);
    else
      bases = utl::getDigits(basis);

    for (const int& b : bases)
      if (!spch->connectPatch(sIdx,*mpch,mIdx,orient,b,coordCheck,thick))
        return false;

    myInterfaces.push_back(ASM::Interface{master, slave, mIdx, sIdx, orient,
                                          dim, basis, thick});
  }
  else
    adm.dd.ghostConnections.insert(ASM::Interface{master, slave,
                                                  mIdx, sIdx, orient,
                                                  dim, basis, thick});

  return true;
}


bool SIM2D::parseGeometryTag (const TiXmlElement* elem)
{
  IFEM::cout <<"  Parsing <"<< elem->Value() <<">"<< std::endl;

  if (!strcasecmp(elem->Value(),"refine") && !isRefined)
  {
    IntVec patches;
    if (!this->parseTopologySet(elem,patches))
      return false;

    ASM2D* pch = nullptr;
    RealArray xi;
    if (!utl::parseKnots(elem,xi))
    {
      int addu = 0, addv = 0;
      utl::getAttribute(elem,"u",addu);
      utl::getAttribute(elem,"v",addv);
      for (int j : patches)
        if ((pch = dynamic_cast<ASM2D*>(this->getPatch(j,true))))
        {
          IFEM::cout <<"\tRefining P"<< j
                     <<" "<< addu <<" "<< addv << std::endl;
          pch->uniformRefine(0,addu);
          pch->uniformRefine(1,addv);
        }
    }
    else
    {
      // Non-uniform (graded) refinement
      int dir = 1;
      double scale = 1.0;
      utl::getAttribute(elem,"dir",dir);
      utl::getAttribute(elem,"scale",scale);
      for (int j : patches)
        if ((pch = dynamic_cast<ASM2D*>(this->getPatch(j,true))))
        {
          IFEM::cout <<"\tRefining P"<< j <<" dir="<< dir
                     <<" with grading "<< elem->FirstChild()->Value() <<":";
          for (size_t i = 0; i < xi.size(); i++)
            IFEM::cout << (i%10 || xi.size() < 11 ? " " : "\n\t") << xi[i];
          IFEM::cout << std::endl;
          pch->refine(dir-1,xi,scale);
        }
    }
  }

  else if (!strcasecmp(elem->Value(),"raiseorder") && !isRefined)
  {
    IntVec patches;
    if (!this->parseTopologySet(elem,patches))
      return false;

    ASM2D* pch = nullptr;
    int addu = 0, addv = 0;
    utl::getAttribute(elem,"u",addu);
    utl::getAttribute(elem,"v",addv);
    for (int j : patches)
      if ((pch = dynamic_cast<ASM2D*>(this->getPatch(j,true))))
      {
        IFEM::cout <<"\tRaising order of P"<< j
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
      int master = 0, slave = 0, mEdge = 0, sEdge = 0, orient = 0, basis = 0, dim = 1;
      bool rever = false, periodic = false;
      utl::getAttribute(child,"master",master);
      if (!utl::getAttribute(child,"midx",mEdge))
        utl::getAttribute(child,"medge",mEdge);
      utl::getAttribute(child,"slave",slave);
      if (!utl::getAttribute(child,"sidx",sEdge))
        utl::getAttribute(child,"sedge",sEdge);
      if (!utl::getAttribute(child,"orient",orient))
        if (utl::getAttribute(child,"reverse",rever) && rever)
          orient = 1;
      utl::getAttribute(child,"basis",basis);
      utl::getAttribute(child,"periodic",periodic);
      utl::getAttribute(child,"dim",dim);

      if (master == slave ||
          master < 1 || master > nGlPatches ||
          slave  < 1 || slave  > nGlPatches)
      {
        std::cerr <<" *** SIM2D::parse: Invalid patch indices "
                  << master <<" "<< slave << std::endl;
        return false;
      }

      if (!this->addConnection(master,slave,mEdge,sEdge,orient,basis,!periodic,dim))
      {
        std::cerr <<" *** SIM2D::parse: Error establishing connection."
                  << std::endl;
        return false;
      }

      if (opt.discretization == ASM::SplineC1)
        top.push_back(Interface(static_cast<ASMs2DC1*>(this->getPatch(master,true)),mEdge,
                                static_cast<ASMs2DC1*>(this->getPatch(slave,true)),sEdge,orient));
    }

    // Second pass for C1-continuous patches, to set up additional constraints
    for (const Interface& iface : top)
      if (!iface.slave.first->connectC1(iface.slave.second,
                                        iface.master.first,
                                        iface.master.second,
                                        iface.reversed)) return false;
  }

  else if (!strcasecmp(elem->Value(),"periodic"))
  {
    if (!this->createFEMmodel()) return false;

    int patch = 0, pedir = 1;
    utl::getAttribute(elem,"patch",patch);
    utl::getAttribute(elem,"dir",pedir);

    if (patch < 1 || patch > nGlPatches)
    {
      std::cerr <<" *** SIM2D::parse: Invalid patch index "
                << patch << std::endl;
      return false;
    }

    ASMs2D* pch = dynamic_cast<ASMs2D*>(this->getPatch(patch,true));
    if (pch)
    {
      IFEM::cout <<"\tPeriodic "<< char('H'+pedir) <<"-direction P"<< patch
                 << std::endl;
      pch->closeEdges(pedir);
#ifdef USE_OPENMP
      // Cannot do multi-threaded assembly with periodicities
      omp_set_num_threads(1);
#endif
    }
  }

  else if (!strcasecmp(elem->Value(),"collapse"))
  {
    if (!this->createFEMmodel()) return false;

    int patch = 0, edge = 1;
    utl::getAttribute(elem,"patch",patch);
    utl::getAttribute(elem,"edge",edge);

    if (patch < 1 || patch > nGlPatches)
    {
      std::cerr <<" *** SIM2D::parse: Invalid patch index "
                << patch << std::endl;
      return false;
    }

    ASMs2D* pch = dynamic_cast<ASMs2D*>(this->getPatch(patch,true));
    if (pch)
    {
      IFEM::cout <<"\tCollapsed edge P"<< patch <<" E"<< edge << std::endl;
      pch->collapseEdge(edge);
    }
  }

  else if (!strcasecmp(elem->Value(),"immersedboundary"))
  {
    int patch = 0;
    utl::getAttribute(elem,"patch",patch);
    utl::getAttribute(elem,"stabilization",Immersed::stabilization);
    if (Immersed::stabilization != 0)
      IFEM::cout <<"\tStabilization option: "<< Immersed::stabilization
                 << std::endl;

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

    IFEM::cout <<"\tConstraining P"<< patch
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
    // Check for triangular/Matlab mesh or immersed boundary calculation.
    // This code must be placed here (and not in parseGeometryTag)
    // due to instantiation of the ASMs2D[Tri|IB|Matlab] class.
    int maxDepth = 0;
    std::string type;
    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"triangular"))
        opt.discretization = ASM::Triangle;
      else if (!strcasecmp(child->Value(),"patchfile") &&
               utl::getAttribute(child,"type",type,true) && type == "matlab")
      {
        opt.discretization = ASM::Lagrange;
        nf.push_back('M');
      }
      else if (!strcasecmp(child->Value(),"immersedboundary"))
        if (utl::getAttribute(child,"max_depth",maxDepth))
        {
          nf.push_back('I');
          nf.push_back(maxDepth);
          IFEM::cout <<"  Parsing <immersedboundary>\n"
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

  if (myGen && result)
    result = myGen->createTopology(*this);

  delete myGen;
  myGen = nullptr;

  return result;
}


bool SIM2D::parse (char* keyWord, std::istream& is)
{
  char* cline = nullptr;
  if (!strncasecmp(keyWord,"REFINE",6))
  {
    int nref = atoi(keyWord+6);
    if (isRefined) // just read through the next lines without doing anything
      for (int i = 0; i < nref && utl::readLine(is); i++);
    else
    {
      ASM2D* pch = nullptr;
      IFEM::cout <<"\nNumber of patch refinements: "<< nref << std::endl;
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
	  int addu = atoi(strtok(nullptr," "));
	  int addv = atoi(strtok(nullptr," "));
	  for (int j = ipatch; j < patch; j++)
	    if ((pch = dynamic_cast<ASM2D*>(myModel[j])))
	    {
              IFEM::cout <<"\tRefining P"<< j+1
                         <<" "<< addu <<" "<< addv << std::endl;
	      pch->uniformRefine(0,addu);
	      pch->uniformRefine(1,addv);
	    }
	}
	else
	{
	  RealArray xi;
	  int dir = atoi(strtok(nullptr," "));
	  if (utl::parseKnots(xi))
	    for (int j = ipatch; j < patch; j++)
	      if ((pch = dynamic_cast<ASM2D*>(myModel[j])))
	      {
		IFEM::cout <<"\tRefining P"<< j+1 <<" dir="<< dir;
		for (size_t i = 0; i < xi.size(); i++)
		  IFEM::cout <<" "<< xi[i];
		IFEM::cout << std::endl;
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
      ASM2D* pch = nullptr;
      IFEM::cout <<"\nNumber of order raise: "<< nref << std::endl;
      for (int i = 0; i < nref && (cline = utl::readLine(is)); i++)
      {
	int patch = atoi(strtok(cline," "));
	int addu  = atoi(strtok(nullptr," "));
	int addv  = atoi(strtok(nullptr," "));
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
            IFEM::cout <<"\tRaising order of P"<< j+1
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
    IFEM::cout <<"\nNumber of patch connections: "<< ntop << std::endl;
    std::vector<Interface> top;
    if (opt.discretization == ASM::SplineC1)
      top.reserve(ntop);

    for (int i = 0; i < ntop && (cline = utl::readLine(is)); i++)
    {
      int master = atoi(strtok(cline," "));
      int mEdge  = atoi(strtok(nullptr," "));
      int slave  = atoi(strtok(nullptr," "));
      int sEdge  = atoi(strtok(nullptr," "));
      bool rever = (cline = strtok(nullptr," ")) ? cline[0] == 'R' : false;
      if (master == slave ||
	  master < 1 || master > (int)myModel.size() ||
	  slave  < 1 || slave  > (int)myModel.size())
      {
	std::cerr <<" *** SIM2D::parse: Invalid patch indices "
		  << master <<" "<< slave << std::endl;
	return false;
      }
      IFEM::cout <<"\tConnecting P"<< slave <<" E"<< sEdge
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
    for (const Interface& iface : top)
      if (!iface.slave.first->connectC1(iface.slave.second,
                                        iface.master.first,
                                        iface.master.second,
                                        iface.reversed)) return false;
  }

  else if (!strncasecmp(keyWord,"PERIODIC",8))
  {
    if (!this->createFEMmodel()) return false;

    int nper = atoi(keyWord+8);
    IFEM::cout <<"\nNumber of periodicities: "<< nper << std::endl;
    for (int i = 0; i < nper && (cline = utl::readLine(is)); i++)
    {
      int patch = atoi(strtok(cline," "));
      int pedir = atoi(strtok(nullptr," "));
      if (patch < 1 || patch > (int)myModel.size())
      {
	std::cerr <<" *** SIM2D::parse: Invalid patch index "
		  << patch << std::endl;
	return false;
      }
      IFEM::cout <<"\tPeriodic "<< char('H'+pedir) <<"-direction P"<< patch
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
    IFEM::cout <<"\nNumber of constraints: "<< ncon << std::endl;
    for (int i = 0; i < ncon && (cline = utl::readLine(is)); i++)
    {
      int patch = atoi(strtok(cline," "));
      int pedge = (cline = strtok(nullptr," ")) ? atoi(cline) : 0;
      int bcode = (cline = strtok(nullptr," ")) ? atoi(cline) : 12;
      double pd = (cline = strtok(nullptr," ")) ? atof(cline) : 0.0;

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

	IFEM::cout << std::endl;
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

    ASM2D* pch = nullptr;
    int nfix = atoi(keyWord+9);
    IFEM::cout <<"\nNumber of fixed points: "<< nfix << std::endl;
    for (int i = 0; i < nfix && (cline = utl::readLine(is)); i++)
    {
      int patch = atoi(strtok(cline," "));
      double rx = atof(strtok(nullptr," "));
      double ry = atof(strtok(nullptr," "));
      int bcode = (cline = strtok(nullptr," ")) ? atoi(cline) : 12;

      int pid = this->getLocalPatchIndex(patch);
      if (pid > 0 && (pch = dynamic_cast<ASM2D*>(myModel[pid-1])))
      {
        IFEM::cout <<"\tConstraining P"<< patch
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
                           int& ngnod, char basis)
{
  if (patch < 1 || patch > (int)myModel.size())
    return constrError("patch index ",patch);

  bool open = ldim < 0; // open means without its end points
  bool project = lndx < -10;
  if (project) lndx += 10;

  IFEM::cout <<"\tConstraining P"<< patch;
  if (abs(ldim) < 2)
    IFEM::cout << (ldim == 0 ? " V" : " E") << abs(lndx);
  IFEM::cout <<" in direction(s) "<< dirs;
  if (lndx < 0) IFEM::cout << (project ? " (local projected)" : " (local)");
  if (code != 0) IFEM::cout <<" code = "<< abs(code);
  if (basis > 1) IFEM::cout <<" basis = "<< (int)basis;
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
	case 1: pch->constrainCorner(-1,-1,dirs,abs(code),basis); break;
	case 2: pch->constrainCorner( 1,-1,dirs,abs(code),basis); break;
	case 3: pch->constrainCorner(-1, 1,dirs,abs(code),basis); break;
	case 4: pch->constrainCorner( 1, 1,dirs,abs(code),basis); break;
	default:
	  IFEM::cout << std::endl;
	  return constrError("vertex index ",lndx);
	}
      break;

    case 1: // Edge constraints
      switch (lndx)
	{
	case  1: pch->constrainEdge(-1,open,dirs,code,basis); break;
	case  2: pch->constrainEdge( 1,open,dirs,code,basis); break;
	case  3: pch->constrainEdge(-2,open,dirs,code,basis); break;
	case  4: pch->constrainEdge( 2,open,dirs,code,basis); break;
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
	  IFEM::cout << std::endl;
	  return constrError("edge index ",lndx);
	}
      break;

    case 2: // Face constraints
      myModel[patch-1]->constrainPatch(dirs,code);
      break;

    case 4: // Explicit nodal constrains
      myModel[patch-1]->constrainNodes(myModel[patch-1]->getNodeSet(lndx),
				       dirs,code);
      break;

    default:
      IFEM::cout << std::endl;
      return constrError("local dimension switch ",ldim);
    }

  return true;
}


ASMbase* SIM2D::readPatch (std::istream& isp, int pchInd,
                           const CharVec& unf) const
{
  const CharVec& uunf = unf.empty() ? nf : unf;
  bool isMixed = uunf.size() > 1 && uunf[1] > 0;
  ASMbase* pch = ASM2D::create(opt.discretization,nsd,uunf,isMixed);
  if (pch)
  {
    if (!pch->read(isp))
      delete pch, pch = nullptr;
    else if (pch->empty() || this->getLocalPatchIndex(pchInd+1) < 1)
      delete pch, pch = nullptr;
    else
      pch->idx = myModel.size();
  }

  if (checkRHSys && pch)
    if (dynamic_cast<ASM2D*>(pch)->checkRightHandSystem())
      IFEM::cout <<"\tSwapped."<< std::endl;

  return pch;
}


bool SIM2D::readPatches (std::istream& isp, PatchVec& patches,
                         const char* whiteSpace) const
{
  ASMbase* pch = nullptr;
  bool isMixed = nf.size() > 1 && nf[1] > 0;
  for (int pchInd = 1; isp.good(); pchInd++)
    if ((pch = ASM2D::create(opt.discretization,nsd,nf,isMixed)))
    {
      if (!pch->read(isp))
      {
        delete pch;
        return false;
      }
      else if (pch->empty() || this->getLocalPatchIndex(pchInd) < 1)
        delete pch;
      else
      {
        if (whiteSpace)
          IFEM::cout << whiteSpace <<"Reading patch "<< pchInd << std::endl;
        pch->idx = patches.size();
        patches.push_back(pch);
        if (checkRHSys)
          if (dynamic_cast<ASM2D*>(pch)->checkRightHandSystem())
            IFEM::cout <<"\tSwapped."<< std::endl;
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
    int pid = this->getLocalPatchIndex(patch+1);
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
  ASM2D* pch = nullptr;
  for (size_t i = 0; i < patches.size(); i++)
    if ((pch = dynamic_cast<ASM2D*>(patches[i])))
      myModel.push_back(pch->clone(nf));

  g2l = &glb2locN;

  if (nGlPatches == 0)
    nGlPatches = myModel.size();
}


ModelGenerator* SIM2D::getModelGenerator (const TiXmlElement* geo) const
{
  return new DefaultGeometry2D(geo);
}


Vector SIM2D::getSolution (const Vector& psol, double u, double v,
                           int deriv, int patch) const
{
  double par[2] = { u, v };
  return this->SIMgeneric::getSolution(psol,par,deriv,patch);
}
