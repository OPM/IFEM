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
#include "FieldFunctions.h"
#include "Functions.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "IFEM.h"
#include "tinyxml.h"
#include <fstream>


/*!
  \brief A struct defining a patch interface for C1-continuous models.
*/

struct Interface
{
  std::pair<ASMs2DC1*,int> master; //!< Patch and edge index of the master
  std::pair<ASMs2DC1*,int> slave;  //!< Patch and edge index of the slave
  bool reversed;                   //!< Relative orientation toggle

  //! \brief Constructor initializing an Interface instance.
  Interface(ASMbase* m, int me, ASMbase* s, int se, bool r)
    : master(std::make_pair(static_cast<ASMs2DC1*>(m),me)),
      slave(std::make_pair(static_cast<ASMs2DC1*>(s),se)), reversed(r) {}
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


bool SIM2D::connectPatches (const ASM::Interface& ifc, bool coordCheck)
{
  if (ifc.master == ifc.slave ||
      ifc.master < 1 || ifc.master > nGlPatches ||
      ifc.slave  < 1 || ifc.slave  > nGlPatches)
  {
    std::cerr <<" *** SIM2D::connectPatches: Invalid patch indices "
              << ifc.master <<" "<< ifc.slave << std::endl;
    return false;
  }

  int lmaster = this->getLocalPatchIndex(ifc.master);
  int lslave  = this->getLocalPatchIndex(ifc.slave);
  if (lmaster > 0 && lslave > 0)
  {
    if (ifc.dim < 1) return true; // ignored in serial

    IFEM::cout <<"\tConnecting P"<< ifc.slave <<" E"<< ifc.sidx
               <<" to P"<< ifc.master <<" E"<< ifc.midx
               <<" reversed? "<< ifc.orient << std::endl;

    ASM2D* spch = dynamic_cast<ASM2D*>(myModel[lslave-1]);
    ASM2D* mpch = dynamic_cast<ASM2D*>(myModel[lmaster-1]);
    if (spch && mpch)
    {
      std::set<int> bases;
      if (ifc.basis == 0)
        for (size_t b = 1; b <= myModel[lslave-1]->getNoBasis(); b++)
          bases.insert(b);
      else
        bases = utl::getDigits(ifc.basis);

      for (int b : bases)
        if (!spch->connectPatch(ifc.sidx,*mpch,ifc.midx,
                                ifc.orient,b,coordCheck,ifc.thick))
          return false;
    }

    myInterfaces.push_back(ifc);
  }
  else
    adm.dd.ghostConnections.insert(ifc);

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

    RealArray xi;
    if (!utl::parseKnots(elem,xi))
    {
      int addu = 0, addv = 0;
      utl::getAttribute(elem,"u",addu);
      utl::getAttribute(elem,"v",addv);
      for (int j : patches)
      {
        IFEM::cout <<"\tRefining P"<< j
                   <<" "<< addu <<" "<< addv << std::endl;
        ASM2D* pch = dynamic_cast<ASM2D*>(this->getPatch(j,true));
        if (pch)
        {
          pch->uniformRefine(0,addu);
          pch->uniformRefine(1,addv);
        }
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
      {
        IFEM::cout <<"\tRefining P"<< j <<" dir="<< dir
                   <<" with grading "<< elem->FirstChild()->Value() <<":";
        for (size_t i = 0; i < xi.size(); i++)
          IFEM::cout << (i%10 || xi.size() < 11 ? " " : "\n\t") << xi[i];
        IFEM::cout << std::endl;
        ASM2D* pch = dynamic_cast<ASM2D*>(this->getPatch(j,true));
        if (pch) pch->refine(dir-1,xi,scale);
      }
    }
  }

  else if (!strcasecmp(elem->Value(),"raiseorder") && !isRefined)
  {
    IntVec patches;
    if (!this->parseTopologySet(elem,patches))
      return false;

    int addu = 0, addv = 0;
    utl::getAttribute(elem,"u",addu);
    utl::getAttribute(elem,"v",addv);
    for (int j : patches)
    {
      IFEM::cout <<"\tRaising order of P"<< j
                 <<" "<< addu <<" "<< addv << std::endl;
      ASM2D* pch = dynamic_cast<ASM2D*>(this->getPatch(j,true));
      if (pch) pch->raiseOrder(addu,addv);
    }
  }

  else if (!strcasecmp(elem->Value(),"topology"))
  {
    if (!this->createFEMmodel()) return false;

    int offset = 0;
    utl::getAttribute(elem,"offset",offset);

    std::vector<Interface> top;
    const TiXmlElement* child = elem->FirstChildElement("connection");
    for (; child; child = child->NextSiblingElement())
    {
      ASM::Interface ifc;
      bool rever = false, periodic = false;
      utl::getAttribute(child,"reverse",rever);
      if (utl::getAttribute(child,"master",ifc.master))
        ifc.master += offset;
      if (!utl::getAttribute(child,"midx",ifc.midx))
        utl::getAttribute(child,"medge",ifc.midx);
      if (utl::getAttribute(child,"slave",ifc.slave))
        ifc.slave += offset;
      if (!utl::getAttribute(child,"sidx",ifc.sidx))
        utl::getAttribute(child,"sedge",ifc.sidx);
      if (!utl::getAttribute(child,"orient",ifc.orient) && rever)
        ifc.orient = 1;
      utl::getAttribute(child,"basis",ifc.basis);
      utl::getAttribute(child,"periodic",periodic);
      if (!utl::getAttribute(child,"dim",ifc.dim))
        ifc.dim = 1;

      if (!this->connectPatches(ifc,!periodic))
        return false;

      if (opt.discretization == ASM::SplineC1)
        top.push_back(Interface(this->getPatch(ifc.master,true),ifc.midx,
                                this->getPatch(ifc.slave,true),ifc.sidx,
                                ifc.orient));
    }

    // Second pass for C1-continuous patches, to set up additional constraints
    for (const Interface& iface : top)
      if (!iface.slave.first->connectC1(iface.slave.second,
                                        iface.master.first,
                                        iface.master.second,
                                        iface.reversed)) return false;
  }

  else if (!strcasecmp(elem->Value(),"periodic"))
    return this->parsePeriodic(elem);

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

    IFEM::cout <<"\tCollapsed edge P"<< patch <<" E"<< edge << std::endl;
    ASMs2D* pch = dynamic_cast<ASMs2D*>(this->getPatch(patch,true));
    if (pch) return pch->collapseEdge(edge);
  }

  else if (!strcasecmp(elem->Value(),"immersedboundary"))
  {
    int patch = 0;
    ASMbase* pch = nullptr;
    if (utl::getAttribute(elem,"patch",patch))
      pch = this->getPatch(patch,true);
    else if (!myModel.empty())
      pch = myModel.front();

    if (!pch)
    {
      std::cerr <<" *** SIM2D::parse: Invalid patch index "
                << patch << std::endl;
      return false;
    }

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
        pch->addHole(R,Xc,Yc);
      }
      else if (!strcasecmp(child->Value(),"Oval"))
      {
        double R = 1.0, X1 = 0.0, Y1 = 0.0, X2 = 0.0, Y2 = 0.0;
        utl::getAttribute(child,"R",R);
        utl::getAttribute(child,"X1",X1);
        utl::getAttribute(child,"Y1",Y1);
        utl::getAttribute(child,"X2",X2);
        utl::getAttribute(child,"Y2",Y2);
        pch->addHole(R,X1,Y1,X2,Y2);
      }
      else if (!strcasecmp(child->Value(),"levelset") && child->FirstChild())
      {
        double power = 1.0, threshold = 0.5;
        std::string type("file");
        utl::getAttribute(child,"power",power);
        utl::getAttribute(child,"threshold",threshold);
        utl::getAttribute(child,"type",type,true);
        const char* val = child->FirstChild()->Value();
        RealFunc* lfunc = nullptr;

        IFEM::cout <<"\tLevel set ("<< type <<")";

        if (type == "file")
        {
          IFEM::cout <<" \""<< val <<"\""<< std::endl;
          std::ifstream ibs(val);
          if (ibs.good())
            lfunc = new FieldFuncStream({pch},ibs);
          else
            std::cerr <<" *** SIM2D::parse: Failed to open file \""<< val
                      <<"\""<< std::endl;
        }
        else
        {
          lfunc = utl::parseRealFunc(val,type);
          IFEM::cout << std::endl;
        }

        if (!pch->setGeometry(lfunc,power,threshold))
          return false;

        if (type == "file")
          myAddScalars[val] = lfunc;
        else
        {
          static std::string fname = "level_set0";
          ++fname[9];
          myAddScalars[fname] = lfunc;
        }
      }
  }

  else if (!strcasecmp(elem->Value(),"projection") && !isRefined)
  {
    const TiXmlElement* child = elem->FirstChildElement();
    if (child && !strncasecmp(child->Value(),"patch",5) && child->FirstChild())
    {
      // Read projection basis from file
      const char* patch = child->FirstChild()->Value();
      std::istream* isp = getPatchStream(child->Value(),patch);
      if (!isp) return false;

      for (ASMbase* pch : myModel)
        pch->createProjectionBasis(false);

      bool ok = true;
      for (int pid = 1; isp->good() && ok; pid++)
      {
        IFEM::cout <<"\tReading projection basis for patch "<< pid << std::endl;
        ASMbase* pch = this->getPatch(pid,true);
        if (pch)
          ok = pch->read(*isp);
        else if ((pch = ASM2D::create(opt.discretization,nsd,nf)))
        {
          // Skip this patch
          ok = pch->read(*isp);
          delete pch;
        }
      }

      delete isp;
      if (!ok) return false;
    }
    else // Generate separate projection basis from current geometry basis
      for (ASMbase* pch : myModel)
        pch->createProjectionBasis(true);

    // Apply refine and/raise-order commands to define the projection basis
    for (; child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"refine") ||
          !strcasecmp(child->Value(),"raiseorder"))
        if (!this->parseGeometryTag(child))
          return false;

    for (ASMbase* pch : myModel)
      if (!pch->createProjectionBasis(false))
      {
        std::cerr <<" *** SIM2D::parseGeometryTag: Failed to create projection"
                  <<" basis, check patch file specification."<< std::endl;
        return false;
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
		for (double u : xi) IFEM::cout <<" "<< u;
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

      if (opt.discretization == ASM::SplineC1)
        top.push_back(Interface(mpch,mEdge,spch,sEdge,rever));
    }

    // Second pass for C1-continuous patches, to set up additional constraints
    for (const Interface& iface : top)
      if (!iface.slave.first->connectC1(iface.slave.second,
                                        iface.master.first,
                                        iface.master.second,
                                        iface.reversed)) return false;
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
      double ry = atof(strtok(nullptr," "));
      int bcode = (cline = strtok(nullptr," ")) ? atoi(cline) : 12;

      ASM2D* pch = dynamic_cast<ASM2D*>(this->getPatch(patch,true));
      if (pch)
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


bool SIM2D::addConstraint (int patch, int lndx, int ldim, int dirs, int code,
                           int& ngnod, char basis)
{
  // Lambda function for error message generation
  auto&& error = [](const char* message, int idx)
  {
    std::cerr <<" *** SIM2D::addConstraint: Invalid "
              << message <<" ("<< idx <<")."<< std::endl;
    return false;
  };

  if (patch < 1 || patch > (int)myModel.size())
    return error("patch index",patch);

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
  if (!pch) return error("2D patch",patch);

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
	  return error("vertex index",lndx);
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
	  return error("edge index",lndx);
	}
      break;

    case 2: // Face constraints
      myModel[patch-1]->constrainPatch(dirs,code);
      break;

    case 4: // Explicit nodal constraints
      myModel[patch-1]->constrainNodes(myModel[patch-1]->getNodeSet(lndx),
				       dirs,code);
      break;

    default:
      IFEM::cout << std::endl;
      return error("local dimension switch",ldim);
    }

  return true;
}


ASMbase* SIM2D::readPatch (std::istream& isp, int pchInd, const CharVec& unf,
                           const char* whiteSpace) const
{
  const CharVec& uunf = unf.empty() ? nf : unf;
  bool isMixed = uunf.size() > 1 && uunf[1] > 0;
  ASMbase* pch = ASM2D::create(opt.discretization,nsd,uunf,isMixed);
  if (pch)
  {
    if (!pch->read(isp) || this->getLocalPatchIndex(pchInd+1) < 1)
    {
      delete pch;
      pch = nullptr;
    }
    else
    {
      if (whiteSpace)
        IFEM::cout << whiteSpace <<"Reading patch "<< pchInd+1 << std::endl;
      if (checkRHSys && dynamic_cast<ASM2D*>(pch)->checkRightHandSystem())
        IFEM::cout <<"\tSwapped."<< std::endl;
      pch->idx = myModel.size();
    }
  }

  return pch;
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
  for (ASMbase* patch : patches)
  {
    ASM2D* pch2D = dynamic_cast<ASM2D*>(patch);
    if (pch2D) myModel.push_back(pch2D->clone(nf));
  }

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


bool SIM2D::writeAddFuncs (int iStep, int& nBlock, int idBlock, double time)
{
  for (const std::pair<std::string,RealFunc*>& func : myAddScalars)
    if (!this->writeGlvF(*func.second, func.first.c_str(), iStep,
                         nBlock, idBlock++, time))
      return false;

  return this->SIMgeneric::writeAddFuncs(iStep,nBlock,idBlock,time);
}
