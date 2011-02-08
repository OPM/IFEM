// $Id: LinearEl.C,v 1.32 2010-12-07 12:56:25 kmo Exp $
//==============================================================================
//!
//! \file LinearEl.C
//!
//! \date Jan 23 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solver for linear elasticity problems using NURBS-based FEM.
//!
//==============================================================================

#include "LinearEl.h"
#include "EigSolver.h"
#include "VolumePatch.h"
#include "ElementBlock.h"
#include "AnalyticSolutions.h"
#include "Functions.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "SAMSpline.h"
#include "Utilities.h"
#include "VTF.h"
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>


/*!
  \brief Local coordinate system for a cylinder along global z-axis.
*/

class CylinderCS : public LocalSystem
{
public:
  //! \brief Constructor printing a message making user aware of its presense.
  CylinderCS()
  {
    std::cout <<"\nLocal coordinate system: Cylindric"<< std::endl;
  }
  //! \brief Computes global-to-local transformation at the point \a X.
  virtual const Tensor& getTmat(const Vec3& X) const
  {
    static Tensor T(3);
    double r = hypot(X.x,X.y);
    T(1,1) = X.x/r;
    T(1,2) = X.y/r;
    T(2,1) = -T(1,2);
    T(2,2) = T(1,1);
    T(3,3) = 1.0;
    return T;
  }
};


LinearEl::LinearEl (const char* fileName, bool checkRHS, bool free)
{
  rCS = 0;
  sam = 0;
  asol = 0;

  if (fileName)
    if (!this->read(fileName,checkRHS,free))
      exit(1);
}


LinearEl::~LinearEl ()
{
  if (rCS) delete rCS;
  if (sam) delete sam;
  if (asol) delete asol;

  for (size_t i = 0; i < model.size(); i++)
    delete model[i];
}


bool LinearEl::read (const char* fileName, bool checkRHS, bool free)
{
  MPCLess::compareSlaveDofOnly = true; // to avoid multiple slave definitions

  model.clear();
  load.clear();

  char* cline = 0;
  std::cout <<"\nReading input file "<< fileName << std::endl;
  std::ifstream is(fileName);
  while (is.good() && (cline = utl::readLine(is)))

    if (!strncasecmp(cline,"PATCHES",7))
    {
      int npatch = atoi(cline+7);
      std::cout <<"\nNumber of patches: "<< npatch << std::endl;
      for (int i = 0; i < npatch && (cline = utl::readLine(is)); i++)
      {
	VolumePatch* vp = new VolumePatch(strtok(cline," "),checkRHS);
	if (vp->empty())
	  delete vp;
	else
	  model.push_back(vp);
      }
      if (model.size() < npatch)
      {
	std::cerr <<" *** LinearEl::read: Expected "<< npatch
		  <<" patches but could read only "<< model.size() << std::endl;
	return false;
      }
    }

    else if (!strncasecmp(cline,"PATCHFILE",9))
    {
      size_t i = 9; while (i < strlen(cline) && isspace(cline[i])) i++;
      std::cout <<"\nReading data file "<< cline+i << std::endl;
      std::ifstream isp(cline+i);
      for (int patchNo = 1; isp.good(); patchNo++)
      {
	std::cout <<"Reading patch "<< patchNo << std::endl;
	VolumePatch* vp = new VolumePatch(isp,checkRHS);
	if (vp->empty())
	  delete vp;
	else
	  model.push_back(vp);
      }
      if (model.empty())
      {
	std::cerr <<" *** LinearEl::read: No patches read"<< std::endl;
	return false;
      }
    }

    else if (!strncasecmp(cline,"REFINE",6))
    {
      int nref = atoi(cline+6);
      std::cout <<"\nNumber of patch refinements: "<< nref << std::endl;
      for (int i = 0; i < nref && (cline = utl::readLine(is)); i++)
      {
	int patch = atoi(strtok(cline," "));
	int addu  = atoi(strtok(NULL," "));
	int addv  = atoi(strtok(NULL," "));
	int addw  = atoi(strtok(NULL," "));
	if (patch == 0 || abs(patch) > model.size())
	{
	  std::cerr <<" *** LinearEl::read: Invalid patch index "
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
	  std::cout <<"\tRefining P"<< j+1
		    <<" "<< addu <<" "<< addv <<" "<< addw << std::endl;
	  if (addu > 0) model[j]->uniformRefine(0,addu);
	  if (addv > 0) model[j]->uniformRefine(1,addv);
	  if (addw > 0) model[j]->uniformRefine(2,addw);
	}
      }
    }

    else if (!strncasecmp(cline,"RAISEORDER",10))
    {
      int nref = atoi(cline+10);
      std::cout <<"\nNumber of order raise: "<< nref << std::endl;
      for (int i = 0; i < nref && (cline = utl::readLine(is)); i++)
      {
	int patch = atoi(strtok(cline," "));
	int addu  = atoi(strtok(NULL," "));
	int addv  = atoi(strtok(NULL," "));
	int addw  = atoi(strtok(NULL," "));
	if (patch == 0 || abs(patch) > model.size())
	{
	  std::cerr <<" *** LinearEl::read: Invalid patch index "
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
	  model[j]->raiseOrder(addu,addv,addw);
	}
      }
    }

    else if (!strncasecmp(cline,"TOPOLOGYFILE",12))
    {
      if (!this->createFEMmodel()) return false;

      size_t i = 12; while (i < strlen(cline) && isspace(cline[i])) i++;
      std::cout <<"\nReading data file "<< cline+i << std::endl;
      std::ifstream ist(cline+i);
      while (cline = utl::readLine(ist))
      {
	int master = atoi(strtok(cline," "));
	int mFace  = atoi(strtok(NULL," "))+1;
	int slave  = atoi(strtok(NULL," "));
	int sFace  = atoi(strtok(NULL," "))+1;
	int swapd  = atoi(strtok(NULL," "));
	int rev_u  = atoi(strtok(NULL," "));
	int rev_v  = atoi(strtok(NULL," "));
	int orient = 4*swapd+2*rev_u+rev_v;
	if (master == slave ||
	    master < 0 || master >= model.size() ||
	    slave  < 0 || slave  >= model.size())
	{
	  std::cerr <<" *** LinearEl::read: Invalid patch indices "
		    << master <<" "<< slave << std::endl;
	  return false;
	}
	std::cout <<"\tConnecting P"<< slave+1 <<" F"<< sFace
		  <<" to P"<< master+1 <<" F"<< mFace
		  <<" orient "<< orient << std::endl;
	if (!model[slave]->connectPatch(sFace,*model[master],mFace,orient))
	  return false;
      }
    }

    else if (!strncasecmp(cline,"TOPOLOGY",8))
    {
      if (!this->createFEMmodel()) return false;

      int ntop = atoi(cline+8);
      std::cout <<"\nNumber of patch connections: "<< ntop << std::endl;
      for (int i = 0; i < ntop && (cline = utl::readLine(is)); i++)
      {
	int master = atoi(strtok(cline," "))-1;
	int mFace  = atoi(strtok(NULL," "));
	int slave  = atoi(strtok(NULL," "))-1;
	int sFace  = atoi(strtok(NULL," "));
	int orient = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
	if (master == slave ||
	    master < 0 || master >= model.size() ||
	    slave  < 0 || slave  >= model.size())
	{
	  std::cerr <<" *** LinearEl::read: Invalid patch indices "
		    << master <<" "<< slave << std::endl;
	  return false;
	}
	std::cout <<"\tConnecting P"<< slave+1 <<" F"<< sFace
		  <<" to P"<< master+1 <<" F"<< mFace
		  <<" orient "<< orient << std::endl;
	if (!model[slave]->connectPatch(sFace,*model[master],mFace,orient))
	  return false;
      }
    }

    else if (!strncasecmp(cline,"PERIODIC",8))
    {
      if (!this->createFEMmodel()) return false;

      int nper = atoi(cline+8);
      std::cout <<"\nNumber of periodicities: "<< nper << std::endl;
      for (int i = 0; i < nper && (cline = utl::readLine(is)); i++)
      {
	int patch = atoi(strtok(cline," "))-1;
	int pfdir = atoi(strtok(NULL," "));
	if (patch < 0 || patch >= model.size())
	{
	  std::cerr <<" *** LinearEl::read: Invalid patch index "
		    << patch << std::endl;
	  return false;
	}
	std::cout <<"\tPeriodic "<< char('H'+pfdir) <<"-direction P"<< patch+1
		  << std::endl;
	model[patch]->closeFaces(pfdir);
      }
    }

    else if (!strncasecmp(cline,"CONSTRAINTS",11))
    {
      if (free) continue; // Ignore all boundary conditions
      if (!this->createFEMmodel()) return false;

      int ncon = atoi(cline+11);
      std::cout <<"\nNumber of constraints: "<< ncon << std::endl;
      for (int i = 0; i < ncon && (cline = utl::readLine(is)); i++)
      {
	int patch = atoi(strtok(cline," "))-1;
	int pface = atoi(strtok(NULL," "));
	int pline = pface > 10 ? pface/10 : 0;
	int bcode = atoi(strtok(NULL," "));
	double pd = (cline = strtok(NULL," ")) ? atof(cline) : 0.0;
	if (pface > 10) pface = pface%10;
	if (patch < 0 || patch >= model.size())
	{
	  std::cerr <<" *** LinearEl::read: Invalid patch index "
		    << patch << std::endl;
	  return false;
	}
	std::cout <<"\tConstraining P"<< patch+1;
	if (pline > 0)
	  std::cout <<" F"<< pface <<" L"<< pline <<" at xi="<< pd
		    <<" with code "<< bcode;
	else
	{
	  std::cout << (pface < 0 ? " C" : " F") << abs(pface)
		    <<" with code "<< bcode;
	  if (pd != 0.0) std::cout <<" value = "<< pd;
	}
	std::cout << std::endl;
	if (pline == 1) switch (pface)
	  {
	  case 1: model[patch]->constrainLine(-1,2,pd,bcode); break;
	  case 2: model[patch]->constrainLine( 1,2,pd,bcode); break;
	  case 3: model[patch]->constrainLine(-2,3,pd,bcode); break;
	  case 4: model[patch]->constrainLine( 2,3,pd,bcode); break;
	  case 5: model[patch]->constrainLine(-3,1,pd,bcode); break;
	  case 6: model[patch]->constrainLine( 3,1,pd,bcode); break;
	  }
	else if (pline == 2) switch (pface)
	  {
	  case 1: model[patch]->constrainLine(-1,3,pd,bcode); break;
	  case 2: model[patch]->constrainLine( 1,3,pd,bcode); break;
	  case 3: model[patch]->constrainLine(-2,1,pd,bcode); break;
	  case 4: model[patch]->constrainLine( 2,1,pd,bcode); break;
	  case 5: model[patch]->constrainLine(-3,2,pd,bcode); break;
	  case 6: model[patch]->constrainLine( 3,2,pd,bcode); break;
	  }
	else switch (pface)
	  {
	  case 1: model[patch]->constrainFace(-1,bcode,pd); break;
	  case 2: model[patch]->constrainFace( 1,bcode,pd); break;
	  case 3: model[patch]->constrainFace(-2,bcode,pd); break;
	  case 4: model[patch]->constrainFace( 2,bcode,pd); break;
	  case 5: model[patch]->constrainFace(-3,bcode,pd); break;
	  case 6: model[patch]->constrainFace( 3,bcode,pd); break;
	  default:
	    switch (-pface)
	      {
	      case 1: model[patch]->constrainCorner(-1,-1,-1,bcode,pd); break;
	      case 2: model[patch]->constrainCorner( 1,-1,-1,bcode,pd); break;
	      case 3: model[patch]->constrainCorner(-1, 1,-1,bcode,pd); break;
	      case 4: model[patch]->constrainCorner( 1, 1,-1,bcode,pd); break;
	      case 5: model[patch]->constrainCorner(-1,-1, 1,bcode,pd); break;
	      case 6: model[patch]->constrainCorner( 1,-1, 1,bcode,pd); break;
	      case 7: model[patch]->constrainCorner(-1, 1, 1,bcode,pd); break;
	      case 8: model[patch]->constrainCorner( 1, 1, 1,bcode,pd); break;
	      default:
		std::cerr <<" *** LinearEl::read: Invalid face index "
			  << pface << std::endl;
	      }
	  }
      }
    }

    else if (!strncasecmp(cline,"FIXPOINTS",9))
    {
      if (free) continue; // Ignore all boundary conditions
      if (!this->createFEMmodel()) return false;

      int nfix = atoi(cline+9);
      std::cout <<"\nNumber of fixed points: "<< nfix << std::endl;
      for (int i = 0; i < nfix && (cline = utl::readLine(is)); i++)
      {
	int patch = atoi(strtok(cline," "))-1;
	double rx = atof(strtok(NULL," "));
	double ry = atof(strtok(NULL," "));
	double rz = atof(strtok(NULL," "));
	int bcode = (cline = strtok(NULL," ")) ? atoi(cline) : 123;
	if (patch < 0 || patch >= model.size())
	{
	  std::cerr <<" *** LinearEl::read: Invalid patch index "
		    << patch << std::endl;
	  return false;
	}
	std::cout <<"\tConstraining P"<< patch+1
		  <<" point at "<< rx <<" "<< ry <<" "<< rz
		  <<" with code "<< bcode << std::endl;
	model[patch]->constrainNode(rx,ry,rz,bcode);
      }
    }

    else if (!strncasecmp(cline,"ANASOL",6))
    {
      cline = strtok(cline+6," ");
      if (!strncasecmp(cline,"HOLE",4))
      {
	double a  = atof(strtok(NULL," "));
	double F0 = atof(strtok(NULL," "));
	double nu = atof(strtok(NULL," "));
	asol = new Hole(a,F0,nu,true);
	std::cout <<"\nAnalytical solution: Hole a="<< a <<" F0="<< F0
		  <<" nu="<< nu << std::endl;
      }
      else if (!strncasecmp(cline,"LSHAPE",6))
      {
	double a  = atof(strtok(NULL," "));
	double F0 = atof(strtok(NULL," "));
	double nu = atof(strtok(NULL," "));
	asol = new Lshape(a,F0,nu,true);
	std::cout <<"\nAnalytical solution: Lshape a="<< a <<" F0="<< F0
		  <<" nu="<< nu << std::endl;
      }
      else if (!strncasecmp(cline,"CANTS",5))
      {
	double L  = atof(strtok(NULL," "));
	double H  = atof(strtok(NULL," "));
	double F0 = atof(strtok(NULL," "));
	asol = new CanTS(L,H,F0,true);
	std::cout <<"\nAnalytical solution: CanTS L="<< L <<" H="<< H
		  <<" F0="<< F0 << std::endl;
      }
      else
	std::cerr <<"  ** LinearEl::read: Unknown analytical solution "
		  << cline << std::endl;
    }

    else if (!strncasecmp(cline,"PRESSURE",8))
    {
      if (free) continue; // Ignore all boundary conditions

      int npres = atoi(cline+8);
      std::cout <<"\nNumber of pressures: "<< npres << std::endl;
      for (int i = 0; i < npres && (cline = utl::readLine(is)); i++)
      {
	int patch = atoi(strtok(cline," "))-1;
	int pface = atoi(strtok(NULL," "));
	if (patch < 0 || patch >= model.size())
	{
	  std::cerr <<" *** LinearEl::read: Invalid patch index "
		    << patch << std::endl;
	  return false;
	}
	int pdir  = -1;
	double pv = 0.0;
	if (asol)
	  std::cout <<"\tTraction on P"<< patch+1 <<" F" << pface << std::endl;
	else
	{
	  pdir = atoi(strtok(NULL," "));
	  pv   = atof(strtok(NULL," "));
	  std::cout <<"\tPressure on P"<< patch+1
		    <<" F" << pface <<" direction "<< pdir
		    <<": "<< pv << std::endl;
	}
	switch (pface)
	  {
	  case 1: load.push_back(PLoad(patch,-1,pdir,pv)); break;
	  case 2: load.push_back(PLoad(patch, 1,pdir,pv)); break;
	  case 3: load.push_back(PLoad(patch,-2,pdir,pv)); break;
	  case 4: load.push_back(PLoad(patch, 2,pdir,pv)); break;
	  case 5: load.push_back(PLoad(patch,-3,pdir,pv)); break;
	  case 6: load.push_back(PLoad(patch, 3,pdir,pv)); break;
	  default:
	    std::cerr <<" *** LinearEl::read: Invalid face index "
		      << pface << std::endl;
	  }
      }
    }

    else if (!strncasecmp(cline,"MATERIAL",8))
    {
      int nmat = atoi(cline+8);
      std::cout <<"\nNumber of materials: "<< nmat << std::endl;
      for (int i = 0; i < nmat && (cline = utl::readLine(is)); i++)
      {
        int last = 0;
	double E = atof(strtok(cline," "));
	double nu = atof(strtok(NULL," "));
	double rho = atof(strtok(NULL," "));
	while (cline = strtok(NULL," "))
	  if (!strncasecmp(cline,"ALL",3))
	  {
	    std::cout <<"\tMaterial for all patches: "
		      << E <<" "<< nu <<" "<< rho << std::endl;
	    for (size_t j = 0; j < model.size(); j++)
	      model[j]->setMaterial(E,nu,rho);
	  }
	  else
	  {
	    int patch = atoi(cline)-1;
	    if (patch < 0 || patch >= model.size())
	    {
	      std::cerr <<" *** LinearEl::read: Invalid patch index "
			<< patch << std::endl;
	      return false;
	    }
	    std::cout <<"\tMaterial for P"<< patch+1
		      <<": "<< E <<" "<< nu <<" "<< rho << std::endl;
	    model[patch]->setMaterial(E,nu,rho);
	  }
      }
    }

    else if (!strncasecmp(cline,"GRAVITY",7))
    {
      gravity.x = atof(strtok(cline+7," "));
      gravity.y = atof(strtok(NULL," "));
      gravity.z = atof(strtok(NULL," "));
      std::cout <<"\nGravitation vector: "<< gravity << std::endl;
    }

    else if (!strncasecmp(cline,"LOCAL_SYSTEM",12))
    {
      size_t i = 12;
      while (i < strlen(cline) && isspace(cline[i])) i++;
      if (!strncasecmp(cline+i,"CYLINDRICZ",10))
	rCS = new CylinderCS;
      else
	std::cerr <<" *** LinearEl::read: Unsupported coordinate system: "
		  << cline+i << std::endl;
    }

    else if (isalpha(cline[0]))
      std::cerr <<" *** LinearEl::read: Unknown keyword: "<< cline << std::endl;

  return this->createFEMmodel();
}


bool LinearEl::createFEMmodel ()
{
  for (size_t i = 0; i < model.size(); i++)
    if (!model[i]->generateFEMTopology())
      return false;

  return true;
}


bool LinearEl::preprocess (const std::vector<int>& ignoredPatches, bool fixDup)
{
  // Erase all patches that should be ignored in the analysis
  std::vector<int>::const_iterator it;
  for (it = ignoredPatches.begin(); it != ignoredPatches.end(); it++)
    if (*it > 0 && *it <= model.size())
      model[*it-1]->clear();

  // Renumber the nodes to account for overlapping nodes and erased patches
  int nnod = VolumePatch::renumberNodes(model);

  if (fixDup && VolumePatch::mergeDuplNodes)
  {
    // Check for duplicated nodes (missing topology)
    int nDupl = 0;
    std::map<Vec3,int> globalNodes;
    for (size_t i = 0; i < model.size(); i++)
    {
      std::map<int,int> old2new;
      for (int node = 1; node <= model[i]->getNoNodes(); node++)
      {
	int globNum = model[i]->getNodeID(node);
	Vec3 X(model[i]->getCoord(node));
	std::map<Vec3,int>::const_iterator xit = globalNodes.find(X);
	if (xit == globalNodes.end())
	  globalNodes.insert(std::make_pair(X,globNum));
	else if (xit->second != globNum)
	{
	  std::cout <<"  ** Merging duplicated nodes "<< xit->second <<" and "
		    << globNum <<" at X="<< X << std::endl;
	  if (model[i]->mergeNodes(node,xit->second))
	    old2new[globNum] = xit->second;
	}
      }
      if (!old2new.empty())
      {
	model[i]->renumberNodes(old2new,true);
	nDupl += old2new.size();
      }
    }
    if (nDupl > 0)
    {
      std::cout <<"   * "<< nDupl <<" duplicated nodes merged."<< std::endl;
      // Renumber the nodes to account for the merged nodes
      nnod = VolumePatch::renumberNodes(model);
    }
  }

  // Resolve possibly chained multi-point constraints
  VolumePatch::resolveMPCchains(model);

  // Initialize data structures for the algebraic system
  sam = new SAMSpline();
  return sam->init(model, VolumePatch::mergeDuplNodes ? nnod : 0);
}


bool LinearEl::assembleKandR (SystemMatrix::Type solver, int nGauss, Vector* R)
{
  if (!sam) return false;
#if SP_DEBUG > 1
  sam->print(std::cout);
#endif

  sys.clear();
  sys.K = SystemMatrix::create(solver);
  if (!sam->initForAssembly(sys,true))
    return true; // no free dofs in system

  // Assemble the system stiffness matrix
  size_t i;
  bool ok = true;
  for (i = 0; i < model.size() && ok; i++)
  {
    std::cout <<"\nAssembling stiffness matrix for P"<< i+1 << std::endl;
    ok = model[i]->assembleSystem(sys,*sam,gravity,nGauss);
  }

  // Assemble right-hand-side contributions from boundary tractions
  for (i = 0; i < load.size() && ok; i++)
  {
    int j = load[i].volp;
    std::cout <<"\nAssembling "<< (load[i].pdir >= 0 ? "pressure" : "traction")
	      <<" forces on face "<< load[i].face <<" for P"<< j+1 << std::endl;

    if (load[i].pdir >= 0)
    {
      // Constant surface pressure
      const RealFunc* p = new ConstFunc(load[i].p);
      ok = model[j]->assembleForces(sys.RHS,*sam,PressureField(p,load[i].pdir),
				    load[i].face,nGauss,&trac);
    }
    else if (asol)
      // Use traction field from the analytical solution
      ok = model[j]->assembleForces(sys.RHS,*sam,TractionField(*asol),
				    load[i].face,nGauss,&trac);
  }

#if SP_DEBUG > 3
  std::cout <<"System stiffness matrix:"<< *sys.K;
  std::cout <<"System force vector:"<< sys.RHS;
#endif

  // Store the right-hand-side vector for visualization
  if (ok && R)
    ok = sam->expandSolution(sys.RHS,*R);

  return ok;
}


bool LinearEl::assembleKandM (SystemMatrix::Type solver, int nGauss)
{
  if (!sam) return false;
#if SP_DEBUG > 1
  sam->print(std::cout);
#endif

  sys.clear();
  sys.K = SystemMatrix::create(solver);
  sys.M = SystemMatrix::create(solver);
  if (!sam->initForAssembly(sys))
    return true; // no free dofs in system

  // Assemble the system stiffness and mass matrices
  bool ok = true;
  for (size_t i = 0; i < model.size() && ok; i++)
  {
    std::cout <<"\nAssembling stiffness/mass matrix for P"<< i+1 << std::endl;
    ok = model[i]->assembleSystem(sys,*sam,Vec3(),nGauss);
  }

#if SP_DEBUG > 3
  std::cout <<"System stiffness matrix:"<< *sys.K;
  std::cout <<"System mass matrix:"<< *sys.M;
#endif
  return ok;
}


bool LinearEl::assembleKandKg (SystemMatrix::Type solver, int nGauss,
			       const Vector& dis)
{
  if (!sam) return false;
#if SP_DEBUG > 1
  sam->print(std::cout);
#endif

  sys.clear();
  sys.K = SystemMatrix::create(solver);
  sys.M = SystemMatrix::create(solver);
  if (!sam->initForAssembly(sys))
    return true; // no free dofs in system

  // Assemble the material- and geometric system stiffness matrices
  bool ok = true;
  Vector displ;
  for (size_t i = 0; i < model.size() && ok; i++)
  {
    std::cout <<"\nAssembling stiffness matrices for P"<< i+1 << std::endl;
    model[i]->extractSolution(dis,displ);
    ok = model[i]->assembleSystem(sys,*sam,Vec3(),nGauss,displ);
  }

#if SP_DEBUG > 3
  std::cout <<"System stiffness matrix:"<< *sys.K;
  std::cout <<"Geometric stiffness matrix:"<< *sys.M;
#endif
  return ok;
}


bool LinearEl::assembleKonly (SystemMatrix::Type solver, int nGauss)
{
  if (!sam) return false;
#if SP_DEBUG > 1
  sam->print(std::cout);
#endif

  sys.clear();
  sys.K = SystemMatrix::create(solver);
  if (!sam->initForAssembly(sys))
    return true; // no free dofs in system

  // Assemble the system stiffness matrix
  bool ok = true;
  for (size_t i = 0; i < model.size() && ok; i++)
  {
    std::cout <<"\nAssembling stiffness matrix for P"<< i+1 << std::endl;
    ok = model[i]->assembleSystem(sys,*sam,Vec3(),nGauss);
  }

#if SP_DEBUG > 3
  std::cout <<"System stiffness matrix:"<< *sys.K;
#endif
  return ok;
}


bool LinearEl::solve (Vector& solution)
{
  if (!sys.K) return false;
  if (!sam) return false;

  // Solve the linear system of equations
  std::cout <<"\nSolving the equation system ..."<< std::endl;
  if (!sys.K->solve(sys.RHS)) return false;
  if (!sam->expandSolution(sys.RHS,solution)) return false;

  size_t iXmax = 0;
  size_t iYmax = 1;
  size_t iZmax = 2;
  double dXmax = solution.normInf(iXmax,3);
  double dYmax = solution.normInf(iYmax,3);
  double dZmax = solution.normInf(iZmax,3);
  double dNorm = solution.norm2();

  std::cout <<"\n >>> Solution summary <<<\n"
	    <<"\nEuclidean norm     : "<< dNorm
	    <<"\nL2-norm            : "<< dNorm/sqrt(solution.size())
	    <<"\nMax X-displacement : "<< dXmax <<" node "<< iXmax
	    <<"\nMax Y-displacement : "<< dYmax <<" node "<< iYmax
	    <<"\nMax Z-displacement : "<< dZmax <<" node "<< iZmax
	    << std::endl;

  if (sam->getNoEquations() < 300)
  {
    std::cout <<"\nSolution vector:";
    for (size_t i = 0; i < solution.size(); i += 3)
      std::cout <<"\nNode "<< i/3+1 <<": "
		<< solution[i] <<" "<< solution[i+1] <<" "<< solution[i+2];
    std::cout << std::endl;
  }

  return true;
}


bool LinearEl::solutionNorms (Matrix& eNorm, int nGauss, const Vector& dis)
{
  Vector gNorm(asol ? 3 : 1);
  eNorm.resize(gNorm.size(),sam->getNoElms(),true);

  size_t i;
  bool ok = true;
  Vector displ;
  for (i = 0; i < model.size() && ok; i++)
  {
    model[i]->extractSolution(dis,displ);
    ok = model[i]->solutionNorms(gNorm,eNorm,nGauss,displ,asol);
  }
  for (i = 1; i <= gNorm.size(); i++)
    gNorm(i) = sqrt(gNorm(i));

  std::cout <<"Energy norm |u^| = a(u^h,u^h)^0.5 : "<< gNorm(1);
  if (asol)
    std::cout <<"\nExact norm  |u|  = a(u,u)^0.5     : "<< gNorm(2)
	      <<"\nExact error a(e,e)^0.5, e=u-u^h   : "<< gNorm(3)
	      <<"\nExact relative error (%) : "<< gNorm(3)/gNorm(2)*100.0;

  std::cout << std::endl;
  return ok;
}


bool LinearEl::modes (int iop, int nev, int ncv, double shift,
		      std::vector<Mode>& solution)
{
  if (nev < 1 || ncv <= nev) return false;
  if (!sam) return false;

  Vector eigVal;
  Matrix eigVec;
  if (nev > sam->getNoEquations()) nev = sam->getNoEquations();
  if (ncv > sam->getNoEquations()) ncv = sam->getNoEquations();

  // If iop < 0, try the LAPack eigensolver first (for dense matrices only).
  // If iop > 0, or the LAPack solver failed or could not be used, use ARPack.
  if (iop > 0 || !eig::solve(sys.K,sys.M,eigVal,eigVec,nev))
    if (!eig::solve(sys.K,sys.M,eigVal,eigVec,nev,ncv,abs(iop),shift))
      return false;

  bool freq = abs(iop) == 3 || abs(iop) == 4 || abs(iop) == 6;

  solution.resize(nev);
  for (int i = 1; i <= nev; i++)
    if (!sam->expandVector(eigVec.getColumn(i),solution[i-1].eigVec))
      return false;
    else if (!freq)
      solution[i-1].eigVal = eigVal(i);
    else if (eigVal(i) < 0.0)
      solution[i-1].eigVal = -sqrt(-eigVal(i))*0.5/M_PI;
    else
      solution[i-1].eigVal = sqrt(eigVal(i))*0.5/M_PI;

  std::cout <<"\n >>> Computed Eigenvalues <<<\n     Mode\t"
	    << (freq ? "Frequency [Hz]" : "Eigenvalue");
  for (int j = 1; j <= nev; j++)
    std::cout <<"\n     "<< j <<"\t\t"<< solution[j-1].eigVal;
  std::cout << std::endl;

  return true;
}


bool LinearEl::writeGlobalGrid (const char* inputFile, const int* n,
				int nenod) const
{
  ElementBlock globalGrid(nenod), blockGrid(nenod);

  int inod, node, j, ngnod, nlnod;
  IntVec nodeNums, firstEl(model.size()+1,1);
  IntMat FEMbc(3);
  for (size_t i = 0; i < model.size(); i++)
  {
    bool ok = true;
    if (i == 0)
      ok = model[i]->convertToElementBlock(globalGrid,n);
    else
    {
      ok = model[i]->convertToElementBlock(blockGrid,n);
      if (ok) globalGrid.merge(&blockGrid,nodeNums);
    }
    ngnod = globalGrid.getNoNodes();
    nlnod = i == 0 ? ngnod : blockGrid.getNoNodes();
    firstEl[i+1] += globalGrid.getNoElms();

    Matrix field, bc(3,model[i]->getNoNodes());
    std::vector<VolumePatch::BC>::const_iterator bit;
    for (bit = model[i]->begin_BC(); bit != model[i]->end_BC(); bit++)
      if (node = model[i]->getNodeIndex(bit->node))
      {
	if (!bit->CX) bc(1,node) = 1.0;
	if (!bit->CY) bc(2,node) = 1.0;
	if (!bit->CZ) bc(3,node) = 1.0;
      }

    if (ok && !model[i]->evalDisplField(field,bc,n))
      ok = false;

    for (j = 0; j < 3 && ok; j++)
      FEMbc[j].resize(ngnod,0);

    for (inod = 0; inod < nlnod && ok; inod++)
    {
      node = i == 0 ? inod : nodeNums[inod];
      for (j = 0; j < 3; j++)
	if (field(j+1,inod+1) == 1.0)
	  FEMbc[j][node] = 1;
    }
  }

  if (ngnod < nenod) return false;

  char* gridName = new char[strlen(inputFile)+5];
  strcpy(gridName,inputFile);
  std::ofstream os(strcat(strtok(gridName,"."),".grid"));
  std::cout <<"\nWriting global grid file "<< gridName << std::endl;
  delete[] gridName;

  os << ngnod <<" "<< globalGrid.getNoElms() <<" "<< nenod <<"\n";
  inod = 0;
  std::vector<Vec3>::const_iterator cit;
  for (cit = globalGrid.begin_XYZ(); cit != globalGrid.end_XYZ(); cit++, inod++)
    os << inod+1 <<" "<< *cit <<" "<< FEMbc[0][inod] <<" "
       << FEMbc[1][inod] <<" "<< FEMbc[2][inod] <<"\n";

  double E, nu, rho;
  size_t ivolp = 0;
  const int* p = globalGrid.getElements();
  for (int iel = 1; iel <= globalGrid.getNoElms(); iel++)
  {
    while (ivolp < model.size() && iel > firstEl[ivolp+1]) ivolp++;
    model[ivolp]->getMaterial(E,nu,rho);
    os << iel <<" "<< E <<" "<< nu <<" "<< rho;
    for (j = 0; j < nenod; j++, p++)
      os <<" "<< 1+(*p);
    os <<"\n";
  }

  return true;
}


bool LinearEl::writeGlv (const char* inputFile, const Result& results,
			 const int* nViz, int format, bool debug) const
{
  // Open VTF-file
  char* vtfName = new char[strlen(inputFile)+4];
  strcpy(vtfName,inputFile);
  VTF vtf(strcat(strtok(vtfName,"."),".vtf"),format);
  std::cout <<"\nWriting VTF-file "<< vtfName << std::endl;
  delete[] vtfName;

  // Write geometry
  int iStep = 1, nBlock = 0;
  if (!this->writeGlvG(vtf,nViz))
    return false;

  // Write boundary tractions, if any
  if (!this->writeGlvT(vtf,iStep,nBlock))
    return false;

  // Write Dirichlet boundary conditions
  if (!this->writeGlvBC(vtf,nViz,nBlock,debug))
    return false;

  // Write load vector
  if (!this->writeGlvR(vtf,results.load,nViz,iStep,nBlock))
    return false;

  // Write solution fields
  if (!this->writeGlvS(vtf,results.displ,nViz,iStep,nBlock))
    return false;

  // Write eigenmodes
  for (size_t i = 0; i < results.modes.size(); i++, iStep++)
    if (!this->writeGlvM(vtf,results.modes[i],results.freq,nViz,iStep,nBlock))
      return false;

  // Write element norms (only when no additional vizualization points are used)
  if (nViz[0] == 2 && nViz[1] == 2 && nViz[2] == 2)
    if (!this->writeGlvN(vtf,results.norms,iStep,nBlock))
      return false;

  return true;
}


bool LinearEl::writeGlvG (VTF& vtf, const int* nViz) const
{
  char pname[16];
  for (size_t i = 0; i < model.size(); i++)
  {
    if (model[i]->empty()) continue; // skip empty patches

    std::cout <<"Writing geometry for patch "<< i+1 << std::endl;
    ElementBlock* lvb = new ElementBlock;
    if (!model[i]->convertToElementBlock(*lvb,nViz))
      return false;

    sprintf(pname,"Patch %ld",i+1);
    if (!vtf.writeGrid(lvb,pname))
      return false;
  }

  return true;
}


bool LinearEl::writeGlvT (VTF& vtf, int iStep, int& nBlock) const
{
  if (trac.empty())
    return true;

  // Write boundary tractions as discrete point vectors
  std::cout <<"Writing boundary tractions" << std::endl;
  if (!vtf.writeVectors(trac,++nBlock))
    return false;

  return vtf.writeVblk(nBlock,"Tractions",1,iStep);
}


bool LinearEl::writeGlvBC (VTF& vtf, const int* nViz,
			   int& nBlock, bool debug) const
{
  size_t i, j;
  int node, nDupl = 0;
  std::vector<int> dID[6], nodeConn, nodeDupl;
  if (debug && model.size() > 1)
  {
    // Generate some additional data for topology debugging
    std::map<Vec3,int> globalNodes;
    nodeConn.resize(1+sam->getNoNodes(),0);
    nodeDupl.resize(1+sam->getNoNodes(),0);
    for (i = 0; i < model.size(); i++)
      for (node = 1; node <= model[i]->getNoNodes(); node++)
      {
	int globNum = model[i]->getNodeID(node);
	nodeConn[globNum]++;
	Vec3 X(model[i]->getCoord(node));
	std::map<Vec3,int>::const_iterator xit = globalNodes.find(X);
	if (xit == globalNodes.end())
	  globalNodes.insert(std::make_pair(X,globNum));
	else if (xit->second != globNum && nodeDupl[globNum] == 0)
	{
	  std::cout <<"  ** Duplicated nodes "<< xit->second <<" and "
		    << globNum <<" at X="<< X << std::endl;
	  nodeDupl[xit->second]++;
	  nodeDupl[globNum] = -xit->second;
	  nDupl++;
	}
      }
    for (i = 1; i < nodeDupl.size(); i++)
      if (nodeDupl[i] < 0) nodeDupl[i] = nodeDupl[-nodeDupl[i]];
    if (nDupl > 0)
      std::cout <<" *** "<< nDupl <<" duplicated nodes detected."<< std::endl;
  }

  Matrix field;
  int geomID = 0;
  for (i = 0; i < model.size(); i++)
  {
    if (model[i]->empty()) continue; // skip empty patches

    geomID++;
    RealArray flag(debug ? 6 : 3, 0.0);
    Matrix bc(debug ? 6 : 3, model[i]->getNoNodes());
    std::vector<VolumePatch::BC>::const_iterator bit;
    for (bit = model[i]->begin_BC(); bit != model[i]->end_BC(); bit++)
      if (node = model[i]->getNodeIndex(bit->node))
      {
	if (!bit->CX) bc(1,node) = flag[0] = 1.0;
	if (!bit->CY) bc(2,node) = flag[1] = 1.0;
	if (!bit->CZ) bc(3,node) = flag[2] = 1.0;
      }

    if (debug)
    {
      for (node = 1; node <= model[i]->getNoNodes(); node++)
      {
	bc(6,node) = model[i]->getNodeID(node);
	if (model.size() > 1)
	{
	  bc(4,node) = nodeConn[model[i]->getNodeID(node)];
	  bc(5,node) = nodeDupl[model[i]->getNodeID(node)];
	  if (bc(5,node) > 0.0) flag[4] = 1.0;
	}
      }
      if (model.size() > 1) flag[3] = 1.0;
      flag[5] = 1.0;
    }
    else if (flag[0]+flag[1]+flag[2] == 0.0)
      continue; // nothing on this patch

    std::cout <<"Writing boundary conditions for patch "<< i+1 << std::endl;
    if (!model[i]->evalDisplField(field,bc,nViz))
      return false;

    // The BC fields should either be 0.0 or 1.0
    if (nViz[0] > 2 || nViz[1] > 2 || nViz[2] > 2)
      for (j = 1; j <= 3; j++)
	if (flag[j-1] == 1.0)
	  for (int n = 1; n <= field.cols(); n++)
	    if (field(j,n) < 0.9999) field(j,n) = 0.0;

    for (j = 0; j < field.rows(); j++)
      if (flag[j] == 1.0)
	if (vtf.writeNres(field.getRow(1+j),++nBlock,geomID))
	  dID[j].push_back(nBlock);
	else
	  return false;
  }

  const char* label[6] = {
    "fix_X", "fix_Y", "fix_Z",
    "Connectivity",
    "Double nodes",
    "Node number"
  };

  for (j = 0; j < 6; j++)
    if (!dID[j].empty())
      if (!vtf.writeSblk(dID[j],label[j],1+j))
	return false;

  return true;
}


bool LinearEl::writeGlvR (VTF& vtf, const Vector& load,
                          const int* nViz, int iStep, int& nBlock) const
{
  if (load.empty())
    return true;

  Vector lovec;
  Matrix field;
  int geomID = 0, idBlock = 2;
  std::vector<int> vID;
  for (size_t i = 0; i < model.size(); i++)
  {
    if (model[i]->empty()) continue; // skip empty patches

    std::cout <<"Writing load vector for patch "<< i+1 << std::endl;
    model[i]->extractSolution(load,lovec);
    if (!model[i]->evalDisplField(field,lovec,nViz))
      return false;

    if (!vtf.writeVres(field,++nBlock,++geomID))
      return false;
    else
      vID.push_back(nBlock);
  }

  return vtf.writeVblk(vID,"Load vector",idBlock,iStep);
}


bool LinearEl::writeGlvS (VTF& vtf, const Vector& solution,
			  const int* nViz, int iStep, int& nBlock) const
{
  if (solution.empty())
    return true;

  Vector displ;
  Matrix field;
  int j, geomID = 0, idBlock = 10;
  std::vector<int> vID, dID[3], sID[14];
  for (size_t i = 0; i < model.size(); i++)
  {
    if (model[i]->empty()) continue; // skip empty patches

    std::cout <<"Writing FE solution for patch "<< i+1 << std::endl;
    model[i]->extractSolution(solution,displ);
    if (!model[i]->evalDisplField(field,displ,nViz))
      return false;

    if (!vtf.writeVres(field,++nBlock,++geomID))
      return false;
    else
      vID.push_back(nBlock);

    for (j = 0; j < 3; j++)
      if (!vtf.writeNres(field.getRow(1+j),++nBlock,geomID))
	return false;
      else
	dID[j].push_back(nBlock);

    if (!model[i]->evalStressField(field,displ,nViz,rCS))
      return false;

    for (j = 0; j < 6; j++)
      if (!vtf.writeNres(field.getRow(1+j),++nBlock,geomID))
	return false;
      else
	sID[j].push_back(nBlock);

    VolumePatch::vonMises(field,displ);
    if (!vtf.writeNres(displ,++nBlock,geomID))
      return false;
    else
      sID[6].push_back(nBlock);

    if (asol)
    {
      std::cout <<"Writing exact solution for patch "<< i+1 << std::endl;
      std::vector<Vec3>::const_iterator cit;
      const ElementBlock* grid = vtf.getBlock(geomID);
      for (j = 0, cit = grid->begin_XYZ(); cit != grid->end_XYZ(); cit++)
	field.fillColumn(++j,(*asol)(*cit));

      for (j = 0; j < 6; j++)
	if (!vtf.writeNres(field.getRow(1+j),++nBlock,geomID))
	  return false;
	else
	  sID[7+j].push_back(nBlock);

      VolumePatch::vonMises(field,displ);
      if (!vtf.writeNres(displ,++nBlock,geomID))
	return false;
      else
	sID[13].push_back(nBlock);
    }
  }

  if (!vtf.writeDblk(vID,"Displacement",idBlock,iStep))
    return false;

  std::string label("u_x");
  for (j = 0; j < 3; j++)
    if (!vtf.writeSblk(dID[j],label.c_str(),++idBlock,iStep))
      return false;
    else
      label[2]++;

  for (int jj = 0; jj <= (asol ? 7 : 0); jj += 7)
  {
    label = "s_xx";
    std::string pfx(asol ? (jj ? "Exact " : "FE ") : "");
    for (j = 0; j < 7; j++)
      if (!vtf.writeSblk(sID[jj+j],(pfx+label).c_str(),++idBlock,iStep))
        return false;
      else if (j < 5)
      {
        label[2]++;
        label[3]++;
        if (j == 2 || j == 3) label[2] = 'x';
        if (j == 2) label[3] = 'y';
        if (j == 4) label[3] = 'z';
      }
      else
        label = "von Mises stress";
  }

  return vtf.writeState(iStep,"Step %g",(double)iStep,2);
}


bool LinearEl::writeGlvM (VTF& vtf, const Mode& mode, bool freq,
			  const int* nViz, int iStep, int& nBlock) const
{
  if (mode.eigVec.empty())
    return true;

  std::cout <<"Writing eigenvector for Mode "<< iStep << std::endl;

  Vector displ;
  Matrix field;
  int geomID = 0, idBlock = 10;
  std::vector<int> vID;

  for (size_t i = 0; i < model.size(); i++)
  {
    if (model[i]->empty()) continue; // skip empty patches
    if (model.size() > 1) std::cout <<"."<< std::flush;

    geomID++;
    model[i]->extractSolution(mode.eigVec,displ);
    if (!model[i]->evalDisplField(field,displ,nViz))
      return false;

    if (!vtf.writeVres(field,++nBlock,geomID))
      return false;
    else
      vID.push_back(nBlock);
  }
  if (model.size() > 1) std::cout << std::endl;

  if (!vtf.writeDblk(vID,"Mode Shape",idBlock,iStep))
    return false;

  return vtf.writeState(iStep, freq ? "Frequency %g" : "Eigenvalue %g",
			mode.eigVal, 1);
}


bool LinearEl::writeGlvN (VTF& vtf, const Matrix& norms,
                          int iStep, int& nBlock) const
{
  if (norms.empty())
    return true;
  else if (norms.rows() > 3)
    return false;

  Matrix field;
  int j, idBlock = 100, geomID = 0;
  std::vector<int> sID[3];

  for (size_t i = 0; i < model.size(); i++)
  {
    if (model[i]->empty()) continue; // skip empty patches

    geomID++;
    std::cout <<"Writing element norms for patch "<< i+1 << std::endl;
    model[i]->extractElmRes(norms,field);

    for (j = 0; j < field.rows(); j++)
      if (!vtf.writeEres(field.getRow(1+j),++nBlock,geomID))
        return false;
      else
        sID[j].push_back(nBlock);
  }

  const char* label[3] = {
    "a(u^h,u^h)^0.5",
    "a(u,u)^0.5",
    "a(e,e)^0.5, e=u-u^h"
  };

  for (j = 0; j < norms.rows(); j++)
    if (!vtf.writeSblk(sID[j],label[j],++idBlock,iStep))
      return false;

  return true;
}


void LinearEl::dumpGeometry (const char* g2file) const
{
  std::cout <<"\nWriting updated g2-file "<< g2file << std::endl;

  std::ofstream os(g2file);
  for (size_t i = 0; i < model.size() && os; i++)
    model[i]->write(os);
}


void LinearEl::dumpSolution (const char* solfile, const Vector& solution) const
{
  std::cout <<"\nWriting solution files "<< solfile << std::endl;

  char* fileName = new char[strlen(solfile)+4];
  std::ofstream osv(strcat(strcpy(fileName,solfile),".vec"));
  std::ofstream os1(strcat(strtok(fileName,"."),".u"));
  std::ofstream os2(strcat(strtok(fileName,"."),".v"));
  std::ofstream os3(strcat(strtok(fileName,"."),".w"));
  utl::nval_per_line = 3;

  Vector displ;
  for (size_t i = 0; i < model.size(); i++)
  {
    model[i]->extractSolution(solution,displ);
    if (displ.empty()) continue;

    osv << displ;
    for (size_t j = 0; j < displ.size();)
    {
      os1 << displ[j++] <<"\n";
      os2 << displ[j++] <<"\n";
      os3 << displ[j++] <<"\n";
    }
  }

  delete[] fileName;
}


void LinearEl::dumpMat (const char* Kfile, const char* Mfile,
			const char* Rfile) const
{
  if (Kfile && sys.K)
  {
    std::cout <<"\nWriting "<< Kfile << std::endl;
    std::ofstream os(Kfile);
    if (os) os << *sys.K;
  }

  if (Mfile && sys.M)
  {
    std::cout <<"\nWriting "<< Mfile << std::endl;
    std::ofstream os(Mfile);
    if (os) os << *sys.M;
  }

  if (Rfile && !sys.RHS.empty())
  {
    std::cout <<"\nWriting "<< Rfile << std::endl;
    std::ofstream os(Rfile);
    if (os) os << sys.RHS;
  }
}
