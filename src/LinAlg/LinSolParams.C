// $Id$
//==============================================================================
//!
//! \file LinSolParams.C
//!
//! \date Jan 29 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Linear solver parameters for PETSc matrices.
//!
//==============================================================================

#include "LinSolParams.h"
#include "PCPerm.h"
#include "Utilities.h"
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <utility>
#include <iterator>
#include "tinyxml.h"


void LinSolParams::setDefault ()
{
#ifdef HAS_PETSC
  // Use GMRES with ILU preconditioner as default
  method      = KSPGMRES;
  prec        = PCILU;
  subprec.resize(1,PCILU);
  hypretype.resize(1,"boomeramg");
  package.resize(1,MATSOLVERPETSC);
  levels.resize(1,0);
  mglevels.resize(1,1);
  presmoother.resize(1,PCILU);
  noPreSmooth.resize(1,1);
  postsmoother.resize(1,PCILU);
  noPostSmooth.resize(1,1);
  noFineSmooth.resize(1,1);
  maxCoarseSize.resize(1,1000);
  overlap.resize(1,1);
  nx.resize(1,0);
  ny.resize(1,0);
  nz.resize(1,0);
  nullspc.resize(1,NONE);
  asmlu.resize(1,false);
  nblock   = 1;
  schur    = false;
  schurPrec = SIMPLE;
  ncomps.resize(2); ncomps[0]   = 2; ncomps[1]   = 1;
    
  atol      = 1.0e-6;
  rtol      = 1.0e-6;
  dtol      = 1.0e3;
  maxIts    = 1000;
#endif
}


void LinSolParams::copy (const LinSolParams& spar)
{
#ifdef HAS_PETSC
  // Copy linear solver parameters
  method        = spar.method;
  prec          = spar.prec;
  subprec       = spar.subprec;
  hypretype     = spar.hypretype;
  package       = spar.package;
  overlap       = spar.overlap;
  levels        = spar.levels;
  mglevels      = spar.mglevels;
  presmoother   = spar.presmoother;
  postsmoother  = spar.postsmoother;
  noPreSmooth   = spar.noPreSmooth;
  noPostSmooth  = spar.noPostSmooth;
  dirsmoother   = spar.dirsmoother;
  dirOrder      = spar.dirOrder;
  maxCoarseSize = spar.maxCoarseSize;
  overlap       = spar.overlap;
  nx            = spar.nx;
  ny            = spar.ny;
  nz            = spar.nz;
  nullspc       = spar.nullspc;
  asmlu         = spar.asmlu;
  nblock        = spar.nblock;
  schur         = spar.schur;
  schurPrec     = spar.schurPrec;
  ncomps        = spar.ncomps;

  atol   = spar.atol;
  rtol   = spar.rtol;
  dtol   = spar.dtol;
  maxIts = spar.maxIts;
#endif
}


bool LinSolParams::read (std::istream& is, int nparam)
{
  char* cline = 0;
  for (int i = 0; i < nparam && (cline = utl::readLine(is)); i++)
#ifdef HAS_PETSC
    if (!strncasecmp(cline,"type",4)) {
      char* c = strchr(cline,'=');
      for (++c; *c == ' '; c++);
      int j = 0;
      while (c[j] != EOF && c[j] != '\n' && c[j] != ' ' && c[j] != '\0') j++;
      method.assign(c,j);
    }

    else if (!strncasecmp(cline,"pc",2)) {
      char* c = strchr(cline,'=');
      for (++c; *c == ' '; c++);
      int j = 0;
      while (c[j] != EOF && c[j] != '\n' && c[j] != ' ' && c[j] != '\0') j++;
      prec.assign(c,j);
      if ((asmlu[0] = (prec.substr(0,5) == "asmlu")))
	prec = "asm";
    }

    else if (!strncasecmp(cline,"hypretype",9)) {
      char* c = strchr(cline,'=');
      for (++c; *c == ' '; c++);
      int j = 0;
      while (c[j] != EOF && c[j] != '\n' && c[j] != ' ' && c[j] != '\0') j++;
      hypretype[0].assign(c,j);
    }

    else if (!strncasecmp(cline,"package",7)) {
      char* c = strchr(cline,'=');
      for (++c; *c == ' '; c++);
      int j = 0;
      while (c[j] != EOF && c[j] != '\n' && c[j] != ' ' && c[j] != '\0') j++;
      package[0].assign(c,j);
    }

    else if (!strncasecmp(cline,"levels",6)) {
      char* c = strchr(cline,'=');
      levels[0] = atoi(++c);
    }

    else if (!strncasecmp(cline,"overlap",7)) {
      char* c = strchr(cline,'=');
      overlap[0] = atoi(++c);
    }

    else if (!strncasecmp(cline,"atol",4)) {
      char* c = strchr(cline,'=');
      atol = atof(++c);
    }

    else if (!strncasecmp(cline,"rtol",4)) {
      char* c = strchr(cline,'=');
      rtol = atof(++c);
    }

    else if (!strncasecmp(cline,"dtol",4)) {
      char* c = strchr(cline,'=');
      dtol = atof(++c);
    }

    else if (!strncasecmp(cline,"maxits",6)) {
      char* c = strchr(cline,'=');
      maxIts = atoi(++c);
    }

    else if (!strncasecmp(cline,"nx",2)) {
      char* c = strchr(cline,'=');
      nx[0] = atoi(++c);
    }

    else if (!strncasecmp(cline,"ny",2)) {
      char* c = strchr(cline,'=');
      ny[0] = atoi(++c);
    }

    else if (!strncasecmp(cline,"nz",2)) {
      char* c = strchr(cline,'=');
      nz[0] = atoi(++c);
    }

    else if (!strncasecmp(cline,"ncomponents",11)) {
      char* c = strchr(cline,'=');
      std::istringstream this_line(c);
      std::istream_iterator<int> begin(this_line), end;
      ncomps.assign(begin, end);
      nblock = ncomps.size();
    }

    else if (!strncasecmp(cline,"nullspace",6)) {
      char* c = strchr(cline,'=');
      if (!strncasecmp(c,"CONSTANT",8))
	nullspc[0] = CONSTANT;
      else if (!strncasecmp(c,"RIGID_BODY",10))
	nullspc[0] = RIGID_BODY;
      else
	nullspc[0] = NONE;
    }
    else {
      std::cerr <<" *** LinSolParams::read: Unknown keyword: "
		<< cline << std::endl;
      return false;
    }
#else
    ;
#endif

  return true;
}


bool LinSolParams::read (const TiXmlElement* child)
{
#ifdef HAS_PETSC
  const char* value = 0;
  if ((value = utl::getValue(child,"type")))
    method = value;
  else if ((value = utl::getValue(child,"pc"))) {
    prec = value;
    if (!strncasecmp(prec.c_str(),"asmlu",5)) {
      prec = "asm";
      asmlu[0] = true;
    }
  }
  else if ((value = utl::getValue(child,"schurpc"))) {
    if (!strncasecmp(value,"SIMPLE",6))
      schurPrec = SIMPLE;
    else if (!strncasecmp(value,"MSIMPLER",8))
      schurPrec = MSIMPLER;
    else if (!strncasecmp(value,"PCD",3))
      schurPrec = PCD;
  }
  else if ((value = utl::getValue(child,"subpc"))) {
    std::istringstream this_line(value);
    std::istream_iterator<std::string> begin(this_line), end;
    subprec.assign(begin,end);
    asmlu.resize(subprec.size());
    for (size_t i = 0;i < subprec.size();i++)
      if (!strncasecmp(subprec[i].c_str(),"asmlu",5)) {
	asmlu[i] = true;
	subprec[i] = "asm";
      }
  }
  else if ((value = utl::getValue(child,"hypretype"))) {
    std::istringstream this_line(value);
    std::istream_iterator<std::string> begin(this_line), end;
    hypretype.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"package"))) {
    std::istringstream this_line(value);
    std::istream_iterator<std::string> begin(this_line), end;
    package.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"levels"))) {
    std::istringstream this_line(value);
    std::istream_iterator<int> begin(this_line), end;
    levels.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"mglevels"))) {
    std::istringstream this_line(value);
    std::istream_iterator<int> begin(this_line), end;
    mglevels.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"noPreSmooth"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscInt> begin(this_line), end;
    noPreSmooth.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"noPostSmooth"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscInt> begin(this_line), end;
    noPostSmooth.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"noFineSmooth"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscInt> begin(this_line), end;
    noFineSmooth.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"presmoother"))) {
    std::istringstream this_line(value);
    std::istream_iterator<std::string> begin(this_line), end;
    presmoother.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"postsmoother"))) {
    std::istringstream this_line(value);
    std::istream_iterator<std::string> begin(this_line), end;
    postsmoother.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"finesmoother"))) {
    std::istringstream this_line(value);
    std::istream_iterator<std::string> begin(this_line), end;
    finesmoother.assign(begin, end);
  }
  else if (!strcasecmp(child->Value(),"dirsmoother")) {
    size_t block, order;
    std::string type;

    if (!utl::getAttribute(child,"block",block))
      return false;
    if (!utl::getAttribute(child,"type",type))
      return false;
    if (!utl::getAttribute(child,"order",order))
      return false;

    if (block < 1)
      return false;

    if (dirsmoother.empty() || (dirsmoother.size() < block))
      dirsmoother.resize(std::max(subprec.size(),block));
    dirsmoother[block-1].push_back(type);

    if (dirOrder.empty() || (dirOrder.size() < block))
      dirOrder.resize(std::max(subprec.size(),block));
    dirOrder[block-1].push_back(order);
  }
  else if ((value = utl::getValue(child,"overlap"))) {
    std::istringstream this_line(value);
    std::istream_iterator<int> begin(this_line), end;
    overlap.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"nx"))) {
    std::istringstream this_line(value);
    std::istream_iterator<int> begin(this_line), end;
    nx.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"ny"))) {
    std::istringstream this_line(value);
    std::istream_iterator<int> begin(this_line), end;
    ny.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"nz"))) {
    std::istringstream this_line(value);
    std::istream_iterator<int> begin(this_line), end;
    nz.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"atol")))
    atol = atof(value);
  else if ((value = utl::getValue(child,"rtol")))
    rtol = atof(value);
  else if ((value = utl::getValue(child,"dtol")))
    dtol = atof(value);
  else if ((value = utl::getValue(child,"maxits")))
    maxIts = atoi(value);
  else if ((value = utl::getValue(child,"maxCoarseSize"))) {
    std::istringstream this_line(value);
    std::istream_iterator<int> begin(this_line), end;
    maxCoarseSize.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"MLCoarsenScheme"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscInt> begin(this_line), end;
    MLCoarsenScheme.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"MLThreshold"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscReal> begin(this_line), end;
    MLThreshold.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"MLDampingFactor"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscReal> begin(this_line), end;
    MLDampingFactor.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"ncomponents"))) {
    std::istringstream this_line(value);
    std::istream_iterator<int> begin(this_line), end;
    ncomps.assign(begin, end);
    nblock = ncomps.size();
  }
  else if ((value = utl::getValue(child,"nullspace"))) {
    std::istringstream this_line(value);
    std::istream_iterator<std::string> begin(this_line), end;
    std::istream_iterator<std::string> it;
    for (it = begin;it != end;it++) {
      if (!strcasecmp(it->c_str(),"constant"))
	nullspc.push_back(CONSTANT);
      else if (!strcasecmp(value,"rigid_body"))
	nullspc.push_back(RIGID_BODY);
    }
  }
  else
  {
    std::cerr <<" *** LinSolParams::read: Unknown keyword: "
	      << child->Value() << std::endl;
    return false;
  }
#endif
  return true;
}


bool LinSolParams::read (const char* filename)
{
  std::ifstream fs(filename);
  return this->read(fs,1000);
}


#ifdef HAS_PETSC
void LinSolParams::setParams (KSP& ksp, PetscIntMat& locSubdDofs,
			      PetscIntMat& subdDofs, PetscRealVec& coords, 
			      PetscInt nsd, ISMat& dirIndexSet) const
{
  // Set linear solver method
  KSPSetType(ksp,method.c_str());
  KSPSetTolerances(ksp,rtol,atol,dtol,maxIts);
  //KSPSetFromOptions(ksp);

  // Set preconditioner
  PC pc;
  KSPGetPC(ksp,&pc);
  if (!strncasecmp(prec.c_str(),"compositedir",12)) {
    Mat mat;
    Mat Pmat;
    MatStructure flag;
    PCGetOperators(pc,&mat,&Pmat,&flag);
    
    this->addDirSmoother(pc,Pmat,dirIndexSet);
  }
  else
    PCSetType(pc,prec.c_str());
  //PCFactorSetMatSolverPackage(pc,package[0].c_str());
  //PCFactorSetLevels(pc,levels[0]);
#if PETSC_HAVE_HYPRE
  if (!strncasecmp(prec.c_str(),"hypre",5))
    PCHYPRESetType(pc,hypretype[0].c_str());
#endif

  if (!strncasecmp(prec.c_str(),"gamg",4)) {
    PetscInt nloc = coords.size()/nsd;
    PCSetCoordinates(pc,nsd,nloc,&coords[0]);
  }

  if (!strncasecmp(prec.c_str(),"asm",3) ||!strncasecmp(prec.c_str(),"gasm",4)) {
    PCASMSetType(pc,PC_ASM_BASIC);
    PCASMSetOverlap(pc,overlap[0]);

    if (!locSubdDofs.empty() && !subdDofs.empty()) {
      const size_t nsubds = subdDofs.size();
      
      IS isLocSubdDofs[nsubds], isSubdDofs[nsubds];
      for (size_t i = 0;i < nsubds;i++) {
	ISCreateGeneral(PETSC_COMM_SELF,locSubdDofs[i].size(),&(locSubdDofs[i][0]),PETSC_USE_POINTER,&(isLocSubdDofs[i]));
	ISCreateGeneral(PETSC_COMM_SELF,subdDofs[i].size(),&(subdDofs[i][0]),PETSC_USE_POINTER,&(isSubdDofs[i]));
      }
      PCASMSetLocalSubdomains(pc,nsubds,isSubdDofs,isLocSubdDofs);
    }

    

    PCSetFromOptions(pc);
    PCSetUp(pc);

    if (asmlu[0]) {
      KSP* subksp;
      PC   subpc;
      PetscInt first, nlocal;
      PCASMGetSubKSP(pc,&nlocal,&first,&subksp);
      
      for (int i = 0; i < nlocal; i++) {
	KSPGetPC(subksp[i],&subpc);
	PCSetType(subpc,PCLU);
	KSPSetType(subksp[i],KSPPREONLY);
      }
    }
  }
  else if (!strncasecmp(prec.c_str(),"ml",2)) {
        PetscInt n;

        //PCMGSetLevels(pc,mglevels[0], PETSC_NULL);
        //PCMGSetNumberSmoothDown(pc,noPreSmooth[0]);
        //PCMGSetNumberSmoothUp(pc,noPostSmooth[0]);
	std::stringstream maxLevel;
	maxLevel << mglevels[0];
        PetscOptionsSetValue("-pc_ml_maxNLevels",maxLevel.str().c_str());
	std::stringstream maxCoarseDof;
	maxCoarseDof << maxCoarseSize[0];
        PetscOptionsSetValue("-pc_ml_maxCoarseSize",maxCoarseDof.str().c_str());
	if (!MLCoarsenScheme.empty()) {
	  std::stringstream coarsenScheme;
	  coarsenScheme << MLCoarsenScheme[0];
	  PetscOptionsSetValue("-pc_ml_CoarsenScheme",coarsenScheme.str().c_str());
	}
	if (!MLThreshold.empty()) {
	  std::stringstream threshold;
	  threshold << MLThreshold[0];
	  PetscOptionsSetValue("-pc_ml_Threshold",threshold.str().c_str());
	}
	if (!MLDampingFactor.empty()) {
	  std::stringstream damping;
	  damping << MLDampingFactor[0];
	  PetscOptionsSetValue("-pc_ml_DampingFactor",damping.str().c_str());
	}
      
        //PCGAMGSetNlevels(pc,mglevels[0]);
        //PCGAMGSetCoarseEqLim(pc,maxCoarseSize[0]);

	PCSetFromOptions(pc);
	PCSetUp(pc);

        PCMGGetLevels(pc,&n);
	// Presmoother settings
	for (int i = 1;i < n;i++) {
          KSP preksp;
          PC  prepc;
          // Not working for some reason
          //PCMGGetSmootherDown(pc,i,&preksp);

	  // Set smoother
	  std::string smoother;
	  PetscInt noSmooth;
	  
	  PCMGGetSmoother(pc,i,&preksp);
	  KSPSetType(preksp,"richardson");
	  
	  if ((i == n-1) && (!finesmoother.empty())) {
            smoother = finesmoother[0];
	    noSmooth = noFineSmooth[0];
	  }
	  else {
	    smoother = presmoother[0];
	    noSmooth = noPreSmooth[0];
	  }
	  
	  KSPSetTolerances(preksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,noSmooth);
	  KSPGetPC(preksp,&prepc);

	  if (smoother == "asm" || smoother == "asmlu") {
	    PCSetType(prepc,"asm");
	    PCASMSetType(prepc,PC_ASM_BASIC);
	    PCASMSetOverlap(prepc,overlap[0]);
	    
	    if (!locSubdDofs.empty() && !subdDofs.empty() && (i==n-1)) {
	      const size_t nsubds = subdDofs.size();

	      IS isLocSubdDofs[nsubds], isSubdDofs[nsubds];
	      for (size_t k = 0;k < nsubds;k++) {
		ISCreateGeneral(PETSC_COMM_SELF,locSubdDofs[k].size(),&(locSubdDofs[k][0]),PETSC_USE_POINTER,&(isLocSubdDofs[k]));
		ISCreateGeneral(PETSC_COMM_SELF,subdDofs[k].size(),&(subdDofs[k][0]),PETSC_USE_POINTER,&(isSubdDofs[k]));
	      }
	      PCASMSetLocalSubdomains(prepc,nsubds,isSubdDofs,isLocSubdDofs);
	    }

	    // If LU factorization is used on each subdomain
	    if (smoother == "asmlu") {
	      KSP* subksp;
	      PC   subpc;
	      PetscInt first, nlocal;
	      PCSetFromOptions(prepc);
	      PCSetUp(prepc);
	      PCASMGetSubKSP(prepc,&nlocal,&first,&subksp);
	      
	      for (int k = 0; k < nlocal; k++) {
		KSPGetPC(subksp[k],&subpc);
		PCSetType(subpc,PCLU);
		PCFactorSetShiftType(prepc,MAT_SHIFT_NONZERO);
		KSPSetType(subksp[k],KSPPREONLY);
	      }
	    }
	  }
	  else if (smoother == "compositedir" && (i==n-1)) {
	    Mat mat;
	    Mat Pmat;
	    MatStructure flag;
	    PCGetOperators(prepc,&mat,&Pmat,&flag);
    
	    this->addDirSmoother(prepc,Pmat,dirIndexSet);
	  }
	  else
	    PCSetType(prepc,smoother.c_str());

	  PCFactorSetLevels(prepc,levels[0]); 
          KSPSetUp(preksp);
	}

	// Postsmoother settings
	// for (int i = 1;i < n;i++) {
        //   KSP postksp;
        //   PC  postpc;
        //   // Not working for some reason
        //   //PCMGGetSmootherUp(pc,i,&postksp);
        //   PCMGGetSmoother(pc,i,&postksp);
        //   KSPSetType(postksp,"richardson");
        //   KSPSetTolerances(postksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,noPostSmooth[0]);
        //   KSPGetPC(postksp,&postpc);

	//   // Set smoother
	//   std::string smoother;
        //   if ((i == n-1) && (!finesmoother.empty())) 	    
        //     smoother = finesmoother[0];
        //   else
	//     smoother = postsmoother[0];

	//   if (smoother == "asm" || smoother == "asmlu") {
	//     PCSetType(postpc,"asm");
	//     PCASMSetType(postpc,PC_ASM_BASIC);
	//     PCASMSetOverlap(postpc,overlap[0]);
	    
	//     if (!locSubdDofs.empty() && !subdDofs.empty()) {
	//       const size_t nsubds = subdDofs.size();

	//       IS isLocSubdDofs[nsubds], isSubdDofs[nsubds];
	//       for (size_t i = 0;i < nsubds;i++) {
	// 	ISCreateGeneral(PETSC_COMM_WORLD,locSubdDofs[i].size(),&(locSubdDofs[i][0]),PETSC_USE_POINTER,&(isLocSubdDofs[i]));
	// 	ISCreateGeneral(PETSC_COMM_WORLD,subdDofs[i].size(),&(subdDofs[i][0]),PETSC_USE_POINTER,&(isSubdDofs[i]));
	//       }
	//       PCASMSetLocalSubdomains(postpc,nsubds,isSubdDofs,isLocSubdDofs);
	//     }

	//     // If LU factorization is used on each subdomain
	//     if (smoother == "asmlu") {
	//       KSP* subksp;
	//       PC   subpc;
	//       PetscInt first, nlocal;
	//       PCSetFromOptions(postpc);
	//       PCSetUp(postpc);
	//       PCASMGetSubKSP(postpc,&nlocal,&first,&subksp);
	      
	//       for (int i = 0; i < nlocal; i++) {
	// 	KSPGetPC(subksp[i],&subpc);
	// 	PCSetType(subpc,PCLU);
	// 	PCFactorSetShiftType(postpc,MAT_SHIFT_NONZERO);
	// 	KSPSetType(subksp[i],KSPPREONLY);
	//       }
	//     }
	//   }
	//   else
	//     PCSetType(postpc,smoother.c_str());

	//   PCFactorSetLevels(postpc,levels[0]); 
        //   KSPSetUp(postksp);
	// }
  }
  else {
    PCSetFromOptions(pc);
    PCSetUp(pc);
  }
  
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);
}

bool LinSolParams::addDirSmoother(PC pc, Mat P, ISMat& dirIndexSet) const
{
  PCSetType(pc,"composite");
  PCCompositeSetType(pc,PC_COMPOSITE_MULTIPLICATIVE);
  for (size_t k = 0;k < dirsmoother[0].size();k++)
    PCCompositeAddPC(pc,"shell");
  for (size_t k = 0;k < dirsmoother[0].size();k++) {
    PC dirpc;
    PCCompositeGetPC(pc,k,&dirpc);
    PCPerm *pcperm;
    PCPermCreate(&pcperm);
    PCShellSetApply(dirpc,PCPermApply);
    PCShellSetContext(dirpc,pcperm);
    PCShellSetDestroy(dirpc,PCPermDestroy);
    PCShellSetName(dirpc,"dir");
    PCPermSetUp(dirpc,&dirIndexSet[0][k],P,dirsmoother[0][k].c_str());
   }
  
  return true;
}

#endif
