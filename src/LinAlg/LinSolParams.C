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
  nResetSolver = 0;

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
  mgKsp.resize(1,"defrichardson");
  GAMGtype.resize(1,"agg");
  GAMGrepartition.resize(1,PETSC_FALSE);
  GAMGuseAggGasm.resize(1,PETSC_FALSE);
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
  else if (value = utl::getValue(child,"noStepBeforeReset"))
    nResetSolver = atoi(value);
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
  else if ((value = utl::getValue(child,"mgksp"))) {
    std::istringstream this_line(value);
    std::istream_iterator<std::string> begin(this_line), end;
    mgKsp.assign(begin, end);
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
  else if ((value = utl::getValue(child,"GAMGtype"))) {
    std::istringstream this_line(value);
    std::istream_iterator<std::string> begin(this_line), end;
    GAMGtype.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"GAMGrepartition"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscInt> begin(this_line), end;
    GAMGrepartition.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"GAMGuseAggGasm"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscInt> begin(this_line), end;
    GAMGuseAggGasm.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"GAMGprocEqLimit"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscInt> begin(this_line), end;
    GAMGprocEqLimit.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"GAMGreuseInterpolation"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscInt> begin(this_line), end;
    GAMGreuseInterp.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"MLCoarsePackage"))) {
    std::istringstream this_line(value);
    std::istream_iterator<std::string> begin(this_line), end;
    MLCoarsePackage.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"MLCoarseSolver"))) {
    std::istringstream this_line(value);
    std::istream_iterator<std::string> begin(this_line), end;
    MLCoarseSolver.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"MLCoarsenScheme"))) {
    std::istringstream this_line(value);
    std::istream_iterator<std::string> begin(this_line), end;
    MLCoarsenScheme.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"MLSymmetrize"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscInt> begin(this_line), end;
    MLSymmetrize.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"MLRepartition"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscInt> begin(this_line), end;
    MLRepartition.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"MLBlockScaling"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscInt> begin(this_line), end;
    MLBlockScaling.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"MLPutOnSingleProc"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscInt> begin(this_line), end;
    MLPutOnSingleProc.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"MLReuseInterpolation"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscInt> begin(this_line), end;
    MLReuseInterp.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"MLReuseable"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscInt> begin(this_line), end;
    MLReusable.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"MLKeepAggInfo"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscInt> begin(this_line), end;
    MLKeepAggInfo.assign(begin, end);
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
  else if ((value = utl::getValue(child,"MLRepartitionRatio"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscReal> begin(this_line), end;
    MLRepartitionRatio.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"MLAux"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscInt> begin(this_line), end;
    MLAux.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"HypreNoAggCoarse"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscInt> begin(this_line), end;
    HypreNoAggCoarse.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"HypreNoPathAggCoarse"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscInt> begin(this_line), end;
    HypreNoPathAggCoarse.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"HypreTruncation"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscReal> begin(this_line), end;
    HypreTruncation.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"HypreThreshold"))) {
    std::istringstream this_line(value);
    std::istream_iterator<PetscReal> begin(this_line), end;
    HypreThreshold.assign(begin, end);
  }
  else if ((value = utl::getValue(child,"HypreCoarsenScheme"))) {
    std::istringstream this_line(value);
    std::istream_iterator<std::string> begin(this_line), end;
    HypreCoarsenScheme.assign(begin, end);
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
    nullspc.clear();
    for (it = begin;it != end;it++) {
      if (!strcasecmp(it->c_str(),"constant"))
	nullspc.push_back(CONSTANT);
      else if (!strcasecmp(value,"rigid_body"))
	nullspc.push_back(RIGID_BODY);
      else
        nullspc.push_back(NONE);
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

int LinSolParams::getLocalPartitioning(size_t dir, size_t i) const
{
  switch (dir) {
  case 0: return nx[i];
  case 1: return ny[i];
  case 2: return nz[i];
  default: return 0;
  }
}


void LinSolParams::setParams (KSP& ksp, PetscIntMat& locSubdDofs,
			      PetscIntMat& subdDofs, PetscRealVec& coords, 
			      ISMat& dirIndexSet) const
{
  // Set linear solver method
  KSPSetType(ksp,method.c_str());
  KSPSetTolerances(ksp,rtol,atol,dtol,maxIts);

  // Set preconditioner
  PC pc;
  KSPGetPC(ksp,&pc);
  if (!strncasecmp(prec.c_str(),"compositedir",12)) {
    Mat mat;
    Mat Pmat;
#if PETSC_VERSION_MINOR < 5
    MatStructure flag;
    PCGetOperators(pc,&mat,&Pmat,&flag);
#else
    PCGetOperators(pc,&mat,&Pmat);
#endif
    
    this->addDirSmoother(pc,Pmat,0,dirIndexSet);
  }
  else
    PCSetType(pc,prec.c_str());
#if PETSC_HAVE_HYPRE
  if (!strncasecmp(prec.c_str(),"hypre",5)) {
    PCHYPRESetType(pc,hypretype[0].c_str());
    setHypreOptions("", 0);
  }
#endif

  if (!strncasecmp(prec.c_str(),"gamg",4) || !strncasecmp(prec.c_str(),"ml",2) ) {
    PetscInt nloc = coords.size()/nsd;
    PCSetCoordinates(pc,nsd,nloc,&coords[0]);
    PCGAMGSetType(pc, PCGAMGAGG); // TODO?
  }

  if (!strncasecmp(prec.c_str(),"asm",3) ||!strncasecmp(prec.c_str(),"gasm",4))
    setupAdditiveSchwarz(pc, overlap[0], asmlu[0], locSubdDofs, subdDofs, false);
  else if (!strncasecmp(prec.c_str(),"ml",2) || !strncasecmp(prec.c_str(),"gamg",4)) {
    if (!strncasecmp(prec.c_str(),"ml",2))
      setMLOptions("", 0);
    else if (!strncasecmp(prec.c_str(),"gamg",4))
      setGAMGOptions("", 0);
    
    PCSetFromOptions(pc);
    PCSetUp(pc);
    
    // Settings for coarse solver
    if (!MLCoarseSolver.empty())
      setupCoarseSolver(pc, "", 0);

    setupSmoothers(pc, 0, dirIndexSet, locSubdDofs, subdDofs);
  } else {
    PCSetFromOptions(pc);
    PCSetUp(pc);
  }
  
  KSPSetFromOptions(ksp);
  KSPSetUp(ksp);

  MPI_Comm comm;
  PetscObjectGetComm((PetscObject)ksp,&comm);
  PCView(pc,PETSC_VIEWER_STDOUT_(comm));
}


bool LinSolParams::addDirSmoother(PC pc, Mat P, int block,
                                  ISMat& dirIndexSet) const
{
  PCSetType(pc,"composite");
  PCCompositeSetType(pc,PC_COMPOSITE_MULTIPLICATIVE);
  for (size_t k = 0;k < dirsmoother[block].size();k++)
    PCCompositeAddPC(pc,"shell");
  for (size_t k = 0;k < dirsmoother[block].size();k++) {
    PC dirpc;
    PCCompositeGetPC(pc,k,&dirpc);
    PCPerm *pcperm;
    PCPermCreate(&pcperm);
    PCShellSetApply(dirpc,PCPermApply);
    PCShellSetContext(dirpc,pcperm);
    PCShellSetDestroy(dirpc,PCPermDestroy);
    PCShellSetName(dirpc,"dir");
    PCPermSetUp(dirpc,&dirIndexSet[block][k],P,dirsmoother[block][k].c_str());
  }
  
  return true;
}


/*! \brief Static helper to optionally add a prefix to a PETSc parameter */
static std::string AddPrefix(const std::string& prefix, const std::string& data)
{
  if (prefix.empty())
    return "-"+data;

  return "-"+prefix+"_"+data;
}


/*! \brief Static helper to convert numbers to a string */
template<class T>
static std::string ToString(T data)
{
  std::stringstream str;
  str << data;
  return str.str();
}


void LinSolParams::setMLOptions(const std::string& prefix, int block) const
{
  PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_maxNLevels").c_str(),
                       ToString(mglevels[block]).c_str());
  PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_maxCoarseSize").c_str(),
                       ToString(maxCoarseSize[block]).c_str());
  if (!MLCoarsenScheme.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_CoarsenScheme").c_str(),
                         MLCoarsenScheme[block].c_str());
  if (!MLThreshold.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_Threshold").c_str(),
                         ToString(MLThreshold[block]).c_str());
  if (!MLDampingFactor.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_DampingFactor").c_str(),
                         ToString(MLDampingFactor[block]).c_str());
  if (!MLRepartitionRatio.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_repartitionMaxMinRatio").c_str(),
                         ToString(MLRepartitionRatio[block]).c_str());
  if (!MLSymmetrize.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_Symmetrize").c_str(),
                         ToString(MLSymmetrize[block]).c_str());
  if (!MLRepartition.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_repartition").c_str(),
                         ToString(MLRepartition[block]).c_str());
  if (!MLBlockScaling.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_BlockScaling").c_str(),
                         ToString(MLBlockScaling[block]).c_str());
  if (!MLPutOnSingleProc.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_repartitionPutOnSingleProc").c_str(),
                         ToString(MLPutOnSingleProc[block]).c_str());
  if (!MLReuseInterp.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_reuse_interpolation").c_str(),
                         ToString(MLReuseInterp[block]).c_str());
  if (!MLKeepAggInfo.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_KeepAggInfo").c_str(),
                         ToString(MLKeepAggInfo[block]).c_str());
  if (!MLReusable.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_Reusable").c_str(),
                         ToString(MLReusable[block]).c_str());
  if (!MLAux.empty()) {
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_Aux").c_str(),
                         ToString(MLAux[block]).c_str());
    if (!MLThreshold.empty())
      PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_AuxThreshold").c_str(),
                           ToString(MLThreshold[block]).c_str());
  }
}


void LinSolParams::setGAMGOptions(const std::string& prefix, int block) const
{
  if (!maxCoarseSize.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_gamg_coarse_eq_limit").c_str(),
                         ToString(maxCoarseSize[block]).c_str());
  if (!GAMGprocEqLimit.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_gamg_process_eq_limit").c_str(),
                         ToString(GAMGprocEqLimit[block]).c_str());
  if (!GAMGtype.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_gamg_type").c_str(),
                         GAMGtype[block].c_str());
  if (!GAMGrepartition.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_gamg_repartition").c_str(),
                         ToString(GAMGrepartition[block]).c_str());
  if (!GAMGuseAggGasm.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_gamg_use_agg_gasm").c_str(),
                         ToString(GAMGuseAggGasm[block]).c_str());
  if (!GAMGreuseInterp.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_gamg_reuse_interpolation").c_str(),
                         ToString(GAMGreuseInterp[block]).c_str());
  if (!MLThreshold.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_gamg_threshold").c_str(),
                         ToString(MLThreshold[block]).c_str());
  if (!mglevels.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_mg_levels").c_str(),
                         ToString(mglevels[block]).c_str());
}


void LinSolParams::setHypreOptions(const std::string& prefix, int block) const
{
  PetscOptionsSetValue(AddPrefix(prefix,"pc_hypre_type").c_str(),
                       hypretype[block].c_str());
  PetscOptionsSetValue(AddPrefix(prefix,"pc_hypre_boomeramg_max_levels").c_str(),
                       ToString(mglevels[block]).c_str());
  PetscOptionsSetValue(AddPrefix(prefix,"pc_hypre_boomeramg_max_iter").c_str(), "1");
  PetscOptionsSetValue(AddPrefix(prefix,"pc_hypre_boomeramg_tol").c_str(), "0.0");;
  if (!HypreThreshold.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_hypre_boomeramg_strong_threshold").c_str(),
                         ToString(HypreThreshold[block]).c_str());
  if (!HypreCoarsenScheme.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_hypre_boomeramg_coarsen_type").c_str(),
                         ToString(HypreCoarsenScheme[block]).c_str());
  if (!HypreNoAggCoarse.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_hypre_boomeramg_agg_nl").c_str(),
                         ToString(HypreNoAggCoarse[block]).c_str());
  if (!HypreNoPathAggCoarse.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_hypre_boomeramg_agg_num_paths").c_str(),
                         ToString(HypreNoPathAggCoarse[block]).c_str());
  if (!HypreTruncation.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_hypre_boomeramg_truncfactor").c_str(),
                         ToString(HypreTruncation[block]).c_str());
}


void LinSolParams::setupCoarseSolver(PC& pc, const std::string& prefix, int block) const
{
  if (MLCoarseSolver[block] == "OneLevelSchwarz" ||
      MLCoarseSolver[block] == "TwoLevelSchwarz") {
    KSP cksp;
    PC  cpc;
    PCMGGetCoarseSolve(pc,&cksp);
    KSPSetType(cksp,"preonly");
    KSPGetPC(cksp,&cpc);
    PCSetType(cpc,"redistribute");
    PCSetUp(cpc);

    KSP sksp;
    PC  spc;
    PCRedistributeGetKSP(cpc,&sksp);
    KSPSetTolerances(sksp,1.0e-2,PETSC_DEFAULT,PETSC_DEFAULT,10);
    KSPGetPC(sksp,&spc);
    if (MLCoarseSolver[block] == "OneLevelSchwarz") {
      PCSetType(spc,PCGASM);
      PCSetUp(spc);
    } else {
      PCSetType(spc,PCML);
      PetscOptionsSetValue(AddPrefix(prefix,"mg_coarse_pc_ml_maxNLevels").c_str(),"2");	
      PCSetFromOptions(spc);
      PCSetUp(spc);

      KSP csksp;
      PC cspc;
      PCMGGetSmoother(spc,1,&csksp);
      KSPSetType(csksp,"richardson");
      KSPGetPC(csksp,&cspc);
      PCSetType(cspc,"asm");
      PCSetUp(cspc);
    }

    KSP* subsksp;
    PC   subspc;
    PetscInt first, nlocal;
    PCGASMGetSubKSP(spc,&nlocal,&first,&subsksp);
    for (int k = 0; k < nlocal; k++) {
      KSPGetPC(subsksp[k],&subspc);
      PCSetType(subspc,PCLU);
      KSPSetType(subsksp[k],KSPPREONLY);
    }
  }

  if (!MLCoarsePackage.empty()) {
    KSP cksp;
    PC  cpc;
    PCMGGetCoarseSolve(pc,&cksp);
    KSPGetPC(cksp,&cpc);
    PCSetType(cpc,PCLU);
    PCFactorSetMatSolverPackage(cpc,MLCoarsePackage[block].c_str());
    PCSetUp(cpc);
  }
}


void LinSolParams::setupSmoothers(PC& pc, int block, ISMat& dirIndexSet,
                                  const PetscIntMat& locSubdDofs,
                                  const PetscIntMat& subdDofs) const
{
  PetscInt n;
  PCMGGetLevels(pc,&n);

  // Presmoother settings
  for (int i = 1;i < n;i++) {
    KSP preksp;
    PC  prepc;

    // Set smoother
    std::string smoother;
    PetscInt noSmooth;

    PCMGGetSmoother(pc,i,&preksp);

    std::string compare = mgKsp[(size_t)block<mgKsp.size()-1?block:0];

    // warn that richardson might break symmetry if the KSP is CG
    if (compare == "defrichardson" && method == KSPCG)
      std::cerr << "WARNING: Using multigrid with Richardson on sublevels.\n"
                << "If you get divergence with KSP_DIVERGED_INDEFINITE_PC, try\n"
                << "adding <mgksp>chebyshev</mgksp. Add <mgksp>richardson</mgksp>\n"
                << "to quell this warning." << std::endl;

    if (compare == "richardson" || compare == "defrichardson")
      KSPSetType(preksp,KSPRICHARDSON);
    else if (compare == "chebyshev")
      KSPSetType(preksp,KSPCHEBYSHEV);

    if ((i == n-1) && (!finesmoother.empty())) {
      smoother = finesmoother[block];
      noSmooth = noFineSmooth[block];
    }
    else {
      smoother = presmoother[block];
      noSmooth = noPreSmooth[block];
    }

    KSPSetTolerances(preksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,noSmooth);
    KSPGetPC(preksp,&prepc);

    if (smoother == "asm" || smoother == "asmlu")
      setupAdditiveSchwarz(prepc, overlap[0], smoother == "asmlu",
                           locSubdDofs, subdDofs, true);
    else if (smoother == "compositedir" && (i==n-1)) {
      Mat mat;
      Mat Pmat;
#if PETSC_VERSION_MINOR < 5
      MatStructure flag;
      PCGetOperators(prepc,&mat,&Pmat,&flag);
#else
      PCGetOperators(prepc,&mat,&Pmat);
#endif

      addDirSmoother(prepc,Pmat,block,dirIndexSet);
    }
    else
      PCSetType(prepc,smoother.c_str());

    PCFactorSetLevels(prepc,levels[block]);
    KSPSetUp(preksp);
  }
}


void LinSolParams::setupAdditiveSchwarz(PC& pc, int overlap, bool asmlu,
                                        const PetscIntMat& locSubdDofs,
                                        const PetscIntMat& subdDofs, bool smoother) const
{
  PCSetType(pc, PCASM);
  if (!smoother)
    PCASMSetType(pc,PC_ASM_BASIC);
  PCASMSetOverlap(pc,overlap);

  if (!locSubdDofs.empty() && !subdDofs.empty()) {
    const size_t nsubds = subdDofs.size();

    IS isLocSubdDofs[nsubds], isSubdDofs[nsubds];
    for (size_t i = 0;i < nsubds;i++) {
      ISCreateGeneral(PETSC_COMM_SELF,locSubdDofs[i].size(),
                      &(const_cast<PetscIntMat&>(locSubdDofs)[i][0]),
                      PETSC_USE_POINTER,&(isLocSubdDofs[i]));
      ISCreateGeneral(PETSC_COMM_SELF,subdDofs[i].size(),
                      &(const_cast<PetscIntMat&>(subdDofs)[i][0]),
                      PETSC_USE_POINTER,&(isSubdDofs[i]));
    }
    PCASMSetLocalSubdomains(pc,nsubds,isSubdDofs,isLocSubdDofs);
  }

  PCSetFromOptions(pc);
  PCSetUp(pc);

  if (asmlu) {
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
#endif
