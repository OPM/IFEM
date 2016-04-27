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
#include "tinyxml.h"
#include <fstream>
#include <sstream>
#include <utility>
#include <iterator>


bool LinSolParams::read (std::istream& is, int nparam)
{
  return false;
}


bool LinSolParams::BlockParams::read(const TiXmlElement* elem)
{
  utl::getAttribute(elem, "basis", basis);
  utl::getAttribute(elem, "components", comps);

  const char* value;
  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if ((value = utl::getValue(child,"pc")))
      prec = value;
    else if ((value = utl::getValue(child,"package")))
      package = value;
    else if ((value = utl::getValue(child,"ilu_fill_level")))
      ilu_fill_level = atoi(value);
    else if ((value = utl::getValue(child,"mglevels")))
      mglevel = atoi(value);
    else if ((value = utl::getValue(child,"noPreSmooth")))
      noPreSmooth = atoi(value);
    else if ((value = utl::getValue(child,"noPostSmooth")))
      noPostSmooth = atoi(value);
    else if ((value = utl::getValue(child,"noFineSmooth")))
      noFineSmooth = atoi(value);
    else if ((value = utl::getValue(child,"presmoother")))
      presmoother = value;
    else if ((value = utl::getValue(child,"postsmoother")))
      postsmoother = value;
    else if ((value = utl::getValue(child,"finesmoother")))
      finesmoother = value;
    else if ((value = utl::getValue(child,"mgksp")))
      mgKSP = value;
    else if (!strcasecmp(child->Value(),"dirsmoother")) {
      size_t order;
      std::string type;

      if (!utl::getAttribute(child,"type",type))
        return false;
      if (!utl::getAttribute(child,"order",order))
        return false;

      dirSmoother.push_back(DirSmoother(order, type));
    }
    else if ((value = utl::getValue(child,"overlap")))
      overlap = atoi(value);
    else if ((value = utl::getValue(child,"nx")))
      subdomains[0] = atoi(value);
    else if ((value = utl::getValue(child,"ny")))
      subdomains[1] = atoi(value);
    else if ((value = utl::getValue(child,"nz")))
      subdomains[2] = atoi(value);
    else if (!strcasecmp(child->Value(),"gamg"))
      gamg.read(child);
    else if (!strcasecmp(child->Value(),"ml"))
      ml.read(child);
    else if (!strcasecmp(child->Value(),"hypre"))
      hypre.read(child);
    else if ((value = utl::getValue(child,"maxCoarseSize")))
      maxCoarseSize = atoi(value);
#ifdef HAS_PETSC
    else if ((value = utl::getValue(child,"nullspace"))) {
       if (!strcasecmp(value,"constant"))
        nullspace = CONSTANT;
      else if (!strcasecmp(value,"rigid_body"))
        nullspace = RIGID_BODY;
      else
        nullspace = NONE;
    }
#endif
    else
      return false;

  return true;
}


void LinSolParams::BlockParams::GAMGSettings::read(const TiXmlElement* elem)
{
  const char* value;
  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if ((value = utl::getValue(child,"type")))
      type = value;
    else if ((value = utl::getValue(child,"repartition")))
      repartition = atoi(value);
    else if ((value = utl::getValue(child,"use_agg_smoother")))
      useAggGasm = atoi(value);
    else if ((value = utl::getValue(child,"proc_eq_limit")))
      procEqLimit = atoi(value);
    else if ((value = utl::getValue(child,"reuse_interpolation")))
      reuseInterp = atoi(value);
    else if ((value = utl::getValue(child,"threshold")))
      threshold = atof(value);
}


void LinSolParams::BlockParams::MLSettings::read(const TiXmlElement* elem)
{
  const char* value;
  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if ((value = utl::getValue(child,"coarse_package")))
      coarsePackage = value;
    else if ((value = utl::getValue(child,"coarse_solver")))
      coarseSolver = value;
    else if ((value = utl::getValue(child,"coarsen_scheme")))
      coarsenScheme = value;
    else if ((value = utl::getValue(child,"symmetrize")))
      symmetrize = atoi(value);
    else if ((value = utl::getValue(child,"repartition")))
      repartition = atoi(value);
    else if ((value = utl::getValue(child,"block_scaling")))
      blockScaling = atoi(value);
    else if ((value = utl::getValue(child,"put_on_single_proc")))
      putOnSingleProc = atoi(value);
    else if ((value = utl::getValue(child,"reuse_interpolation")))
      reuseInterp = atoi(value);
    else if ((value = utl::getValue(child,"reusable")))
      reusable = atoi(value);
    else if ((value = utl::getValue(child,"keep_agg_info")))
      keepAggInfo = atoi(value);
    else if ((value = utl::getValue(child,"threshold")))
      threshold = atof(value);
    else if ((value = utl::getValue(child,"damping_factor")))
      dampingFactor = atof(value);
    else if ((value = utl::getValue(child,"repartition_ratio")))
      repartitionRatio = atof(value);
    else if ((value = utl::getValue(child,"aux")))
      aux = atoi(value);
    else if ((value = utl::getValue(child,"aux_threshold")))
      auxThreshold = atof(value);
}


void LinSolParams::BlockParams::HypreSettings::read(const TiXmlElement* elem)
{
  const char* value;
  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if ((value = utl::getValue(child,"type")))
      type = value;
    else if ((value = utl::getValue(child,"no_agg_coarse")))
      noAggCoarse = atoi(value);
    else if ((value = utl::getValue(child,"no_path_agg_coarse")))
      noPathAggCoarse = atof(value);
    else if ((value = utl::getValue(child,"truncation")))
      truncation = atof(value);
    else if ((value = utl::getValue(child,"threshold")))
      threshold = atof(value);
    else if ((value = utl::getValue(child,"coarsen_scheme")))
      coarsenScheme = value;
}


bool LinSolParams::read (const TiXmlElement* elem)
{
  const char* value = 0;

  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if ((value = utl::getValue(child,"type")))
      method = value;
    else if ((value = utl::getValue(child,"noStepBeforeReset")))
      nResetSolver = atoi(value);
    else if ((value = utl::getValue(child,"pc")))
      prec = value;
#ifdef HAS_PETSC
    else if ((value = utl::getValue(child,"schurpc"))) {
      if (!strncasecmp(value,"SIMPLE",6))
        schurPrec = SIMPLE;
      else if (!strncasecmp(value,"MSIMPLER",8))
        schurPrec = MSIMPLER;
      else if (!strncasecmp(value,"PCD",3))
        schurPrec = PCD;
    }
#endif
    else if (!strcasecmp(child->Value(),"block")) {
      blocks.resize(blocks.size()+1);
      blocks.back().read(child);
    }
    else if ((value = utl::getValue(child,"atol")))
      atol = atof(value);
    else if ((value = utl::getValue(child,"rtol")))
      rtol = atof(value);
    else if ((value = utl::getValue(child,"dtol")))
      dtol = atof(value);
    else if ((value = utl::getValue(child,"maxits")))
      maxIts = atoi(value);

  if (blocks.size() == 0) {
    blocks.resize(1);
    blocks.back().read(elem);
  }

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
    PCSetType(pc,this->getPreconditioner());

#if PETSC_HAVE_HYPRE
  if (!strncasecmp(prec.c_str(),"hypre",5)) {
    PCHYPRESetType(pc,blocks[0].hypre.type.c_str());
    setHypreOptions("", 0);
  }
#endif

  if (!strncasecmp(prec.c_str(),"gamg",4) || !strncasecmp(prec.c_str(),"ml",2) ) {
    PetscInt nloc = coords.size()/nsd;
    PCSetCoordinates(pc,nsd,nloc,&coords[0]);
    PCGAMGSetType(pc, PCGAMGAGG); // TODO?
  }

  if (!strncasecmp(prec.c_str(),"asm",3) || !strncasecmp(prec.c_str(),"gasm",4))
    setupAdditiveSchwarz(pc, blocks[0].overlap, !strncasecmp(prec.c_str(),"asmlu",5),
                         locSubdDofs, subdDofs, false);
  else if (!strncasecmp(prec.c_str(),"ml",2) || !strncasecmp(prec.c_str(),"gamg",4)) {
    if (!strncasecmp(prec.c_str(),"ml",2))
      setMLOptions("", 0);
    else if (!strncasecmp(prec.c_str(),"gamg",4))
      setGAMGOptions("", 0);

    PCSetFromOptions(pc);
    PCSetUp(pc);

    // Settings for coarse solver
    if (!blocks[0].ml.coarseSolver.empty())
      setupCoarseSolver(pc, "", 0);

    setupSmoothers(pc, 0, dirIndexSet, locSubdDofs, subdDofs);
  }
  else {
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
  for (size_t k = 0;k < blocks[block].dirSmoother.size();k++)
    PCCompositeAddPC(pc,"shell");
  for (size_t k = 0;k < blocks[block].dirSmoother.size();k++) {
    PC dirpc;
    PCCompositeGetPC(pc,k,&dirpc);
    PCPerm *pcperm;
    PCPermCreate(&pcperm);
    PCShellSetApply(dirpc,PCPermApply);
    PCShellSetContext(dirpc,pcperm);
    PCShellSetDestroy(dirpc,PCPermDestroy);
    PCShellSetName(dirpc,"dir");
    PCPermSetUp(dirpc,&dirIndexSet[block][k],P,blocks[block].dirSmoother[k].type.c_str());
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
#if PETSC_VERSION_MINOR > 6
#define PetscOptionsSetValue(x,y) PetscOptionsSetValue(nullptr, x, y)
#endif
  PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_maxNLevels").c_str(),
                       ToString(blocks[block].mglevel).c_str());
  PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_maxCoarseSize").c_str(),
                       ToString(blocks[block].maxCoarseSize).c_str());
  if (!blocks[block].ml.coarsenScheme.empty())
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_CoarsenScheme").c_str(),
                         blocks[block].ml.coarsenScheme.c_str());
  if (blocks[block].ml.threshold > -1.0)
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_Threshold").c_str(),
                         ToString(blocks[block].ml.threshold).c_str());
  if (blocks[block].ml.dampingFactor > -1.0)
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_DampingFactor").c_str(),
                         ToString(blocks[block].ml.dampingFactor).c_str());
  if (blocks[block].ml.repartitionRatio > -1.0)
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_repartitionMaxMinRatio").c_str(),
                         ToString(blocks[block].ml.repartitionRatio).c_str());
  if (blocks[block].ml.symmetrize > -1)
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_Symmetrize").c_str(),
                         ToString(blocks[block].ml.symmetrize).c_str());
  if (blocks[block].ml.repartition > -1)
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_repartition").c_str(),
                         ToString(blocks[block].ml.repartition).c_str());
  if (blocks[block].ml.blockScaling > -1)
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_BlockScaling").c_str(),
                         ToString(blocks[block].ml.blockScaling).c_str());
  if (blocks[block].ml.putOnSingleProc > -1)
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_repartitionPutOnSingleProc").c_str(),
                         ToString(blocks[block].ml.putOnSingleProc).c_str());
  if (blocks[block].ml.reuseInterp > -1)
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_reuse_interpolation").c_str(),
                         ToString(blocks[block].ml.reuseInterp).c_str());
  if (blocks[block].ml.keepAggInfo> -1)
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_KeepAggInfo").c_str(),
                         ToString(blocks[block].ml.keepAggInfo).c_str());
  if (blocks[block].ml.reusable > -1)
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_Reusable").c_str(),
                         ToString(blocks[block].ml.reusable).c_str());
  if (blocks[block].ml.aux > -1) {
    PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_Aux").c_str(),
                         ToString(blocks[block].ml.aux).c_str());
    if (blocks[block].ml.auxThreshold > -1)
      PetscOptionsSetValue(AddPrefix(prefix,"pc_ml_AuxThreshold").c_str(),
                           ToString(blocks[block].ml.auxThreshold).c_str());
  }
}


template<class T>
static void condSetup(const std::string& prefix, const std::string& option, const T& var)
{
  if (var > -1)
    PetscOptionsSetValue(AddPrefix(prefix,option).c_str(), ToString(var).c_str());
};


template<>
void condSetup(const std::string& prefix, const std::string& option, const std::string& var)
{
  if (!var.empty())
    PetscOptionsSetValue(AddPrefix(prefix,option).c_str(), var.c_str());
};


void LinSolParams::setGAMGOptions(const std::string& prefix, size_t iblock) const
{
  if (iblock >= blocks.size())
    return;
  const BlockParams& block = blocks[iblock];

  condSetup(prefix, "pc_gamg_type", block.gamg.type);
  condSetup(prefix, "pc_gamg_coarse_eq_limit", block.maxCoarseSize);
  condSetup(prefix, "pc_gamg_process_eq_limit", block.gamg.procEqLimit);
  condSetup(prefix, "pc_gamg_repartition", block.gamg.repartition);
  condSetup(prefix, "pc_gamg_use_agg_gasm", block.gamg.useAggGasm);
  condSetup(prefix, "pc_gamg_reuse_interpolation", block.gamg.reuseInterp);
  condSetup(prefix, "pc_gamg_threshold", block.gamg.threshold);
  condSetup(prefix, "pc_mg_levels", block.mglevel);
}


void LinSolParams::setHypreOptions(const std::string& prefix, size_t iblock) const
{
  if (iblock >= blocks.size())
    return;
  const BlockParams& block = blocks[iblock];

  condSetup(prefix,"pc_hypre_type", block.hypre.type);
  condSetup(prefix, "pc_hypre_boomeramg_max_levels", block.mglevel);
  // TODO: Investigate why these are hardcoded.
  PetscOptionsSetValue(AddPrefix(prefix,"pc_hypre_boomeramg_max_iter").c_str(), "1");
  PetscOptionsSetValue(AddPrefix(prefix,"pc_hypre_boomeramg_tol").c_str(), "0.0");;
  condSetup(prefix, "pc_hypre_boomeramg_strong_threshold", block.hypre.threshold);
  condSetup(prefix, "pc_hypre_boomeramg_coarsen_type", block.hypre.coarsenScheme);
  condSetup(prefix, "pc_hypre_boomeramg_agg_nl", block.hypre.noAggCoarse);
  condSetup(prefix, "pc_hypre_boomeramg_agg_num_paths", block.hypre.noPathAggCoarse);
  condSetup(prefix, "pc_hypre_boomeramg_truncfactor", block.hypre.truncation);
}


void LinSolParams::setupCoarseSolver(PC& pc, const std::string& prefix, size_t iblock) const
{
  if (iblock >= blocks.size())
    return;
  const BlockParams& block = blocks[iblock];

  if (block.ml.coarseSolver == "OneLevelSchwarz" ||
      block.ml.coarseSolver == "TwoLevelSchwarz") {
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
    if (block.ml.coarseSolver == "OneLevelSchwarz") {
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

  if (!block.ml.coarsePackage.empty()) {
    KSP cksp;
    PC  cpc;
    PCMGGetCoarseSolve(pc,&cksp);
    KSPGetPC(cksp,&cpc);
    PCSetType(cpc,PCLU);
    PCFactorSetMatSolverPackage(cpc,block.ml.coarsePackage.c_str());
    PCSetUp(cpc);
  }
}


void LinSolParams::setupSmoothers(PC& pc, size_t iblock, ISMat& dirIndexSet,
                                  const PetscIntMat& locSubdDofs,
                                  const PetscIntMat& subdDofs) const
{
  if (iblock >= blocks.size())
    return;
  const BlockParams& block = blocks[iblock];

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

    // warn that richardson might break symmetry if the KSP is CG
    if (block.mgKSP == "defrichardson" && method == KSPCG)
      std::cerr << "WARNING: Using multigrid with Richardson on sublevels.\n"
                << "If you get divergence with KSP_DIVERGED_INDEFINITE_PC, try\n"
                << "adding <mgksp>chebyshev</mgksp. Add <mgksp>richardson</mgksp>\n"
                << "to quell this warning." << std::endl;

    if (block.mgKSP == "richardson" || block.mgKSP == "defrichardson")
      KSPSetType(preksp,KSPRICHARDSON);
    else if (block.mgKSP == "chebyshev")
      KSPSetType(preksp,KSPCHEBYSHEV);

    if ((i == n-1) && (!block.finesmoother.empty())) {
      smoother = block.finesmoother;
      noSmooth = block.noFineSmooth;
    }
    else {
      smoother = block.presmoother;
      noSmooth = block.noPreSmooth;
    }

    KSPSetTolerances(preksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,noSmooth);
    KSPGetPC(preksp,&prepc);

    if (smoother == "asm" || smoother == "asmlu")
      setupAdditiveSchwarz(prepc, block.overlap, smoother == "asmlu",
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

      addDirSmoother(prepc,Pmat,iblock,dirIndexSet);
    }
    else
      PCSetType(prepc,smoother.c_str());

    PCFactorSetLevels(prepc,block.ilu_fill_level);
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
