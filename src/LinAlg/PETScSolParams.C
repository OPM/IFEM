// $Id$
//==============================================================================
//!
//! \file PETScSolParams.C
//!
//! \date Mar 10 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Linear solver parameters for PETSc.
//!
//==============================================================================

#include "PETScSolParams.h"
#include "PETScPCPerm.h"
#include "ASMstruct.h"
#include "LinSolParams.h"
#include "ProcessAdm.h"
#include "SAMpatch.h"
#include "Utilities.h"
#include "tinyxml.h"
#include <fstream>
#include <sstream>
#include <utility>
#include <iterator>


void PETScSolParams::setupPC(PC& pc,
                             size_t block,
                             const std::string& prefix,
                             const std::set<int>& blockEqs)
{
  // Set preconditioner
  std::string prec = params.getBlock(block).getStringValue("pc");
  if (prec == "default")
    prec = adm.getNoProcs() > 1 ? "asm" : "ilu";

  if (prec == "compositedir") {
    Mat mat;
    Mat Pmat;
#if PETSC_VERSION_MINOR < 5
    MatStructure flag;
    PCGetOperators(pc,&mat,&Pmat,&flag);
#else
    PCGetOperators(pc,&mat,&Pmat);
#endif
    //TODO: Dir smoothers are complicated by condensated DOFs
    //addDirSmoother(pc,Pmat,params.getBlock(block),block,dirIndexSet);
  }
  else if (prec == "asmlu")
    PCSetType(pc,"asm");
  else
    PCSetType(pc,prec.c_str());

#if PETSC_HAVE_HYPRE
  if (prec == "hypre") {
    PCHYPRESetType(pc,params.getBlock(block).getStringValue("hypre_type").c_str());
    setHypreOptions(prefix, params.getBlock(block));
  }
#endif

  // TODO: coordinate based agglomeration won't work with condensated dofs.
//  if (prec == "gamg" || prec == "ml") {
//    PetscInt nloc = coords.size()/nsd;
//    PCSetCoordinates(pc,nsd,nloc,&coords[0]);
//    PCGAMGSetType(pc, PCGAMGAGG);
//  }

  if (prec == "asm" || prec == "gasm" || prec == "asmlu")
    setupAdditiveSchwarz(pc, block, prec == "asmlu", false, blockEqs);
  else if (prec == "ml")
    setMLOptions(prefix, params.getBlock(block));
  else if (prec == "gamg")
    setGAMGOptions(prefix, params.getBlock(block));
  else if (prec == "ilu")
    PCFactorSetLevels(pc,params.getBlock(block).getIntValue("ilu_fill_level"));

  PCSetFromOptions(pc);
  PCSetUp(pc);

  // Settings for coarse solver
  if ((prec == "ml" || prec == "gamg")) {
    if (params.getBlock(block).hasValue("multigrid_coarse_solver"))
      setupCoarseSolver(pc, prefix, params.getBlock(block));
    // TODO: dir smoothers
    setupSmoothers(pc, block, ISMat(), blockEqs);
  }
}


bool PETScSolParams::addDirSmoother(PC pc, const Mat& P,
                                    int iBlock, const ISMat& dirIndexSet)
{
  const LinSolParams::BlockParams& block = params.getBlock(iBlock);
  PCSetType(pc,"composite");
  PCCompositeSetType(pc,PC_COMPOSITE_MULTIPLICATIVE);
  for (size_t k = 0;k < block.dirSmoother.size();k++)
    PCCompositeAddPC(pc,"shell");
  for (size_t k = 0;k < block.dirSmoother.size();k++) {
    PC dirpc;
    PCCompositeGetPC(pc,k,&dirpc);
    PCPerm *pcperm;
    PCPermCreate(&pcperm);
    PCShellSetApply(dirpc,PCPermApply);
    PCShellSetContext(dirpc,pcperm);
    PCShellSetDestroy(dirpc,PCPermDestroy);
    PCShellSetName(dirpc,"dir");
    PCPermSetUp(dirpc,const_cast<IS*>(&dirIndexSet[iBlock][k]),P,block.dirSmoother[k].type.c_str());
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


//! \brief Conditionally add a setting from map to PETsc.
//! \param[in] prefix Prefix for petsc setting
//! \param[in] petsc_option Name of option in PETsc
//! \param[in] map_option Name of option in setting map
//! \param[in] map The setting map
static void condSetup(const std::string& prefix, const std::string& petsc_option,
                      const std::string& map_option, const SettingMap& map)
{
#if PETSC_VERSION_MINOR > 6
#define PetscOptionsSetValue(x,y) PetscOptionsSetValue(nullptr, x, y)
#endif
  if (map.hasValue(map_option))
    PetscOptionsSetValue(AddPrefix(prefix,petsc_option).c_str(), map.getStringValue(map_option).c_str());
}


void PETScSolParams::setMLOptions(const std::string& prefix, const SettingMap& map)
{
  condSetup(prefix, "pc_ml_maxNLevels", "multigrid_levels", map);
  condSetup(prefix, "pc_ml_maxCoarseSize", "multigrid_max_coarse_size", map);
  condSetup(prefix, "pc_ml_maxCoarsenScheme", "ml_coarsen_scheme", map);
  condSetup(prefix, "pc_ml_Threshold", "ml_threshold", map);
  condSetup(prefix, "pc_ml_DampingFactor", "ml_damping_factor", map);
  condSetup(prefix, "pc_ml_repartitionMaxMinRatio", "ml_repartition_max_min_ratio", map);
  condSetup(prefix, "pc_ml_Symmetrize", "ml_symmetrize", map);
  condSetup(prefix, "pc_ml_repartition", "ml_repartition", map);
  condSetup(prefix, "pc_ml_BlockScaling", "ml_block_scaling", map);
  condSetup(prefix, "pc_ml_repartitionPutOnSingleProc", "ml_put_on_single_proc", map);
  condSetup(prefix, "pc_ml_reuse_interpolation", "ml_reuse_interpolation", map);
  condSetup(prefix, "pc_ml_KeepAggInfo", "ml_keep_agg_info", map);
  condSetup(prefix, "pc_ml_Reusable", "ml_reusable", map);
  condSetup(prefix, "pc_ml_Aux", "ml_aux", map);
  condSetup(prefix, "pc_ml_AuxThreshold", "ml_aux_threshold", map);
}


void PETScSolParams::setGAMGOptions(const std::string& prefix, const SettingMap& map)
{
  condSetup(prefix, "pc_gamg_type", "gamg_type", map);
  condSetup(prefix, "pc_gamg_coarse_eq_limit", "gamg_coarse_eq_limit", map);
  condSetup(prefix, "pc_gamg_process_eq_limit", "gamg_process_eq_limit", map);
  condSetup(prefix, "pc_gamg_repartition", "gamg_repartition", map);
  condSetup(prefix, "pc_gamg_use_agg_gasm", "gamg_use_agg_gasm", map);
  condSetup(prefix, "pc_gamg_reuse_interpolation", "gamg_reuse_interpolation", map);
  condSetup(prefix, "pc_gamg_threshold", "gamg_threshold", map);
  condSetup(prefix, "pc_mg_levels", "multigrid_levels", map);
}


void PETScSolParams::setHypreOptions(const std::string& prefix, const SettingMap& map)
{
  condSetup(prefix,"pc_hypre_type", "hypre_type", map);
  condSetup(prefix, "pc_hypre_boomeramg_max_levels", "multigrid_levels", map);
  // TODO: Investigate why these are hardcoded.
  PetscOptionsSetValue(AddPrefix(prefix,"pc_hypre_boomeramg_max_iter").c_str(), "1");
  PetscOptionsSetValue(AddPrefix(prefix,"pc_hypre_boomeramg_tol").c_str(), "0.0");;
  condSetup(prefix, "pc_hypre_boomeramg_strong_threshold", "hypre_threshold", map);
  condSetup(prefix, "pc_hypre_boomeramg_coarsen_type", "hypre_boomeramg_coarsen_type", map);
  condSetup(prefix, "pc_hypre_boomeramg_agg_nl", "hypre_boomeramg_agg_nl", map);
  condSetup(prefix, "pc_hypre_boomeramg_agg_num_paths", "hypre_boomeramg_agg_num_paths", map);
  condSetup(prefix, "pc_hypre_boomeramg_truncfactor", "hypre_boomeramg_truncation_factor", map);
}


void PETScSolParams::setupCoarseSolver(PC& pc, const std::string& prefix, const SettingMap& map)
{
  std::string coarseSolver = map.getStringValue("multigrid_coarse_solver");
  std::string coarsePackage = map.getStringValue("multigrid_coarse_package");
  if (coarseSolver == "OneLevelSchwarz" ||
      coarseSolver == "TwoLevelSchwarz") {
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
    if (coarseSolver == "OneLevelSchwarz") {
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
  } else if (coarseSolver == "lu") {
    KSP cksp;
    PC  cpc;
    PCMGGetCoarseSolve(pc,&cksp);
    KSPSetType(cksp,"preonly");
    KSPGetPC(cksp,&cpc);
    PCSetType(cpc,"lu");
    PCSetUp(cpc);
  }

  if (!coarsePackage.empty()) {
    KSP cksp;
    PC  cpc;
    PCMGGetCoarseSolve(pc,&cksp);
    KSPGetPC(cksp,&cpc);
    PCSetType(cpc,PCLU);
#if PETSC_VERSION_MINOR >= 9
    PCFactorSetMatSolverType(cpc,coarsePackage.c_str());
#else
    PCFactorSetMatSolverPackage(cpc,coarsePackage.c_str());
#endif
    PCSetUp(cpc);
  }
}


void PETScSolParams::setupSmoothers(PC& pc, size_t iBlock,
                                    const ISMat& dirIndexSet,
                                    const std::set<int>& blockEqs)
{
  PetscInt n;
  PCMGGetLevels(pc,&n);

  std::string mgKSP = params.getBlock(iBlock).getStringValue("multigrid_ksp");

  // warn that richardson might break symmetry if the KSP is CG
  if (mgKSP == "defrichardson" && params.getStringValue("type") == KSPCG)
    std::cerr << "**PETScMatrix** WARNING: Using multigrid with Richardson on sublevels.\n"
              << "If you get divergence with KSP_DIVERGED_INDEFINITE_PC, try\n"
              << "adding ksp=\"chebyshev\" to the multigrid tag.\n"
              << "Add ksp=\"richardson\" to quell this warning." << std::endl;

  // Presmoother settings
  for (int i = 1;i < n;i++) {
    KSP preksp;
    PC  prepc;

    // Set smoother
    std::string smoother;
    PetscInt noSmooth;

    PCMGGetSmoother(pc,i,&preksp);

    if (mgKSP == "defrichardson")
      KSPSetType(preksp,KSPRICHARDSON);
    else
      KSPSetType(preksp,mgKSP.c_str());

    std::string finesmoother = params.getBlock(iBlock).getStringValue("multigrid_finesmoother");
    if ((i == n-1) && !finesmoother.empty()) {
      smoother = finesmoother;
      noSmooth = params.getBlock(iBlock).getIntValue("multigrid_no_fine_smooth");
    }
    else {
      smoother = params.getBlock(iBlock).getStringValue("multigrid_smoother");
      noSmooth = params.getBlock(iBlock).getIntValue("multigrid_no_smooth");;
    }

    if (smoother.empty()) {
      if (adm.getNoProcs() > 1)
        smoother = "asm";
      else
        smoother = "ilu";
    }

    KSPSetTolerances(preksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,noSmooth);
    KSPGetPC(preksp,&prepc);

    if (smoother == "asm" || smoother == "asmlu")
      setupAdditiveSchwarz(prepc, iBlock, smoother == "asmlu", true, blockEqs);
    else if (smoother == "compositedir" && (i==n-1)) {
      Mat mat;
      Mat Pmat;
#if PETSC_VERSION_MINOR < 5
      MatStructure flag;
      PCGetOperators(prepc,&mat,&Pmat,&flag);
#else
      PCGetOperators(prepc,&mat,&Pmat);
#endif

      //TODO: dir smoothers are complicated by condensated DOFs
      //addDirSmoother(prepc,Pmat,params.getBlock(iBlock),iBlock,dirIndexSet);
    }
    else
      PCSetType(prepc,smoother.c_str());

    PCFactorSetLevels(prepc,params.getBlock(0).getIntValue("ilu_fill_level"));
    KSPSetUp(preksp);
  }
}


void PETScSolParams::setupAdditiveSchwarz(PC& pc, size_t block,
                                          bool asmlu,
                                          bool smoother,
                                          const std::set<int>& blockEqs)
{
  int overlap = params.getBlock(block).getIntValue("asm_overlap");
  int nx = params.getBlock(block).getIntValue("asm_nx");
  int ny = params.getBlock(block).getIntValue("asm_ny");
  int nz = params.getBlock(block).getIntValue("asm_nz");
  std::vector<std::set<int>> subdDofs;
  if (nx+ny+nz > 0)
    subdDofs = adm.dd.getSubdomains(nx, ny, nz, overlap, block);

  PCSetType(pc, PCASM);
  if (!smoother)
    PCASMSetType(pc,PC_ASM_BASIC);
  PCASMSetOverlap(pc,overlap);

  if (!subdDofs.empty()) {
    std::vector<IS> isLocSubdDofs(subdDofs.size());
    std::vector<IS> isSubdDofs(subdDofs.size());
    size_t i = 0;
    for (auto& it : subdDofs) {
      IntVec subdofs, locSubdDofs;
      locSubdDofs.reserve(it.size());
      for (auto& it2 : it) {
        int geq = adm.dd.getGlobalEq(it2,block);
        if (it2 >= adm.dd.getMinEq(block))
          locSubdDofs.push_back(geq-1);
        subdofs.push_back(geq-1);
      }

      ISCreateGeneral(PETSC_COMM_SELF,locSubdDofs.size(),
                      &locSubdDofs[0],
                      PETSC_USE_POINTER,&isLocSubdDofs[i]);
      ISCreateGeneral(PETSC_COMM_SELF,subdofs.size(),
                      &subdofs[0],
                      PETSC_USE_POINTER,&isSubdDofs[i++]);
    }
    PCASMSetLocalSubdomains(pc,subdDofs.size(),isSubdDofs.data(),isLocSubdDofs.data());
  }

  PCSetFromOptions(pc);
  PCSetUp(pc);

  KSP* subksp;
  PC   subpc;
  PetscInt first, nlocal;
  PCASMGetSubKSP(pc,&nlocal,&first,&subksp);

  int fill_level = params.getBlock(block).getIntValue("ilu_fill_level");
  for (int i = 0; i < nlocal; i++) {
    KSPGetPC(subksp[i],&subpc);
    if (asmlu) {
      PCSetType(subpc,PCLU);
      KSPSetType(subksp[i],KSPPREONLY);
    } else
      PCFactorSetLevels(subpc,fill_level);
  }
}
