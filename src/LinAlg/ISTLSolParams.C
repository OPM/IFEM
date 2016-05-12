// $Id$
//==============================================================================
//!
//! \file ISTLSolParams.C
//!
//! \date Mar 10 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Linear solver parameters for ISTL.
//!
//==============================================================================

#include "ISTLSolParams.h"
#include "ASMstruct.h"
#include "LinSolParams.h"
#include "ProcessAdm.h"
#include "SAMpatch.h"
#include <dune/istl/overlappingschwarz.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/matrixmatrix.hh>
#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 3)
#include <dune/istl/paamg/fastamg.hh>
#include <dune/istl/paamg/twolevelmethod.hh>
#endif


namespace ISTL {

BlockPreconditioner::BlockPreconditioner(const ISTL::Mat& A,
                                         const DomainDecomposition& dd_,
                                         const std::string& schurType) :
  dd(dd_)
{
  blocks.resize(dd.getNoBlocks()*dd.getNoBlocks()+1);
  size_t k = 0;
  for (size_t i = 0; i < dd.getNoBlocks(); ++i)
    for (size_t j = 0; j < dd.getNoBlocks(); ++j)
      extractBlock(blocks[k++], A, dd.getBlockEqs(i), dd.getBlockEqs(j), dd.getG2LEQ(j));

  // Scale blocks[1] with inverse diagonal of blocks[0]
  blocks.back() = blocks[1];
  for (auto row = blocks.back().begin(); row != blocks.back().end(); ++row)
    for (auto it = row->begin(); it != row->end(); ++it)
      *it /= blocks[0][row.index()][row.index()];

  // blocks[1] = D - C*diag(A)^-1*B
  ISTL::Mat SP;
  Dune::matMultMat(SP, blocks[2], blocks.back());
  blocks[2] = blocks[1];
  subtractMatrices(blocks[1], blocks[3], SP);

  // clear memory for unused blocks
  blocks.resize(schurType == "diag" ? 2 : 3);

  blockPre.resize(dd.getNoBlocks());
  blockOp.resize(dd.getNoBlocks());
}


void BlockPreconditioner::pre(ISTL::Vec& x, ISTL::Vec& b)
{
  ISTL::Vec tempx, tempb;
  for (size_t block = 0; block < dd.getNoBlocks(); ++block) {
    ISTL::Vec tempx, tempb;
    tempx.resize(dd.getBlockEqs(block).size());
    tempb.resize(dd.getBlockEqs(block).size());
    size_t i = 0;
    for (auto& it : dd.getBlockEqs(block)) {
      tempx[i] = x[it-1];
      tempb[i++] = b[it-1];
    }

    blockPre[block]->pre(tempx, tempb);

    i = 0;
    for (auto& it : dd.getBlockEqs(block)) {
      x[it-1] = tempx[i];
      b[it-1] = tempb[i++];
    }
  }
}


void BlockPreconditioner::apply(ISTL::Vec& v, const ISTL::Vec& d)
{
  // backwards due to upper schur complement
  ISTL::Vec tempx, tempb;
  for (int block = dd.getNoBlocks()-1; block >= 0; --block) {
    tempb.resize(dd.getBlockEqs(block).size());
    size_t i = 0;
    for (auto& it : dd.getBlockEqs(block))
      tempb[i++] = d[it-1];

    if (block == 0 && blocks.size() > 2)
      blocks.back().mmv(tempx, tempb);

    tempx.resize(dd.getBlockEqs(block).size());
    tempx = 0;
    blockPre[block]->apply(tempx, tempb);

    i = 0;
    for (auto& it : dd.getBlockEqs(block))
      v[it-1] = tempx[i++];
  }
}


void BlockPreconditioner::post(ISTL::Vec& x)
{
  ISTL::Vec tempx, tempb;
  for (size_t block = 0; block < dd.getNoBlocks(); ++block) {
    ISTL::Vec tempx;
    tempx.resize(dd.getBlockEqs(block).size());
    size_t i = 0;
    for (auto& it : dd.getBlockEqs(block))
      tempx[i++] = x[it-1];

    blockPre[block]->post(tempx);

    i = 0;
    for (auto& it : dd.getBlockEqs(block))
      x[it-1] = tempx[i++];
  }
}


void BlockPreconditioner::extractBlock(ISTL::Mat& B, const ISTL::Mat& A,
                                       const std::set<int>& eqs_row1,
                                       const std::set<int>& eqs_col1,
                                       const std::map<int,int>& eqs_col_g2l)
{
  std::vector<int> eqs_row(eqs_row1.begin(), eqs_row1.end());
  std::vector<int> eqs_col(eqs_col1.begin(), eqs_col1.end());
  size_t sum=0;
  std::vector<std::set<int>> adj;
  adj.resize(eqs_row.size());
  size_t i = 0;

  for (auto& it3 : eqs_row) {
    auto it = A.begin()+it3-1;
    for (auto it2 = it->begin(); it2 != it->end(); ++it2) {
      auto itc = eqs_col_g2l.find(it2.index()+1);
      if (itc != eqs_col_g2l.end()) {
        adj[i].insert(itc->second-1);
        ++sum;
      }
    }
    ++i;
  }

  B.setSize(eqs_row.size(), eqs_col.size(), sum);
  B.setBuildMode(ISTL::Mat::random);

  for (size_t i = 0; i < adj.size(); i++)
    B.setrowsize(i,adj[i].size());
  B.endrowsizes();

  for (size_t i = 0; i < adj.size(); i++)
    for (auto& it : adj[i])
      B.addindex(i, it);

  B.endindices();
  B = 0;
  for (auto it = B.begin(); it != B.end(); ++it)
    for (auto it2 = it->begin(); it2 != it->end(); ++it2)
      *it2 = A[eqs_row[it.index()]-1][eqs_col[it2.index()]-1];
}


void BlockPreconditioner::subtractMatrices(ISTL::Mat& A, const ISTL::Mat& B,
                                           const ISTL::Mat& C)
{
  size_t sum=0;
  std::vector<std::set<int>> adj;
  adj.resize(B.N());
  for (auto row = B.begin(); row != B.end(); ++row) {
    for (auto it = row->begin(); it != row->end(); ++it)
      adj[row.index()].insert(it.index());
  }
  for (auto row = C.begin(); row != C.end(); ++row) {
    for (auto it = row->begin(); it != row->end(); ++it)
      adj[row.index()].insert(it.index());
    sum += adj[row.index()].size();
  }

  A.setSize(B.N(), B.M(), sum);
  A.setBuildMode(ISTL::Mat::random);

  for (size_t i = 0; i < adj.size(); i++)
    A.setrowsize(i,adj[i].size());
  A.endrowsizes();

  for (size_t i = 0; i < adj.size(); i++)
    for (auto& it : adj[i])
      A.addindex(i, it);

  A.endindices();
  A = 0;
  for (auto row = B.begin(); row != B.end(); ++row)
    for (auto it = row->begin(); it != row->end(); ++it)
      A[row.index()][it.index()] = *it;

  for (auto row = C.begin(); row != C.end(); ++row)
    for (auto it = row->begin(); it != row->end(); ++it)
      A[row.index()][it.index()] -= *it;
}


} // namespace ISTL


/*! \brief Helper template for setting up a solver with the appropriate
           preconditioner type.
    \details We cannot instance using dynamic polymorphism, the solver need
             access to the real type for the preconditioner. We can however
             call the solver in the interface class scope afterwards.
 */
template<class Prec>
static ISTL::InverseOperator* setupWithPreType(const LinSolParams& solParams,
                                               ISTL::Operator& op,
                                               Prec& pre)
{
  std::string type = solParams.getStringValue("type");
  double rtol = solParams.getDoubleValue("rtol");
  int maxits = solParams.getIntValue("maxits");
  int verbosity = solParams.getIntValue("verbosity");
  int restart = solParams.getIntValue("gmres_restart_iterations");
  if (type == "bcgs")
    return new Dune::BiCGSTABSolver<ISTL::Vec>(op, pre, rtol, maxits, verbosity);
  else if (type == "cg")
    return new Dune::CGSolver<ISTL::Vec>(op, pre, rtol, maxits, verbosity);
  else if (type == "minres")
    return new Dune::MINRESSolver<ISTL::Vec>(op, pre, rtol, maxits, verbosity);
  else if (type == "gmres")
    return new Dune::RestartedGMResSolver<ISTL::Vec>(op, pre, rtol, restart,
                                                     maxits, verbosity);

  return nullptr;
}


/*! \brief Helper template for setting up an AMG preconditioner with a given smoother */

template<class Smoother>
static ISTL::Preconditioner* setupAMG(const LinSolParams& params,
                                      size_t block, ISTL::Operator& op,
                                      std::unique_ptr<ISTL::InverseOperator>* solver)
{
  // The coupling metric used in the AMG
  typedef Dune::Amg::FirstDiagonal CouplingMetric;
  // The coupling criterion used in the AMG
  typedef Dune::Amg::SymmetricCriterion<ISTL::Mat, CouplingMetric> CritBase;
  // The coarsening criterion used in the AMG
  typedef Dune::Amg::CoarsenCriterion<CritBase> Criterion;
  Criterion crit;
  typedef typename Dune::Amg::AMG<ISTL::Operator, ISTL::Vec, Smoother> AMG;
  typename AMG::SmootherArgs args;
  args.relaxationFactor = 1.0;
  args.iterations = std::max(1, params.getBlock(block).getIntValue("multigrid_no_smooth"));

  if (params.getBlock(block).hasValue("multigrid_max_coarse_size"))
    crit.setCoarsenTarget(params.getBlock(block).getIntValue("multigrid_max_coarse_size"));

  if (params.getBlock(block).hasValue("multigrid_no_smooth")) {
    int val = params.getBlock(block).getIntValue("multigrid_no_smooth");
    crit.setNoPreSmoothSteps(val);
    crit.setNoPostSmoothSteps(val);
  }

  if (params.getBlock(block).hasValue("multigrid_max_coarse_size"))
    crit.setCoarsenTarget(params.getBlock(block).getIntValue("multigrid_max_coarse_size"));

  auto result = new AMG(op, crit, args);
  if (solver)
    solver->reset(setupWithPreType(params, op, *result));

  return result;
}


#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 3)
/*! \brief Helper template for setting up a AMG preconditioner
 *!        with fine smoother differing from the other levels */

template<class Smoother, class FineSmoother>
static ISTL::Preconditioner* setupAMG2_full(const LinSolParams& params, size_t block,
                                            ISTL::Operator& op,
                                            std::unique_ptr<ISTL::InverseOperator>* solver,
                                            FineSmoother* fsmooth)
{
  typedef Dune::Amg::FirstDiagonal CouplingMetric;
  typedef Dune::Amg::SymmetricCriterion<ISTL::Mat,CouplingMetric> CritBase;
  typedef Dune::Amg::CoarsenCriterion<CritBase> Criterion;
  typedef Dune::Amg::AggregationLevelTransferPolicy<ISTL::Operator,Criterion> TransferPolicy;
  typedef Dune::Amg::OneStepAMGCoarseSolverPolicy<ISTL::Operator,Smoother,Criterion> CoarsePolicy;
  typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
  typedef typename Dune::Amg::TwoLevelMethod<ISTL::Operator, CoarsePolicy, FineSmoother> AMG2;

  SmootherArgs args;
  args.relaxationFactor = 1.0;

  Criterion crit;
  if (params.getBlock(block).hasValue("multigrid_max_coarse_size"))
    crit.setCoarsenTarget(params.getBlock(block).getIntValue("multigrid_max_coarse_size"));

  if (params.getBlock(block).hasValue("multigrid_no_smooth")) {
    int val = params.getBlock(block).getIntValue("multigrid_no_smooth");
    crit.setNoPreSmoothSteps(val);
    crit.setNoPostSmoothSteps(val);
  }

  if (params.getBlock(block).hasValue("multigrid_max_coarse_size"))
    crit.setCoarsenTarget(params.getBlock(block).getIntValue("multigrid_max_coarse_size"));

  CoarsePolicy coarsePolicy(args, crit);
  TransferPolicy policy(crit);
  Dune::shared_ptr<FineSmoother> fsp(fsmooth);
  auto result = new AMG2(op, fsp, policy, coarsePolicy);
  if (solver)
    solver->reset(setupWithPreType(params, op, *result));

  return result;
}


template<class FineSmoother>
static ISTL::Preconditioner* setupAMG2_smoother(const LinSolParams& params, size_t block,
                                                ISTL::Operator& op,
                                                std::unique_ptr<ISTL::InverseOperator>* solver,
                                                FineSmoother* fsmooth)
{
  std::string smoother = params.getBlock(block).getStringValue("multigrid_smoother");
  if (smoother == "ilu")
    return setupAMG2_full<ISTL::ILU0>(params, block, op, solver, fsmooth);
  else if (smoother == "sor")
    return setupAMG2_full<ISTL::SOR>(params, block, op, solver, fsmooth);
  else if (smoother == "ssor")
    return setupAMG2_full<ISTL::SSOR>(params, block, op, solver, fsmooth);
  else if (smoother == "jacobi")
    return setupAMG2_full<ISTL::GJ>(params, block, op, solver, fsmooth);
  else {
    std::cerr << "** ISTLSolParams ** Invalid smoother " << smoother << "." << std::endl;
    return nullptr;
  }
}


static ISTL::Preconditioner* setupAMG2(const LinSolParams& params, size_t block,
                                       ISTL::Operator& op, std::unique_ptr<ISTL::InverseOperator>* solver)
{
  std::string smoother = params.getBlock(block).getStringValue("multigrid_finesmoother");
  if (params.getBlock(block).getStringValue("multigrid_finesmoother") == "ilu") {
    auto fsmooth = new ISTL::ILU0(op.getmat(), 1.0);
    return setupAMG2_smoother(params, block, op, solver, fsmooth);
  } else  if (params.getBlock(block).getStringValue("multigrid_finesmoother") == "sor") {
    auto fsmooth = new ISTL::SOR(op.getmat(), 1, 1.0);
    return setupAMG2_smoother(params, block, op, solver, fsmooth);
  } else  if (params.getBlock(block).getStringValue("multigrid_finesmoother") == "ssor") {
    auto fsmooth = new ISTL::SSOR(op.getmat(), 1, 1.0);
    return setupAMG2_smoother(params, block, op, solver, fsmooth);
  } else  if (params.getBlock(block).getStringValue("multigrid_finesmoother") == "jacobi") {
    auto fsmooth = new ISTL::GJ(op.getmat(), 1, 1.0);
    return setupAMG2_smoother(params, block, op, solver, fsmooth);
  } else {
    std::cerr << "** ISTLSolParams ** Invalid fine smoother " << smoother << "." << std::endl;
    return nullptr;
  }
}
#endif


/*! \brief Conditionally setup a KSP */
template<class Type>
static ISTL::Preconditioner*
setupSolver(Type* pre,
            ISTL::Operator& op,
            const LinSolParams& solParams,
            std::unique_ptr<ISTL::InverseOperator>* solver)
{
  if (solver)
    solver->reset(setupWithPreType(solParams, op, *pre));

  return pre;
}

ISTL::Preconditioner* ISTLSolParams::setupPCInternal(ISTL::Mat& A,
                                                     ISTL::Operator& op,
                                                     size_t block,
                                                     std::unique_ptr<ISTL::InverseOperator>* solver)
{
  std::string prec = solParams.getBlock(block).getStringValue("pc");
  if (prec == "ilu") {
    int fill_level = solParams.getBlock(block).getIntValue("ilu_fill_level");
    if (fill_level == 0)
      return setupSolver(new ISTL::ILU0(A, 1.0), op, solParams, solver);
    else
      return setupSolver(new ISTL::ILUn(A, fill_level, 1.0),
                         op, solParams, solver);
  } else if (prec == "lu")
    return setupSolver(new ISTL::LU(new ISTL::LUType(A)),
                       op, solParams, solver);
  else if (prec == "sor")
    return setupSolver(new ISTL::SOR(A, 1, 1.0), op, solParams, solver);
  else if (prec == "ssor")
    return setupSolver(new ISTL::SSOR(A, 1, 1.0), op, solParams, solver);
  else if (prec == "jacobi")
    return setupSolver(new ISTL::GJ(A, 1, 1.0), op, solParams, solver);
 else if (prec == "gs")
    return setupSolver(new ISTL::GS(A, 1, 1.0), op, solParams, solver);
  else if (prec == "asm" || prec == "asmlu") {
    size_t nx = solParams.getBlock(block).getIntValue("asm_nx");
    nx = std::max(1ul, nx);
    size_t ny = solParams.getBlock(block).getIntValue("asm_ny");
    ny = std::max(1ul, ny);
    size_t nz = solParams.getBlock(block).getIntValue("asm_nz");
    nz = std::max(1ul, nz);
    int overlap = solParams.getBlock(block).getIntValue("asm_overlap");
    auto locSubdDofs = adm.dd.getSubdomains(nx, ny, nz, overlap, block);

    if (prec == "asmlu") {
      ISTL::ASMLU::subdomain_vector ddofs(locSubdDofs.size());
      for (size_t i = 0; i < locSubdDofs.size(); ++i)
        for (const auto& it : locSubdDofs[i])
          ddofs[i].insert(it-1);

      return setupSolver(new ISTL::ASMLU(A, ddofs), op, solParams, solver);
    } else {
      ISTL::ASM::subdomain_vector ddofs(locSubdDofs.size());
      for (size_t i = 0; i < locSubdDofs.size(); ++i) {
        for (size_t i = 0; i < locSubdDofs.size(); ++i)
          for (const auto& it : locSubdDofs[i])
            ddofs[i].insert(it-1);
      }

      return setupSolver(new ISTL::ASM(A, ddofs), op, solParams, solver);
    }
  } else if (prec == "amg" || prec == "gamg") {
    if (solParams.getBlock(block).getStringValue("multigrid_smoother") !=
        solParams.getBlock(block).getStringValue("multigrid_finesmoother")) {
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 3)
      return setupAMG2(solParams, block, op, solver);
#else
      std::cerr << "** ISTLSolParams ** Separate fine smoother not implemented." << std::endl;
      return nullptr;
#endif
    } else {
      std::string smoother = solParams.getBlock(block).getStringValue("multigrid_smoother");
      if (smoother.empty()) {
        std::cerr << "** ISTLSolParams ** No smoother defined for AMG, defaulting to ILU" << std::endl;
        smoother = "ilu";
      }
      if (smoother == "ilu")
        return setupAMG<ISTL::ILU0>(solParams, block, op, solver);
      else if (smoother == "sor")
        return setupAMG<ISTL::SOR>(solParams, block, op, solver);
      else if (smoother == "ssor")
        return setupAMG<ISTL::SSOR>(solParams, block, op, solver);
      else if (smoother == "jacobi")
        return setupAMG<ISTL::GJ>(solParams, block, op, solver);
    }
  }

  return nullptr;
}


std::tuple<std::unique_ptr<ISTL::InverseOperator>,
           std::unique_ptr<ISTL::Preconditioner>,
           std::unique_ptr<ISTL::Operator>> ISTLSolParams::setupPC(ISTL::Mat& A)
{
  std::unique_ptr<ISTL::InverseOperator> solver;
  std::unique_ptr<ISTL::Preconditioner> pre;
  std::unique_ptr<ISTL::Operator> op;

  if (solParams.getNoBlocks() > 2) {
    std::cerr << "*** ISTLSolParams ** More than two blocks are not implemented." << std::endl;
    return std::make_tuple(nullptr, nullptr, nullptr);
  }

  if (solParams.getNoBlocks() > 1) {
    op.reset(new ISTL::Operator(A));
    ISTL::BlockPreconditioner* bpre = new ISTL::BlockPreconditioner(A, adm.dd, solParams.getStringValue("schur"));
    pre.reset(bpre);
    for (size_t i = 0; i < adm.dd.getNoBlocks(); ++i)
      bpre->getBlockPre(i).reset(setupPCInternal(bpre->getBlock(i),
                                                 bpre->getBlockOp(i), i, nullptr));
    solver.reset(setupWithPreType(solParams, *op, *bpre));
  } else {
#ifdef HAVE_MPI
    if (adm.isParallel()) {
        Dune::IndexInfoFromGrid<int, int> index;
      for (size_t i = 0; i < adm.dd.getMLGEQ().size(); ++i)
        index.addLocalIndex(std::make_tuple(adm.dd.getGlobalEq(i+1), i, 0));
      Dune::OwnerOverlapCopyCommunication<int,int> comm(index, *adm.getCommunicator());
      //ISTL::ParMatrixAdapter* nop = new ISTL::ParMatrixAdapter(A, comm);
      //op.reset(nop);
      //pre.reset(setupPCInternal(A, *nop, 0, &solver));
    } else
#endif
    {
      ISTL::Operator* nop = new ISTL::Operator(A);
      op.reset(nop);
      pre.reset(setupPCInternal(A, *nop, 0, &solver));
    }
  }

  return std::make_tuple(std::move(solver), std::move(pre), std::move(op));
}
