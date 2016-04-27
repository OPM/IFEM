// $Id$
//==============================================================================
//!
//! \file LinSolParams.h
//!
//! \date Jan 29 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Linear solver parameters for PETSc matrices.
//! \details Includes linear solver method, preconditioner
//! and convergence criteria.
//!
//==============================================================================

#ifndef _LINSOLPARAMS_H
#define _LINSOLPARAMS_H

#include "PETScSupport.h"

#include <array>
#include <iostream>
#include <string>
#include <vector>

#ifdef HAS_PETSC

typedef std::vector<int>         IntVec;       //!< Integer vector
typedef std::vector<IntVec>      IntMat;       //!< Integer matrix
typedef std::vector<std::string> StringVec;    //!< String vector
typedef std::vector<StringVec>   StringMat;    //!< String matrix
typedef std::vector<IS>          ISVec;        //!< Index set vector
typedef std::vector<ISVec>       ISMat;        //!< Index set matrix

//! \brief Null-space for matrix operator
enum NullSpace { NONE, CONSTANT, RIGID_BODY };

//! \brief Schur preconditioner methods
enum SchurPrec { SIMPLE, MSIMPLER, PCD };

#endif

class TiXmlElement;


#ifdef HAS_PETSC
  #define BLANK_IF_NO_PETSC(a) a
#else
  #define BLANK_IF_NO_PETSC(a) ""
#endif


/*!
  \brief Class for linear solver parameters.
  \details It contains information about solver method, preconditioner
  and convergence criteria.
*/

class LinSolParams
{
public:
  //! \brief Default constructor.
  LinSolParams(PetscInt n, int m = 1) :
    atol(1e-6),
    rtol(1e-6),
    dtol(1e3),
    maxIts(1000),
    nResetSolver(0),
#ifdef HAS_PETSC
    schurPrec(SIMPLE),
#endif
    blocks(1),
    nsd(n),
    msgLev(m)
  {}

  //! \brief Read linear solver parameters from stream.
  bool read(std::istream& is, int nparams = 10);

  //! \brief Read linear solver parameters from XML document.
  bool read(const TiXmlElement* elem);

  //! \brief Read linear solver parameters from file.
  bool read(const char* filename);

  //! \brief Linear solver settings for a block of the linear system
  struct BlockParams {

    BlockParams() :
      package(BLANK_IF_NO_PETSC(MATSOLVERPETSC)),
      basis(1),
      comps(0),
      overlap(1),
      ilu_fill_level(0),
      mglevel(-1),
      noPreSmooth(-1),
      noFineSmooth(-1),
      noPostSmooth(-1),
      presmoother(BLANK_IF_NO_PETSC(PCILU)),
      postsmoother(BLANK_IF_NO_PETSC(PCILU)),
      mgKSP("defrichardson"),
      maxCoarseSize(-1),
      subdomains({{0,0,0}})
#ifdef HAS_PETSC
      , nullspace(NONE)
#endif
    {}

    //! \brief Settings for a directional smoother
    struct DirSmoother {
      int order; //!< Ordering of DOFs
      std::string type; //!< Directional smoother types

      DirSmoother(int o, const std::string& t) : order(o), type(t) {}
    };

    //! \brief Structure holding settings for a GAMG preconditioner
    struct GAMGSettings {
      //! \brief Default constructor
      GAMGSettings() :
        procEqLimit(-1),
        repartition(-1),
        useAggGasm(-1),
        reuseInterp(-1),
        threshold(-1)
      {}

      std::string type;     //!< Type for GAMG multigrid method ("geo" or "agg")
      PetscInt procEqLimit; //!< Limit number of equations on each process on coarse grid
      PetscInt repartition; //!< Repartition coarse grid for GAMG preconditioner
      PetscInt useAggGasm;  //!< Use aggregation aggregates for GASM smoother
      PetscInt reuseInterp; //!< Reuse interpolation
      PetscReal threshold;  //!< Smoother drop toleranse

      //! \brief Read linear solver parameters from XML document.
      void read(const TiXmlElement* elem);
    };

    //! \brief Structure holding settings for a ML preconditioner
    struct MLSettings {
      //! \brief Default constructor
      MLSettings() :
        symmetrize(-1),
        reuseInterp(-1),
        blockScaling(-1),
        putOnSingleProc(-1),
        threshold(-1),
        dampingFactor(-1),
        repartitionRatio(-1),
        reusable(-1),
        repartition(-1),
        keepAggInfo(-1),
        aux(-1),
        auxThreshold(-1)
      {}

      std::string coarsePackage; //!< Coarse matrix solver package
      std::string coarseSolver; //!< DD type coarse solver for ML
      std::string coarsenScheme; //!< Coarsening scheme for ML
      PetscInt symmetrize; //!< Symmetrize aggregation
      PetscInt reuseInterp; //!< Reuse interpolation operators between solves
      PetscInt blockScaling; //!< Scale all dofs at each node together
      PetscInt putOnSingleProc; //!< If below threshold, assign to a single processor
      PetscReal threshold; //!< Smoother drop toleranse for ML
      PetscReal dampingFactor; //!< Damping factor
      PetscInt repartitionRatio; //!< Max-min ratio for repartitioning
      PetscInt reusable; //!< Store intermediate data structures such that the MG hierarchy is reusable
      PetscInt repartition; //!< Repartitioning
      PetscInt keepAggInfo;  //!< Allow ML preconditioner to be reused
      PetscInt aux; //!< Use auxillary coordinate based laplacian for aggregation
      PetscReal auxThreshold; //!< Tolerance for auxillary laplacian aggregation

      //! \brief Read linear solver parameters from XML document.
      void read(const TiXmlElement* elem);
    };

    //! \brief Structure holding settings for a Hypre preconditioner
    struct HypreSettings {
      //! \brief Default constructor
      HypreSettings() :
        type("boomeramg"),
        noAggCoarse(-1),
        noPathAggCoarse(-1),
        truncation(-1),
        threshold(-1)
      {}

      std::string type;          //!< Type of hypre preconditioner
      PetscInt noAggCoarse;      //!< Number of levels of agressive coarsening
      PetscInt noPathAggCoarse;  //!< Number of paths for aggressive coarsening
      PetscReal truncation;      //!< Truncation factor for interpolation
      PetscReal threshold;       //!< Drop tolerance for Hypre
      std::string coarsenScheme; //!< Coarsening scheme for Hypre

      //! \brief Read linear solver parameters from XML document.
      void read(const TiXmlElement* elem);
    };

    //! \brief Read settings from XML block
    //! \param[in] child XML block
    bool read(const TiXmlElement* child);

    std::string prec; //!< Preconditioner for block
    std::string package; //!< Package providing solver
    size_t basis; //!< Basis for block
    size_t comps; //!< Components from basis (1, 2, 3, 12, 13, 23, 123, ..., 0 = all)
    PetscInt overlap; //!< Overlap in ASM
    PetscInt ilu_fill_level; //!< Fill level for ILU
    PetscInt mglevel; //!< Levels for multigrid
    PetscInt noPreSmooth; //!< Number of presmoothings for AMG
    PetscInt noFineSmooth; //!< Number of fine smoothings for AMG
    PetscInt noPostSmooth; //!< Number of postsmoothings for AMG

    std::string presmoother;  //!< Presmoother for AMG
    std::string postsmoother; //!< Postsmoother for AMG
    std::string finesmoother; //!< Smoother on finest grid
    std::string mgKSP; //!< KSP type on MG levels
    int maxCoarseSize; //!< Max number of DOFs for coarse AMG system
    GAMGSettings gamg;      //!< Settings for GAMG preconditioner
    MLSettings ml;          //!< Settings for ML preconditioner
    HypreSettings hypre;    //!< Settings for Hypre preconditioner
    std::array<size_t,3> subdomains; //!< Number of local subdomains in each parameter direction

    std::vector<DirSmoother> dirSmoother; //!< Direction smoothers
#ifdef HAS_PETSC
    NullSpace nullspace; //!< Nullspace for matrix
#endif
  };

  //! \brief Get linear solver method
  const std::string& getMethod() const { return method; }

  //! \brief Get absolute convergence tolerance
  Real getAbsTolerance() const { return atol; }

  //! \brief Get relative convergence tolerance
  Real getRelTolerance() const { return rtol; }

  //! \brief Get divergence tolerance
  Real getDivTolerance() const { return dtol; }

  //! \brief Get number of iterations before reset (GMRES)
  PetscInt getResetSolver() const { return nResetSolver; }

  //! \brief Get maximum number of iterations
  int getMaxIterations() const { return maxIts; }

  //! \brief Number of blocks in matrix system
  size_t getNoBlocks() const { return blocks.size(); }

  //! \brief Obtain settings for a given block
  const BlockParams& getBlock(size_t i) const { return blocks[i]; }

  //! \brief Get main preconditioner
  const char* getPreconditioner() const
  {
    return (prec == "asmlu" ? "asm" : prec.c_str());
  }

#ifdef HAS_PETSC
  //! \brief Get schur complement type
  SchurPrec getSchurType() const { return schurPrec; }
#endif

  //! \brief Get dimensionality of node coordinates
  size_t getDimension() const { return nsd; }

  //! \brief Get the message level
  int getMessageLevel() const { return msgLev; }

#ifdef HAS_PETSC
  //! \brief Set linear solver parameters for KSP object
  void setParams(KSP& ksp, PetscIntMat& locSubdDofs,
                 PetscIntMat& subdDofs, PetscRealVec& coords,
                 ISMat& dirIndexSet) const;

  //! \brief Set directional smoother
  //! \param[in] PC The preconditioner to add smoother for
  //! \param[in] P The preconditioner matrix
  //! \param[in] block The index of the block to add smoother to
  //! \param[in] dirIndexSet The index set for the smoother
  bool addDirSmoother(PC pc, Mat P, int block, ISMat& dirIndexSet) const;

  //! \brief Set ML options
  //! \param[in] prefix The prefix of the block to set parameters for
  //! \param[in] block The index of the block to set parameters for
  void setMLOptions(const std::string& prefix, int block) const;

  //! \brief Set GAMG options
  //! \param[in] prefix The prefix of the block to set parameters for
  //! \param[in] block The index of the block to set parameters for
  void setGAMGOptions(const std::string& prefix, size_t block) const;

  //! \brief Set Hypre options
  //! \param[in] prefix The prefix of the block to set parameters for
  //! \param[in] block The index of the block to set parameters for
  void setHypreOptions(const std::string& prefix, size_t block) const;

  //! \brief Setup the coarse solver in a multigrid
  //! \param[in] PC The preconditioner to set coarse solver for
  //! \param[in] prefix The prefix of the block to set parameters for
  //! \param[in] block The index of the block to set parameters for
  void setupCoarseSolver(PC& pc, const std::string& prefix, size_t block) const;

  //! \brief Setup the smoothers in a multigrid
  //! \param[in] PC The preconditioner to set coarse solver for
  //! \param[in] block The index of the  block to set parameters for
  //! \param[in] dirIndexSet The index set for direction smoothers
  //! \param[in] locSubdDofs Local subdomain DOFs for ASM preconditioners
  //! \param[in] subdDofs Subdomain DOFs for ASM preconditioners
  void setupSmoothers(PC& pc, size_t block, ISMat& dirIndexSet,
                      const PetscIntMat& locSubdDofs,
                      const PetscIntMat& subdDofs) const;

  //! \brief Setup an additive Schwarz preconditioner
  //! \param[in] PC The preconditioner to set coarse solver for
  //! \param[in] overlap The overlap
  //! \param[in] asmlu True to use LU subdomain solvers
  //! \param[in] locSubdDofs Local subdomain DOFs for ASM preconditioners
  //! \param[in] subdDofs Subdomain DOFs for ASM preconditioners
  //! \param[in] smoother True if this is a smoother in multigrid
  void setupAdditiveSchwarz(PC& pc, int overlap, bool asmlu,
                            const PetscIntMat& locSubdDofs,
                            const PetscIntMat& subdDofs, bool smoother) const;
#endif
private:
  PetscReal              atol;              // Absolute tolerance
  PetscReal              rtol;              // Relative tolerance
  PetscReal              dtol;              // Divergence tolerance
  PetscInt               maxIts;            // Maximum number of iterations
  PetscInt               nResetSolver;      // Number of linear solves before reset of KSP/PC
#ifdef HAS_PETSC
  SchurPrec              schurPrec;         // Preconditioner for Schur system
#endif
  std::string            method;            // Linear solver method
  std::string            prec;              // Preconditioner

  std::vector<BlockParams> blocks; //!< Parameters for each block
  PetscInt               nsd;               //!< Number of space dimensions
  int                    msgLev;            //!< Flag for extra output
};

#endif
