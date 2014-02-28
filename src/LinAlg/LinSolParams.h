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

#include <iostream>
#include <string>
#include <vector>
#include "PETScSupport.h"

#ifdef HAS_PETSC

typedef std::vector<int>         IntVec;       //!< Integer vector
typedef std::vector<IntVec>      IntMat;       //!< Integer matrix
typedef std::vector<std::string> StringVec;    //!< String vector
typedef std::vector<StringVec>   StringMat;    //!< String matrix
typedef std::vector<bool>        BoolVec;      //!< Boolean vector
typedef std::vector<PetscBool>   PetscBoolVec; //!< Petsc boolean vector
typedef std::vector<IS>          ISVec;        //!< Index set vector
typedef std::vector<ISVec>       ISMat;        //!< Index set matrix

//! \brief Null-space for matrix operator
enum NullSpace { NONE, CONSTANT, RIGID_BODY };

//! \brief Schur preconditioner methods
enum SchurPrec { SIMPLE, MSIMPLER, PCD };

#endif

class TiXmlElement;


/*!
  \brief Class for linear solver parameters.
  \details It contains information about solver method, preconditioner
  and convergence criteria.
*/

class LinSolParams
{
public:
  //! \brief Default constructor.
  LinSolParams() { this->setDefault(); }
  //! \brief Copy constructor.
  LinSolParams(const LinSolParams& spar) { this->copy(spar); }
  
  //! \brief Virtual destructor
  ~LinSolParams() {}
  

  //! \brief Set default values.
  void setDefault();

  //! \brief Copy linear solver parameters.
  void copy(const LinSolParams& spar);

  //! \brief Read linear solver parameters from stream.
  bool read(std::istream& is, int nparams = 10);

  //! \brief Read linear solver parameters from XML document.
  bool read(const TiXmlElement* elem);

  //! \brief Read linear solver parameters from file.
  bool read(const char* filename);

#ifdef HAS_PETSC
  //! \brief Set linear solver method
  void setMethod(KSPType type) { method.assign(type); }

  //! \brief Set preconditioner
  void setPreconditioner(PCType type) { prec.assign(type); }

  //! \brief Set preconditioner for sub-block
  void setSubPreconditioner(PCType type, size_t i = 0) { subprec[i].assign(type); }

  //! \brief Set linear solver package
  void setPackage(MatSolverPackage stype, size_t i = 0) { package[i].assign(stype); }

  //! \brief Set absolute convergence tolerance
  void setAbsTolerance(Real eps) { atol = eps; }

  //! \brief Set relative convergence tolerance
  void setRelTolerance(Real eps) { rtol = eps; }

  //! \brief Set divergence tolerance
  void setDivTolerance(Real eps) { dtol = eps; }

  //! \brief Set maximum number of iterations
  void setMaxIterations(int its) { maxIts = its; }

  //! \brief Set number of overlap
  void setOverlap(int olap, size_t i = 0) { overlap[i] = olap; }

  //! \brief Set number of local subdomains for each patch
  void setLocalPartitioning(size_t NX = 0, size_t NY = 0, size_t NZ = 0, size_t i = 0)
  { nx[i] = NX; ny[i] = NY; nz[i] = NZ; }

  //! \brief Set null-space of matrix
  void setNullSpace(NullSpace nspc, size_t i = 0) { nullspc[i] = nspc; }

  //! \brief Get linear solver method
  const char* getMethod() const { return method.c_str(); }

  //! \brief Get preconditioner
  const char* getPreconditioner(size_t i = 0) const { return prec.c_str(); }

  //! \brief Get preconditioner
  const char* getSubPreconditioner(size_t i = 0) const { return subprec[i].c_str(); }

  //! \brief Get linear solver package
  const char* getPackage(size_t i = 0) const { return package[i].c_str(); }

  //! \brief Get absolute convergence tolerance
  Real getAbsTolerance() const { return atol; }

  //! \brief Get relative convergence tolerance
  Real getRelTolerance() const { return rtol; }

  //! \brief Get divergence tolerance
  Real getDivTolerance() const { return dtol; }

  //! \brief Get maximum number of iterations
  int getMaxIterations() const { return maxIts; }

  //! \brief Get number of overlaps
  int getOverlap(int i = 0) const { return overlap[i]; }

  //! \brief Get local partitioning
  int getLocalPartitioning(size_t dir = 0, size_t i = 0) const
  {
    switch (dir)
    {
      case 0:
	return nx[i];
      case 1:
	return ny[i];
      case 2:
	return nz[i];
      default:
	return 0;
    }
  }

  //! \brief Number of blocks in matrix system
  int getNoBlocks() const { return nblock; }

  //! \brief Number of components in a matrix block
  const std::vector<int>& getComponents() const { return ncomps; }

  //! \brief Get number of overlaps
  NullSpace getNullSpace(size_t i = 0) const { return nullspc[i]; }

  //! \brief Set linear solver parameters for KSP object
  void setParams(KSP& ksp, PetscIntMat& locSubdDofs,
                 PetscIntMat& subdDofs, PetscRealVec& coords,
                 PetscInt nsd, ISMat& dirIndexSet) const;

  //! \brief Set directional smoother
  bool addDirSmoother(PC pc, Mat P, ISMat& dirIndexSet) const;

private:
  int                    nLinSolves;        // Counter for linear solver

  PetscReal              atol;              // Absolute tolerance
  PetscReal              rtol;              // Relative tolerance
  PetscReal              dtol;              // Divergence tolerance
  PetscInt               maxIts;            // Maximum number of iterations
  PetscInt               nResetSolver;      // Number of linear solves before reset of KSP/PC
  int                    nblock;            // Number of block
  bool                   schur;             // Schur complement solver
  SchurPrec              schurPrec;         // Preconditioner for Schur system
  std::string            method;            // Linear solver method
  std::string            prec;              // Preconditioner
  StringVec              subprec;           // Preconditioners for block-system
  StringVec              hypretype;         // Type of hypre preconditioner
  StringVec              package;           // Linear software package (petsc, superlu_dist, ...)
  BoolVec                asmlu;             // Use lu as subdomain solver
  IntVec                 ncomps;            // Components for each fields in block-vector
  PetscIntVec            overlap;           // Number of overlap in ASM
  PetscIntVec            levels;            // Number of levels of fill to use
  PetscIntVec            mglevels;          // Number of levels for MG
  PetscIntVec            noPreSmooth;       // Number of presmoothings for AMG
  PetscIntVec            noPostSmooth;      // Number of postsmoothings for AMG
  PetscIntVec            noFineSmooth;      // Number of fine grid smoothings for AMG
  StringVec              presmoother;       // Presmoother for AMG
  StringVec              postsmoother;      // Postsmoother for AMG
  StringVec              finesmoother;      // Smoother on finest grid
  IntVec                 maxCoarseSize;     // Max number of DOFS for coarse AMG system
  StringVec              MLCoarsePackage;   // Coarse matrix solver package
  StringVec              MLCoarseSolver;    // DD type coarse solver for ML
  StringVec              MLCoarsenScheme;   // Coarsening scheme for ML
  PetscIntVec            MLSymmetrize;      // Symmetrize aggregation
  PetscIntVec            MLReuseInterp;     // Reuse interpolation operators between solves
  PetscIntVec            MLBlockScaling;    // Scale all dofs at each node together
  PetscIntVec            MLPutOnSingleProc; // If below assign to a sigle processor 
  PetscRealVec           MLThreshold;       // Smoother drop toleranse for ML
  PetscRealVec           MLDampingFactor;   // Damping factor
  PetscRealVec           MLRepartitionRatio;// Max-min ratio for repartitioning
  PetscIntVec            MLRepartition;     // Repartitioning
  PetscIntVec            MLReusable;        // Store intermediate datastructures such that the MG hierarchy is reusable
  PetscIntVec            MLKeepAggInfo;     // Allow ML preconditioner to be reused
  PetscIntVec            MLAux;             // Use auxillary coordinate based laplacian for aggregation
  StringVec              GAMGtype;          // Type for GAMG multigrid method ("geo" or "agg") ("geo")
  PetscIntVec            GAMGprocEqLimit;   // Limit number of equation on each process on coarse grid
  PetscIntVec            GAMGrepartition;   // Repartition coarse grid for GAMG preconditioner (0 = false, 1 = true)
  PetscIntVec            GAMGuseAggGasm;    // Use aggregation aggregates for GASM smoother     (0 = false, 1 = true)
  PetscIntVec            GAMGreuseInterp;   // Reuse interpolation
  PetscIntVec            HypreNoAggCoarse; // Number of levels of aggressive coarsening
  PetscIntVec            HypreNoPathAggCoarse; // Number of paths for aggressive coarsening
  PetscRealVec           HypreTruncation;   // Truncation factor for interpolation
  PetscRealVec           HypreThreshold;    // Drop tolerance for Hypre
  StringVec              HypreCoarsenScheme;// Coarsening scheme for Hypre
  IntVec                 nx;                // Number of local subdomains in first parameter direction
  IntVec                 ny;                // Number of local subdomains in second parameter direction
  IntVec                 nz;                // Number of local subdomains in third parameter direction
  StringMat              dirsmoother;       // Directional smoother types
  IntMat                 dirOrder;          // Direction smoother orders/numberings
  std::vector<NullSpace> nullspc;           // Null-space for matrix

  friend class PETScMatrix;
  friend class PETScBlockMatrix;
#endif // HAS_PETSC
};

#endif
