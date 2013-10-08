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

#ifndef _LINSOLPARAMS_MATRIX_H
#define _LINSOLPARAMS_MATRIX_H

#include <iostream>
#include <string>
#include <vector>
#ifdef HAS_PETSC
#include "PETScSupport.h"

typedef std::vector<PetscInt>    PetscIntVec;  //!< PETSc integer vector
typedef std::vector<PetscIntVec> PetscIntMat;  //!< PETSc integer matrix
typedef std::vector<PetscReal>   PetscRealVec; //!< PETSc real vector
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

  TODO: Why are all the member functions declared virtual when no sub-class?
  What a waste of vtable's. Consider removing.
*/

class LinSolParams
{
public:
  //! \brief Default constructor.
  LinSolParams() { this->setDefault(); }
  //! \brief Copy constructor.
  LinSolParams(const LinSolParams& spar) { this->copy(spar); }

  //! \brief Set default values.
  virtual void setDefault();

  //! \brief Copy linear solver parameters.
  virtual void copy(const LinSolParams& spar);

  //! \brief Read linear solver parameters from stream.
  virtual bool read(std::istream& is, int nparams = 10);

  //! \brief Read linear solver parameters from XML document.
  virtual bool read(const TiXmlElement* elem);

  //! \brief Read linear solver parameters from file.
  virtual bool read(const char* filename);

#ifdef HAS_PETSC
  //! \brief Set linear solver method
  virtual void setMethod(KSPType type) { method.assign(type); }

  //! \brief Set preconditioner
  virtual void setPreconditioner(PCType type) { prec.assign(type); }

  //! \brief Set preconditioner for sub-block
  virtual void setSubPreconditioner(PCType type, size_t i = 0) { subprec[i].assign(type); }

  //! \brief Set linear solver package
  virtual void setPackage(MatSolverPackage stype, size_t i = 0) { package[i].assign(stype); }

  //! \brief Set absolute convergence tolerance
  virtual void setAbsTolerance(Real eps) { atol = eps; }

  //! \brief Set relative convergence tolerance
  virtual void setRelTolerance(Real eps) { rtol = eps; }

  //! \brief Set divergence tolerance
  virtual void setDivTolerance(Real eps) { dtol = eps; }

  //! \brief Set maximum number of iterations
  virtual void setMaxIterations(int its) { maxIts = its; }

  //! \brief Set number of overlap
  virtual void setOverlap(int olap, size_t i = 0) { overlap[i] = olap; }

  //! \brief Set number of local subdomains for each patch
  virtual void setLocalPartitioning(size_t NX = 0, size_t NY = 0, size_t NZ = 0, size_t i = 0)
  { nx[i] = NX; ny[i] = NY; nz[i] = NZ; }

  //! \brief Set null-space of matrix
  virtual void setNullSpace(NullSpace nspc, size_t i = 0) { nullspc[i] = nspc; }

  //! \brief Get linear solver method
  virtual const char* getMethod() const { return method.c_str(); }

  //! \brief Get preconditioner
  virtual const char* getPreconditioner(size_t i = 0) const { return prec.c_str(); }

  //! \brief Get preconditioner
  virtual const char* getSubPreconditioner(size_t i = 0) const { return subprec[i].c_str(); }

  //! \brief Get linear solver package
  virtual const char* getPackage(size_t i = 0) const { return package[i].c_str(); }

  //! \brief Get absolute convergence tolerance
  virtual Real getAbsTolerance() const { return atol; }

  //! \brief Get relative convergence tolerance
  virtual Real getRelTolerance() const { return rtol; }

  //! \brief Get divergence tolerance
  virtual Real getDivTolerance() const { return dtol; }

  //! \brief Get maximum number of iterations
  virtual int getMaxIterations() const { return maxIts; }

  //! \brief Get number of overlaps
  virtual int getOverlap(int i = 0) const { return overlap[i]; }

  //! \brief Get local partitioning
  virtual int getLocalPartitioning(size_t dir = 0, size_t i = 0) const
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
  virtual int getNoBlocks() const { return nblock; }

  //! \brief Number of components in a matrix block
  virtual const std::vector<int>& getComponents() const { return ncomps; }

  //! \brief Get number of overlaps
  virtual NullSpace getNullSpace(size_t i = 0) const { return nullspc[i]; }

  //! \brief Set linear solver parameters for KSP object
  virtual void setParams(KSP& ksp, PetscIntMat& locSubdDofs,
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
  PetscIntVec            MLCoarsenScheme;   // Coarsening scheme for ML
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
