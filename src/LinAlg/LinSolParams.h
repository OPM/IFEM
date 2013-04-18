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
#ifdef HAS_PETSC
#include <string>
#include <vector>
#include "petscksp.h"
#endif


class TiXmlElement;


//! \brief Null-space for matrix operator
enum NullSpace { NONE, CONSTANT, RIGID_BODY };

//! \brief Schur preconditioner methods
enum SchurPrec { SIMPLE, MSIMPLER, PCD };

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
  virtual void setParams(KSP& ksp, std::vector<std::vector<PetscInt> >& locSubdDofs,
			 std::vector<std::vector<PetscInt> >& subdDofs) const;

private:
  PetscReal atol;          // Absolute tolerance
  PetscReal rtol;          // Relative tolerance
  PetscReal dtol;          // Divergence tolerance
  PetscInt  maxIts;        // Maximum number of iterations
  int       nblock;        // Number of block
  bool      schur;         // Schur complement solver 
  SchurPrec schurPrec;     // Preconditioner for Schur system 
  std::string                method;           // Linear solver method
  std::string                prec;             // Preconditioner
  std::vector<std::string>   subprec;          // Preconditioners for block-system
  std::vector<std::string>   hypretype;        // Type of hypre preconditioner
  std::vector<std::string>   package;          // Linear software package (petsc, superlu_dist, ...)
  std::vector<bool>          asmlu;            // Use lu as subdomain solver
  std::vector<int>           ncomps;           // Components for each fields in block-vector
  std::vector<PetscInt>      overlap;          // Number of overlap in ASM
  std::vector<PetscInt>      levels;           // Number of levels of fill to use
  std::vector<PetscInt>      mglevels;         // Number of levels for MG
  std::vector<PetscInt>      noPreSmooth;      // Number of presmoothings for AMG
  std::vector<PetscInt>      noPostSmooth;     // Number of postsmoothings for AMG
  std::vector<PetscInt>      noFineSmooth;     // Number of fine grid smoothings for AMG
  std::vector<std::string>   presmoother;      // Presmoother for AMG
  std::vector<std::string>   postsmoother;     // Postsmoother for AMG
  std::vector<std::string>   finesmoother;     // Smoother on finest grid
  std::vector<int>           maxCoarseSize;    // Max number of DOFS for coarse AMG system
  std::vector<PetscInt>      MLCoarsenScheme;  // Coarsening scheme for ML
  std::vector<PetscReal>     MLThreshold;      // Smoother drop toleranse for ML
  std::vector<PetscReal>     MLDampingFactor;  // Damping factor
  std::vector<int>           nx;               // Number of local subdomains in first parameter direction
  std::vector<int>           ny;               // Number of local subdomains in second parameter direction
  std::vector<int>           nz;               // Number of local subdomains in third parameter direction
  std::vector<NullSpace>     nullspc;           // Null-space for matrix

  friend class PETScMatrix;
  friend class PETScBlockMatrix;
#endif // HAS_PETSC
};

#endif
