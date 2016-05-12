// $Id$
//==============================================================================
//!
//! \file PETScBlockMatrix.h
//!
//! \date Oct 20 2012
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Representation of the system matrix in block format with interface
//! to PETSc routines for assembling and solving linear equation systems.
//!
//==============================================================================

#ifndef _PETSC_BLOCK_MATRIX_H
#define _PETSC_BLOCK_MATRIX_H

#include "PETScMatrix.h"
#ifdef HAS_PETSC
#include "PETScPCProd.h"
#endif

typedef std::vector<int>         IntVec;        //!< 3D PETSc integer matrix
typedef std::vector<PetscIntMat> PetscIntMat3D; //!< 3D PETSc integer matrix


/*!
  \brief Class for representing the system matrix in PETSc block format.
  \details It is an interface to PETSc modules for assembling and solving
  linear systems of equations.
*/

class PETScBlockMatrix : public PETScMatrix
{
public:
#ifdef HAS_PETSC
  //! \brief Constructor.
  PETScBlockMatrix(const ProcessAdm& padm, const LinSolParams& spar);
  //! \brief Destructor
  virtual ~PETScBlockMatrix();
#else
  //! \brief Constructor.
 PETScBlockMatrix(const ProcessAdm& padm, const LinSolParams& spar) : PETScMatrix(padm,spar) {}
#endif

  //! \brief Returns the matrix type.
  virtual Type getType() const { return PETSCBLOCK; }

#ifdef HAS_PETSC
  //! \brief Creates a copy of the system matrix and returns a pointer to it.
  virtual SystemMatrix* copy() const { return new PETScBlockMatrix(*this); }

  //! \brief Initialize blocks.
  //! \param[in] nc Number of components for blocks
  void setBlocks(const IntVec& nc) { ncomps = nc; }

  //! \brief Initializes the element assembly process.
  //! \details Must be called once before the element assembly loop.
  //! The PETSc data structures are initialized and the all symbolic operations
  //! that are needed before the actual assembly can start are performed.
  //! \param[in] sam Auxiliary data describing the FE model topology, etc.
  virtual void initAssembly(const SAM& sam, bool);

  //! \brief Initializes the matrix to zero assuming it is properly dimensioned.
  virtual void init();

  //! \brief Adds an element matrix into the associated system matrix.
  //! \details When multi-point constraints are present, contributions from
  //! these are also added into the system right-hand-side vector.
  //! \param[in] eM  The element matrix
  //! \param[in] sam Auxiliary data describing the FE model topology,
  //!                nodal DOF status and constraint equations
  //! \param     B   The system right-hand-side vector
  //! \param[in] e   Identifier for the element that \a eM belongs to
  //! \return \e true on successful assembly, otherwise \e false
  virtual bool assemble(const Matrix& eM, const SAM& sam,
                        SystemVector& B, int e);

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param B Right-hand-side vector on input, solution vector on output
  //! \param[in] newLHS \e true if the left-hand-side matrix has been updated
  virtual bool solve(SystemVector& B, bool newLHS, Real*);

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param[in] B Right-hand-side vector
  //! \param[out] x Solution vector
  //! \param[in] newLHS \e true if the left-hand-side matrix has been updated
  virtual bool solve(const SystemVector& B, SystemVector& x, bool newLHS);

  //! \brief Returns matrix block (for assignment).
  virtual Mat& getMatrixBlock(size_t i, size_t j)
  { return matvec[i*nblocks+j]; }
  //! \brief Returns matrix block (for read access).
  virtual const Mat& getMatrixBlock(size_t i, size_t j) const
  { return matvec[i*nblocks+j]; }

protected:
  Mat    Sp;                          //!< Preconditioner for Schur block
  Vec    QpL;                         //!< Vector with lumoed pressure mass 
  PC     S, Fp;                       //!< Preconditioners for pressure-convection-diffusion pc
  PCProd*  pcprod;                      //!< PCD preconditioner
  size_t   nblocks;                     //!< Number of blocks
  IntVec   ncomps;                      //!< Number of components
  IS*      isvec;                       //!< Index set for blocks
  Mat*     matvec;                      //!< Vector of matrix blocks
  PetscIntMat3D locSubdDofsBlock; //!< Dofs for subdomains for block matrix
  PetscIntMat3D subdDofsBlock;    //!< Dofs for subdomains for block matrix

  //! \brief Auxiliary method for assembly of block matrices.
  void assemPETSc(const Matrix& eM, PETScVector& SV,
		  const std::vector<int>& meen, const int* meqn,
		  const int* mpmceq, const int* mmceq, const Real* ttcc);

  //! \brief Decomposes an element matrix into different blocks.
  //! \param[in] eM    Element system for coupled problem
  //! \param[in] l2g   Local to global mapping for element degrees of freedom
  //! \param[out] eMb  Element block matrices
  //! \param[out] l2gb Local to global mapping for block element DOFs
  void getBlockElmMatData(const Matrix& eM, const PetscIntVec& l2g,
			  std::vector< std::vector<Matrix> >& eMb,
			  PetscIntMat& l2gb) const;

  //! \brief Renumbers a given right-hand-side vector.
  void renumberRHS(const Vec& b, Vec& bnew, bool renum2block = true);

  //! \brief Initializes the block matrix.
  bool setParameters(PETScMatrix* P2 = nullptr, PETScVector* Pb = nullptr);
#endif
};

#endif
