// $Id$
//==============================================================================
//!
//! \file SparseMatrix.h
//!
//! \date Jan 8 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of the system matrix on an unstructured sparse format.
//!
//==============================================================================

#ifndef _SPARSE_MATRIX_H
#define _SPARSE_MATRIX_H

#include "SystemMatrix.h"
#include <iostream>
#include <map>
#include <set>

typedef std::vector<int>         IntVec;    //!< General integer vector
typedef std::pair<size_t,size_t> IJPair;    //!< 1-based matrix indices
typedef std::map<IJPair,Real>    ValueMap;  //!< Index to matrix value mapping
typedef ValueMap::const_iterator ValueIter; //!< Iterator over matrix elements

struct SuperLUdata;


/*!
  \brief Class for representing a system matrix on an unstructured sparse form.
  \details The sparse matrix is editable in the sense that non-zero entries may
  be added at arbitrary locations. The class comes with methods for solving a
  linear system of equations based on the current matrix and a given RHS-vector,
  using either the commercial SAMG package or the public domain SuperLU package.
*/

class SparseMatrix : public SystemMatrix
{
public:
  //! \brief Available equation solvers for this matrix type.
  enum SparseSolver { NONE, SUPERLU, S_A_M_G, UMFPACK };

  //! \brief Default constructor creating an empty matrix.
  SparseMatrix(SparseSolver eqSolver = NONE, int nt = 1);
  //! \brief Constructor creating a \f$m \times n\f$ matrix.
  SparseMatrix(size_t m, size_t n = 0);
  //! \brief Copy constructor.
  SparseMatrix(const SparseMatrix& B);
  //! \brief The destructor frees the dynamically allocated arrays.
  virtual ~SparseMatrix();

  //! \brief Returns the matrix type.
  virtual Type getType() const { return solver == S_A_M_G ? SAMG : SPARSE; }

  //! \brief Creates a copy of the system matrix and returns a pointer to it.
  virtual SystemMatrix* copy() const { return new SparseMatrix(*this); }

  //! \brief Locks or unlocks the sparsity pattern.
  //! \param[in] doLock If \e true, lock pattern, otherwise unlock it
  //! \return \e true, if the pattern already was locked
  virtual bool lockPattern(bool doLock);

  //! \brief Resizes the matrix to dimension \f$r \times c\f$.
  //! \details Will erase previous content, also if the size is unchanged.
  //! If the size is not changed, the sparsity pattern stored in the members
  //! \a IA and \a JA is retained such that subsequent matrix assembly steps
  //! may be performed directly into the optimized storage format, unless
  //! \a forceEditable is \e true.
  void resize(size_t r, size_t c = 0, bool forceEditable = false);

  //! \brief Resizes the matrix to dimension \f$r \times c\f$.
  //! \details Will preserve existing matrix content within the new dimension.
  bool redim(size_t r, size_t c);

  //! \brief Query number of matrix rows.
  size_t rows() const { return nrow; }
  //! \brief Query number of matrix columns.
  size_t cols() const { return ncol; }
  //! \brief Query total matrix size in terms of number of non-zero elements.
  size_t size() const { return editable ? elem.size() : A.size(); }

  //! \brief Returns the dimension of the system matrix.
  //! \param[in] idim Which direction to return the dimension in.
  virtual size_t dim(int idim = 1) const;

  //! \brief Index-1 based element access.
  Real& operator()(size_t r, size_t c);
  //! \brief Index-1 based element reference.
  const Real& operator()(size_t r, size_t c) const;

  //! \brief For traversal of the non-zero elements of an editable matrix.
  const ValueMap& getValues() const { return elem; }

  //! \brief Print sparsity pattern - for inspection purposes.
  void printSparsity(std::ostream& os) const;

  //! \brief Print the matrix in full rectangular form.
  //! \warning Not recommended for matrices of nontrivial size.
  void printFull(std::ostream& os) const;

  //! \brief Dumps the system matrix on a specified format.
  virtual void dump(std::ostream&, char, const char* = nullptr);

  //! \brief Initializes the element assembly process.
  //! \details Must be called once before the element assembly loop.
  //! \param[in] sam Auxiliary data describing the FE model topology, etc.
  //! \param[in] delayLocking If \e true, do not lock the sparsity pattern yet
  virtual void initAssembly(const SAM& sam, bool delayLocking);

  //! \brief Initializes the element sparsity pattern based on node connections.
  //! \param[in] sam Auxiliary data describing the FE model topology, etc.
  //! \param[in] delayLocking If \e true, do not lock the sparsity pattern yet
  void preAssemble(const SAM& sam, bool delayLocking);

  //! \brief Initializes the element sparsity pattern based on node connections.
  //! \param[in] MMNPC Matrix of matrices of nodal point correspondances
  //! \param[in] nel Number of elements
  void preAssemble(const std::vector<IntVec>& MMNPC, size_t nel);

  //! \brief Initializes the matrix to zero assuming it is properly dimensioned.
  virtual void init();

  //! \brief Adds an element matrix into the associated system matrix.
  //! \param[in] eM  The element matrix
  //! \param[in] sam Auxiliary data describing the FE model topology,
  //!                nodal DOF status and constraint equations
  //! \param[in] e   Identifier for the element that \a eM belongs to
  //! \return \e true on successful assembly, otherwise \e false
  virtual bool assemble(const Matrix& eM, const SAM& sam, int e);
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
  //! \brief Adds an element matrix into the associated system matrix.
  //! \details When multi-point constraints are present, contributions from
  //! these are also added into the system right-hand-side vector.
  //! \param[in] eM   The element matrix
  //! \param[in] sam  Auxiliary data describing the FE model topology,
  //!                 nodal DOF status and constraint equations
  //! \param     B    The system right-hand-side vector
  //! \param[in] meen Matrix of element equation numbers
  //! \return \e true on successful assembly, otherwise \e false
  virtual bool assemble(const Matrix& eM, const SAM& sam,
			SystemVector& B, const IntVec& meen);

  //! \brief Adds a nodal vector into columns of a non-symmetric sparse matrix.
  //! \param[in] V   The nodal vector
  //! \param[in] sam Auxiliary data describing the FE model topology,
  //!                nodal DOF status and constraint equations
  //! \param[in] n   Identifier for the node that \a V belongs to
  //! \param[in] col Index of first column which should receive contributions
  //! \return \e true on successful assembly, otherwise \e false
  //!
  //! \details This method can be used for rectangular matrices whose rows
  //! correspond to the equation ordering og the provided \a sam object.
  bool assembleCol(const RealArray& V, const SAM& sam, int n, size_t col);

  //! \brief Adds a scalar value into columns of a non-symmetric sparse matrix.
  //! \param[in] val The value to add for each DOF of the specified node
  //! \param[in] sam Auxiliary data describing the FE model topology,
  //!                nodal DOF status and constraint equations
  //! \param[in] n   Identifier for the node that \a val belongs to
  //! \param[in] col Index of first column which should receive contributions
  //! \return \e true on successful assembly, otherwise \e false
  //!
  //! \details This method can be used for rectangular matrices whose rows
  //! correspond to the equation ordering og the provided \a sam object.
  bool assembleCol(Real val, const SAM& sam, int n, size_t col)
  {
    return this->assembleCol(RealArray(1,val),sam,n,col);
  }

  //! \brief Augments a similar matrix symmetrically to the current matrix.
  //! \param[in] B  The matrix to be augmented
  //! \param[in] r0 Row offset for the augmented matrix
  //! \param[in] c0 Column offset for the augmented matrix
  virtual bool augment(const SystemMatrix& B, size_t r0, size_t c0);

  //! \brief Truncates all small matrix elements to zero.
  //! \param[in] threshold Zero tolerance relative to largest diagonal element
  virtual bool truncate(Real threshold = Real(1.0e-16));

  //! \brief Adds a matrix with similar sparsity pattern to the current matrix.
  //! \param[in] B     The matrix to be added
  //! \param[in] alpha Scale factor for matrix \b B
  virtual bool add(const SystemMatrix& B, Real alpha = Real(1));

  //! \brief Adds the diagonal matrix \f$\sigma\f$\b I to the current matrix.
  virtual bool add(Real sigma);

  //! \brief Performs the matrix-vector multiplication \b C = \a *this * \b B.
  virtual bool multiply(const SystemVector& B, SystemVector& C) const;

  using SystemMatrix::solve;
  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param B Right-hand-side vector on input, solution vector on output
  //! \param[in] newLHS \e true if the left-hand-side matrix has been updated
  //! \param[out] rc Reciprocal condition number of the LHS-matrix (optional)
  virtual bool solve(SystemVector& B, bool newLHS = true, Real* rc = nullptr);

  //! \brief Calculate compressed-sparse-row arrays from element map.
  //! \param[out] IA Start index of each row in JA
  //! \param[out] JA Column indices
  //! \param[in] nrow Number of rows
  //! \param[in] elem Matrix element indices
  static void calcCSR(IntVec& IA, IntVec& JA,
                      size_t nrow, const ValueMap& elem);

protected:
  //! \brief Converts the matrix to an optimized row-oriented format.
  //! \details The optimized format is suitable for the SAMG equation solver.
  bool optimiseSAMG(bool transposed = false);

  //! \brief Converts the matrix to an optimized column-oriented format.
  //! \details The optimized format is suitable for the SuperLU equation solver.
  bool optimiseSLU();

  //! \brief Converts the matrix to an optimized column-oriented format.
  //! \param[in] dofc Set of free DOFs coupled to each free DOF
  //!
  //! \details The optimized format is suitable for the SuperLU equation solver.
  bool optimiseSLU(const std::vector< std::set<int> >& dofc);

  //! \brief Invokes the SAMG equation solver for a given right-hand-side.
  //! \param B Right-hand-side vector on input, solution vector on output
  bool solveSAMG(Vector& B);

  //! \brief Invokes the SuperLU equation solver for a given right-hand-side.
  //! \details This method uses the simple driver \a dgssv.
  //! \param B Right-hand-side vector on input, solution vector on output
  bool solveSLU(Vector& B);

  //! \brief Invokes the SuperLU equation solver for a given right-hand-side.
  //! \details This method uses the expert driver \a dgssvx.
  //! \param B Right-hand-side vector on input, solution vector on output
  //! \param[out] rcond Reciprocal condition number of the LHS-matrix (optional)
  bool solveSLUx(Vector& B, Real* rcond);

  //! \brief Invokes the UMFPACK equation solver for a given right-hand-side.
  //! \param B Right-hand-side vector on input, solution vector on output
  //! \param[out] rcond Reciprocal condition number of the LHS-matrix (optional)
  bool solveUMF(Vector& B, Real* rcond);

  //! \brief Writes the system matrix to the given output stream.
  virtual std::ostream& write(std::ostream& os) const;

  //! \brief Returns the L-infinity norm of the matrix.
  virtual Real Linfnorm() const;

public:
  static bool printSLUstat; //!< Print solution statistics for SuperLU?

private:
  //! Flag for the editability of the matrix elements:
  //!  'V' : Values may be edited but the pattern is temporarily locked
  //!  'P' : Both values and pattern may be edited
  //! '\0' : Values may be edited but the pattern is permanently locked
  char editable;

  bool factored; //!< \e true when the matrix is factorized
  size_t nrow;   //!< Number of matrix rows
  size_t ncol;   //!< Number of matrix columns

  ValueMap       elem; //!< Stores nonzero matrix elements with index pairs
  SparseSolver solver; //!< Which equation solver to use
  SuperLUdata*    slu; //!< Matrix data for the SuperLU equation solver
  int      numThreads; //!< Number of threads to use for the SuperLU_MT solver

#ifdef HAS_UMFPACK
  void* umfSymbolic; //!< Symbolically factored matrix for UMFPACK
#endif

protected:
  IntVec IA; //!< Identifies the beginning of each row or column
  IntVec JA; //!< Specifies column/row index of each nonzero element
  Vector  A; //!< Stores the nonzero matrix elements
};

#endif
