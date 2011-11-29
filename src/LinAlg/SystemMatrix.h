// $Id$
//==============================================================================
//!
//! \file SystemMatrix.h
//!
//! \date Jan 4 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief General representation of system matrices and vectors.
//!
//==============================================================================

#ifndef _SYSTEM_MATRIX_H
#define _SYSTEM_MATRIX_H

#include "MatVec.h"

class SAM;
class LinSolParams;


/*!
  \brief Base class for representing a system vector on different formats.
*/

class SystemVector
{
public:
  //! \brief The available system vector formats.
  enum Type { STD = 0, PETSC = 1 };

  //! \brief Static method creating a vector of the given \a vectorType
  static SystemVector* create(Type vectorType = STD);

protected:
  //! \brief Default constructor.
  SystemVector() {}

public:
  //! \brief Empty destructor
  virtual ~SystemVector() {}

  //! \brief Returns the vector type.
  virtual Type getType() const = 0;

  //! \brief Creates a copy of the system vector and returns a pointer to it.
  virtual SystemVector* copy() const = 0;

  //! \brief Returns the dimension of the system vector.
  virtual size_t dim() const = 0;

  //! \brief Sets the dimension of the system vector.
  virtual void redim(size_t n) = 0;

  //! \brief Resize the vector to length \a n.
  virtual void resize(size_t n, bool = false) { this->redim(n); }

  //! \brief Checks if the vector is empty.
  virtual bool empty() const { return this->dim() == 0; }

  //! \brief Returns the size of the system vector.
  size_t size() const { return this->dim(); }

  //! \brief Access through pointer.
  virtual real* getPtr() = 0;
  //! \brief Reference through pointer.
  virtual const real* getRef() const = 0;

  //! \brief Restores the vector contents from an array.
  //! \details This method must be implemented only by subclasses for which the
  //! getPtr and getRef methods do not return a pointer to the actual internal
  //! memory segment containing the actual vector data.
  virtual void restore(const real*) {}

  //! \brief Initializes the vector assuming it is properly dimensioned.
  virtual void init(real value = real(0)) = 0;

  //! \brief Copies entries from input vector
  virtual void copy(const SystemVector& x);

  //! \brief Begins communication step needed in parallel vector assembly.
  virtual bool beginAssembly() { return true; }
  //! \brief Ends communication step needed in parallel vector assembly.
  virtual bool endAssembly() { return true; }

  //! \brief Multiplication with a scalar.
  virtual void mult(real alpha) = 0;

  //! \brief L1-norm of the vector.
  virtual real L1norm() const = 0;

  //! \brief L2-norm of the vector.
  virtual real L2norm() const = 0;

  //! \brief Linfinity-norm of the vector.
  virtual real Linfnorm() const = 0;

protected:
  //! \brief Writes the system vector to the given output stream.
  virtual std::ostream& write(std::ostream& os) const { return os; }

  //! \brief Global stream operator printing the vector contents.
  friend std::ostream& operator<<(std::ostream& os, const SystemVector& X)
  {
    return X.write(os);
  }
};


/*!
  \brief Standard system vector stored as a single continuous array.
*/

class StdVector : public SystemVector, public utl::vector<real>
{
public:
  //! \brief Constructor creating an empty vector.
  StdVector() {}
  //! \brief Constructor creating a vector of length \a n.
  StdVector(size_t n) : utl::vector<real>(n) {}
  //! \brief Constructor creating a vector from an array.
  StdVector(const real* values, size_t n) : utl::vector<real>(values,n) {}
  //! \brief Overloaded copy constructor.
  StdVector(const std::vector<real>& vec)
  { this->insert(this->end(),vec.begin(),vec.end()); }

  //! \brief Returns the vector type.
  virtual Type getType() const { return STD; }

  //! \brief Creates a copy of the system vector and returns a pointer to it.
  virtual SystemVector* copy() const { return new StdVector(*this); }

  //! \brief Checks if the vector is empty.
  virtual bool empty() const { return this->std::vector<real>::empty(); }

  //! \brief Returns the dimension of the system vector.
  virtual size_t dim() const { return this->std::vector<real>::size(); }

  //! \brief Sets the dimension of the system vector.
  virtual void redim(size_t n) { this->std::vector<real>::resize(n,real(0)); }

  //! \brief Resize the vector to length \a n.
  //! \details Will erase the previous content, but only if the size changed,
  //! unless \a forceClear is \e true.
  virtual void resize(size_t n, bool forceClear = false)
  { this->utl::vector<real>::resize(n,forceClear); }

  //! \brief Access through pointer.
  virtual real* getPtr() { return this->ptr(); }
  //! \brief Reference through pointer.
  virtual const real* getRef() const { return this->ptr(); }

  //! \brief Initializes the vector to a given scalar value.
  virtual void init(real value = real(0)) { this->fill(value); }

  //! \brief Begins communication step needed in parallel vector assembly.
  virtual bool beginAssembly() { return true; }
  //! \brief Ends communication step needed in parallel vector assembly.
  virtual bool endAssembly() { return true; }

  //! \brief Multiplication with a scalar.
  virtual void mult(real alpha) { this->operator*=(alpha); }

  //! \brief L1-norm of the vector.
  virtual real L1norm() const { return this->asum(); }

  //! \brief L2-norm of the vector.
  virtual real L2norm() const { return this->norm2(); }

  //! \brief Linfinity-norm of the vector.
  virtual real Linfnorm() const { size_t off = 0; return this->normInf(off); }

protected:
  //! \brief Writes the system vector to the given output stream.
  virtual std::ostream& write(std::ostream& os) const
  { return os << static_cast<const utl::vector<real>&>(*this); }
};


/*!
  \brief Base class for representing a system matrix on different formats.
  \details The purpose of this class is to define a clean interface for the
  system matrix and underlying equation solvers, such that the finite element
  implementation can be performed without paying attention to the actual solver
  being used.
*/

class SystemMatrix
{
public:
  //! \brief The available system matrix formats.
  enum Type { DENSE = 0, SPR = 1, SPARSE = 2, SAMG = 3, PETSC = 4 };

  //! \brief Static method creating a matrix of the given \a matrixType.
  static SystemMatrix* create(Type matrixType, int num_thread_SLU = 1);
  //! \brief Static method creating a matrix of the given \a matrixType.
  static SystemMatrix* create(Type matrixType, const LinSolParams& spar);

protected:
  //! \brief Default constructor.
  SystemMatrix() {}

public:
  //! \brief Empty destructor.
  virtual ~SystemMatrix() {}

  //! \brief Returns the matrix type.
  virtual Type getType() const = 0;

  //! \brief Creates a copy of the system matrix and returns a pointer to it.
  virtual SystemMatrix* copy() const = 0;

  //! \brief Checks if the matrix is empty.
  virtual bool empty() const { return this->dim(0) == 0; }

  //! \brief Returns the dimension of the system matrix.
  virtual size_t dim(int idim = 1) const = 0;

  //! \brief Initializes the element assembly process.
  //! \details Must be called once before the element assembly loop.
  //! \param[in] sam Auxilliary data describing the FE model topology, etc.
  virtual void initAssembly(const SAM& sam) = 0;

  //! \brief Initializes the matrix to zero assuming it is properly dimensioned.
  virtual void init() = 0;

  //! \brief Begins communication step needed in parallel matrix assembly.
  virtual bool beginAssembly() { return true; }
  //! \brief Ends communication step needed in parallel matrix assembly.
  virtual bool endAssembly() { return true; }

  //! \brief Adds an element stiffness matrix into the system stiffness matrix.
  //! \param[in] eM  The element stiffness matrix
  //! \param[in] sam Auxilliary data describing the FE model topology,
  //!                nodal DOF status and constraint equations
  //! \param[in] e   Identifier for the element that \a eM belongs to
  //! \return \e true on successful assembly, otherwise \e false
  virtual bool assemble(const Matrix& eM, const SAM& sam, int e) = 0;
  //! \brief Adds an element stiffness matrix into the system stiffness matrix.
  //! \details When multi-point constraints are present, contributions from
  //! these are also added into the system right-hand-side load vector.
  //! \param[in] eM  The element stiffness matrix
  //! \param[in] sam Auxilliary data describing the FE model topology,
  //!                nodal DOF status and constraint equations
  //! \param     B   The system right-hand-side load vector
  //! \param[in] e   Identifier for the element that \a eM belongs to
  //! \return \e true on successful assembly, otherwise \e false
  virtual bool assemble(const Matrix& eM, const SAM& sam,
			SystemVector& B, int e) = 0;

  //! \brief Augments a similar matrix symmetrically to the current matrix.
  virtual bool augment(const SystemMatrix&, size_t, size_t) { return false; }

  //! \brief Truncates all small matrix elements to zero.
  virtual bool truncate(real = real(1.0e-16)) { return false; }

  //! \brief Adds a matrix with similar structure to the current matrix.
  virtual bool add(const SystemMatrix&, real = real(1)) { return false; }

  //! \brief Adds a constant diagonal matrix to the current matrix.
  virtual bool add(real) { return false; }

  //! \brief Performs a matrix-vector multiplication.
  virtual bool multiply(const SystemVector&, SystemVector&) { return false; }

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param newLHS \e true if the left-hand-side matrix has been updated
  virtual bool solve(SystemVector&, bool newLHS = true) { return false; }

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param b Right-hand-side vector
  //! \param x Solution vector
  //! \param newLHS \e true if the left-hand-side matrix has been updated
  virtual bool solve(const SystemVector& b, SystemVector& x, bool newLHS = true);

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param b Right-hand-side vector
  //! \param P Preconditioning matrix (if different than system matrix)
  //! \param newLHS \e true if the left-hand-side matrix has been updated
  virtual bool solve(SystemVector& b, SystemMatrix& P, bool newLHS = true)
  { return false; }

  //! \brief Returns the L-infinity norm of the matrix.
  virtual real Linfnorm() const = 0;

protected:
  //! \brief Writes the system matrix to the given output stream.
  virtual std::ostream& write(std::ostream& os) const { return os; }

  //! \brief Global stream operator printing the matrix contents.
  friend std::ostream& operator<<(std::ostream& os, const SystemMatrix& A)
  {
    return A.write(os);
  }
};

#endif
