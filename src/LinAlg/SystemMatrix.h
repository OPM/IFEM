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
#include "LinAlgenums.h"

class SAM;
class LinSolParams;
class ProcessAdm;


/*!
  \brief Base class for representing a system vector on different formats.
*/

class SystemVector
{
public:
  //! \brief Static method creating a vector of the given type.
  static SystemVector* create(const ProcessAdm* adm, LinAlg::MatrixType vtype);

protected:
  //! \brief The default constructor is protected to allow sub-classes only.
  SystemVector() {}

public:
  //! \brief Empty destructor.
  virtual ~SystemVector() {}

  //! \brief Returns the vector type.
  virtual LinAlg::MatrixType getType() const = 0;

  //! \brief Creates a copy of the system vector and returns a pointer to it.
  virtual SystemVector* copy() const = 0;

  //! \brief Returns the dimension/size of the system vector.
  virtual size_t dim() const = 0;

  //! \brief Sets the dimension of the system vector.
  virtual void redim(size_t n) = 0;

  //! \brief Resizes the vector to length \a n.
  virtual void resize(size_t n, bool = false) { this->redim(n); }

  //! \brief Checks if the vector is empty.
  bool empty() const { return this->dim() == 0; }

  //! \brief Access through pointer.
  virtual Real* getPtr() = 0;
  //! \brief Reference through pointer.
  virtual const Real* getRef() const = 0;

  //! \brief Initializes the vector assuming it is properly dimensioned.
  virtual void init(Real value = Real(0)) = 0;

  //! \brief Copies entries from input vector \b x into \a *this.
  SystemVector& copy(const SystemVector& x);

  //! \brief Begins communication step needed in parallel vector assembly.
  virtual bool beginAssembly() { return true; }
  //! \brief Ends communication step needed in parallel vector assembly.
  virtual bool endAssembly() { return true; }

  //! \brief Multiplication with a scalar.
  virtual void mult(Real alpha) = 0;

  //! \brief Addition of another system vector to this one.
  virtual void add(const SystemVector& vec, Real scale = Real(1)) = 0;

  //! \brief L1-norm of the vector.
  virtual Real L1norm() const = 0;

  //! \brief L2-norm of the vector.
  virtual Real L2norm() const = 0;

  //! \brief Linfinity-norm of the vector.
  virtual Real Linfnorm() const = 0;

  //! \brief Dumps the system vector on a specified format.
  virtual void dump(std::ostream&, LinAlg::StorageFormat,
                    const char* = nullptr) const {}

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

class StdVector : public SystemVector, public utl::vector<Real>
{
public:
  //! \brief Default constructor creating a vector of length \a n.
  explicit StdVector(size_t n = 0) : utl::vector<Real>(n) {}
  //! \brief Constructor creating a vector from an array.
  StdVector(const Real* values, size_t n) : utl::vector<Real>(values,n) {}
  //! \brief Overloaded copy constructor.
  explicit StdVector(const std::vector<Real>& vec) : utl::vector<Real>(vec) {}

  //! \brief Returns the vector type.
  virtual LinAlg::MatrixType getType() const { return LinAlg::DENSE; }

  //! \brief Creates a copy of the system vector and returns a pointer to it.
  virtual SystemVector* copy() const { return new StdVector(*this); }

  //! \brief Returns the dimension of the system vector.
  virtual size_t dim() const { return this->size(); }

  //! \brief Sets the dimension of the system vector while retaining content.
  virtual void redim(size_t n)
  {
    this->utl::vector<Real>::resize(n,utl::RETAIN);
  }

  //! \brief Resizes the vector to length \a n.
  //! \details Will erase the previous content, but only if the size changed,
  //! unless \a forceClear is \e true.
  virtual void resize(size_t n, bool forceClear = false)
  {
    this->utl::vector<Real>::resize(n,forceClear);
  }

  //! \brief Access through pointer.
  virtual Real* getPtr() { return this->ptr(); }
  //! \brief Reference through pointer.
  virtual const Real* getRef() const { return this->ptr(); }

  //! \brief Initializes the vector to a given scalar value.
  virtual void init(Real value = Real(0)) { this->fill(value); }

  //! \brief Multiplication with a scalar.
  virtual void mult(Real alpha) { this->operator*=(alpha); }

  //! \brief Addition of another system vector to this one.
  virtual void add(const SystemVector& vec, Real scale)
  {
    this->utl::vector<Real>::add(static_cast<const StdVector&>(vec),scale);
  }

  //! \brief L1-norm of the vector.
  virtual Real L1norm() const { return this->asum(); }

  //! \brief L2-norm of the vector.
  virtual Real L2norm() const { return this->norm2(); }

  //! \brief Linfinity-norm of the vector.
  virtual Real Linfnorm() const { size_t off = 0; return this->normInf(off); }

  //! \brief Dumps the system vector on a specified format.
  virtual void dump(std::ostream& os, LinAlg::StorageFormat format,
                    const char* label) const { dump(*this,label,format,os); }

  //! \brief Dumps a standard vector to given output stream on specified format.
  static void dump(const utl::vector<Real>& x, const char* label,
                   LinAlg::StorageFormat format, std::ostream& os);

protected:
  //! \brief Writes the system vector to the given output stream.
  virtual std::ostream& write(std::ostream& os) const
  {
    return os << static_cast<const utl::vector<Real>&>(*this);
  }
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
  //! \brief Static method creating a matrix of the given type.
  static SystemMatrix* create(const ProcessAdm* adm, LinAlg::MatrixType mType,
                              int num_thread_SLU = 1);
  //! \brief Static method creating a matrix of the given type.
  static SystemMatrix* create(const ProcessAdm* adm, LinAlg::MatrixType mType,
                              const LinSolParams& spar);

protected:
  //! \brief Default constructor.
  SystemMatrix() : haveContributions(false) {}
  //! \brief Copy constructor.
  SystemMatrix(const SystemMatrix& a) : haveContributions(a.haveContributions){}

public:
  //! \brief Empty destructor.
  virtual ~SystemMatrix() {}

  //! \brief Returns the matrix type.
  virtual LinAlg::MatrixType getType() const = 0;

  //! \brief Creates a copy of the system matrix and returns a pointer to it.
  virtual SystemMatrix* copy() const = 0;

  //! \brief Locks or unlocks the sparsity pattern.
  virtual bool lockPattern(bool) { return false; }

  //! \brief Checks if the matrix is empty.
  virtual bool empty() const { return this->dim(0) == 0; }
  //! \brief Checks if the matrix have no non-zero contributions.
  bool isZero() const { return !haveContributions && this->dim(1) > 0; }

  //! \brief Returns the dimension of the system matrix.
  virtual size_t dim(int idim = 1) const = 0;

  //! \brief Initializes the element assembly process.
  //! \param[in] sam Auxiliary data describing the FE model topology, etc.
  //! \param[in] delayLocking If \e true, do not lock the sparsity pattern yet
  //!
  //! \details This method must be called once before the element assembly loop.
  virtual void initAssembly(const SAM& sam, bool delayLocking = false) = 0;

  //! \brief Initializes the matrix to zero assuming it is properly dimensioned.
  virtual void init() = 0;

  //! \brief Begins communication step needed in parallel matrix assembly.
  virtual bool beginAssembly() { return true; }
  //! \brief Ends communication step needed in parallel matrix assembly.
  virtual bool endAssembly() { return true; }

  //! \brief Adds an element matrix into the associated system matrix.
  //! \param[in] eM  The element matrix
  //! \param[in] sam Auxiliary data for FE assembly management
  //! \param[in] e   Identifier for the element that \a eM belongs to
  //! \return \e true on successful assembly, otherwise \e false
  virtual bool assemble(const Matrix& eM, const SAM& sam, int e) = 0;
  //! \brief Adds an element matrix into the associated system matrix.
  //! \param[in] eM  The element matrix
  //! \param[in] sam Auxiliary data for FE assembly management
  //! \param     B   The system right-hand-side vector
  //! \param[in] e   Identifier for the element that \a eM belongs to
  //! \return \e true on successful assembly, otherwise \e false
  //!
  //! \details When multi-point constraints are present, contributions from
  //! these are also added into the system right-hand-side vector, \a B.
  virtual bool assemble(const Matrix& eM, const SAM& sam,
                        SystemVector& B, int e) = 0;
  //! \brief Adds an element matrix into the associated system matrix.
  //! \param[in] eM  The element matrix
  //! \param[in] sam Auxiliary data for FE assembly management
  //! \param     B   The system right-hand-side vector
  //! \param[in] meq Matrix of element equation numbers
  //! \return \e true on successful assembly, otherwise \e false
  //!
  //! \details When multi-point constraints are present, contributions from
  //! these are also added into the system right-hand-side vector, \a B.
  virtual bool assemble(const Matrix& eM, const SAM& sam,
                        SystemVector& B, const std::vector<int>& meq) = 0;

  //! \brief Augments a similar matrix symmetrically to the current matrix.
  virtual bool augment(const SystemMatrix&, size_t, size_t) { return false; }

  //! \brief Truncates all small matrix elements to zero.
  virtual bool truncate(Real = Real(1.0e-16)) { return false; }

  //! \brief Multiplication with a scalar.
  virtual void mult(Real alpha) = 0;

  //! \brief Adds a matrix with similar structure to the current matrix.
  virtual bool add(const SystemMatrix&, Real = Real(1)) { return false; }

  //! \brief Adds a constant diagonal matrix to the current matrix.
  virtual bool add(Real) { return false; }

  //! \brief Performs a matrix-vector multiplication.
  virtual bool multiply(const SystemVector&, SystemVector&) const
  { return false; }

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param b Right-hand-side vector on input, solution vector on output
  //! \param[out] rc Reciprocal condition number of the LHS-matrix (optional)
  virtual bool solve(SystemVector& b, Real* rc = nullptr) = 0;

  //! \brief Solves the linear system of equations for a given right-hand-side.
  //! \param[in] b Right-hand-side vector
  //! \param[out] x Solution vector
  virtual bool solve(const SystemVector& b, SystemVector& x)
  {
    return this->solve(x.copy(b));
  }

  //! \brief Returns the L-infinity norm of the matrix.
  virtual Real Linfnorm() const = 0;

  //! \brief Dumps the system matrix on a specified format.
  virtual void dump(std::ostream&, LinAlg::StorageFormat,
                    const char* = nullptr) {}

  //! \brief Calculates a matrix-vector product.
  StdVector operator*(const SystemVector& b) const;

  //! \brief Solves a linear equation system.
  StdVector operator/(const SystemVector& b);

protected:
  //! \brief Writes the system matrix to the given output stream.
  virtual std::ostream& write(std::ostream& os) const { return os; }

  //! \brief Global stream operator printing the matrix contents.
  friend std::ostream& operator<<(std::ostream& os, const SystemMatrix& A)
  {
    return A.write(os);
  }

  bool haveContributions; //!< If \e true, the matrix have some non-zero terms
};

#endif
