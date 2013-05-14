// $Id$
//==============================================================================
//!
//! \file matrix.h
//!
//! \date Oct 1 2007
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Simple template classes for dense rectangular matrices and vectors.
//! \details The classes have some algebraic operators defined, such that the
//! class type \a T has to be of a numerical type, i.e., \a float or \a double.
//! The multiplication methods are implemented based on the CBLAS library if
//! either of the macro symbols USE_CBLAS or USE_MKL are defined.
//! Otherwise, inlined methods are used.
//!
//==============================================================================

#ifndef UTL_MATRIX_H
#define UTL_MATRIX_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cmath>
#ifdef USE_MKL
#include <mkl_cblas.h>
#elif defined(USE_CBLAS)
extern "C"
{
#include <cblas.h>
}
#endif

#ifdef INDEX_CHECK
#if INDEX_CHECK > 1
#define ABORT_ON_INDEX_CHECK abort()
#else
#define ABORT_ON_INDEX_CHECK
#endif
#define CHECK_INDEX(label,i,n) if (i < 1 || i > n) { \
    std::cerr << label << i <<" is out of range [1,"<< n <<"]"<< std::endl; \
    ABORT_ON_INDEX_CHECK; }
#else
#define CHECK_INDEX(label,i,n)
#define ABORT_ON_INDEX_CHECK
#endif

#ifdef SING_CHECK
#define ABORT_ON_SINGULARITY abort()
#else
#define ABORT_ON_SINGULARITY
#endif


namespace utl //! General utility classes and functions.
{
  /*!
    \brief Sub-class of std::vector with some added algebraic operations.
    \details The class type \a T has to be of a numerical type, i.e.,
    \a float or \a double.
  */

  template<class T> class vector : public std::vector<T>
  {
  public:
    //! \brief Constructor creating an empty vector.
    vector<T>() {}
    //! \brief Constructor creating a vector of length \a n.
    vector<T>(size_t n) { this->resize(n); }
    //! \brief Constructor creating a vector from an array.
    vector<T>(const T* values, size_t n) { this->fill(values,n); }

    //! \brief Overloaded assignment operator.
    vector<T>& operator=(const std::vector<T>& X)
    {
      this->std::vector<T>::operator=(X);
      return *this;
    }

    //! \brief Access through pointer.
    T* ptr() { return &this->front(); }
    //! \brief Reference through pointer.
    const T* ptr() const { return &this->front(); }

    //! \brief Index-1 based element access.
    T& operator()(size_t i)
    {
      CHECK_INDEX("vector::operator(): Index ",i,this->size());
      return this->operator[](i-1);
    }

    //! \brief Index-1 based element reference.
    const T& operator()(size_t i) const
    {
      CHECK_INDEX("vector::operator(): Index ",i,this->size());
      return this->operator[](i-1);
    }

    //! \brief Fill the vector with a scalar value.
    void fill(const T& s) { std::fill(this->begin(),this->end(),s); }
    //! \brief Fill the vector with data from an array.
    void fill(const T* values, size_t n = 0)
    {
      if (n > this->size()) this->std::vector<T>::resize(n);
      memcpy(&this->front(),values,this->size()*sizeof(T));
    }

    //! \brief Multiplication with a scalar.
    vector<T>& operator*=(const T& c);
    //! \brief Division by a scalar.
    vector<T>& operator/=(const T& d) { return this->operator*=(T(1)/d); }

    //! \brief Add the given vector \b X to \a *this.
    vector<T>& operator+=(const vector<T>& X) { return this->add(X); }
    //! \brief Subtract the given vector \b X from \a *this.
    vector<T>& operator-=(const vector<T>& X) { return this->add(X,T(-1)); }
    //! \brief Add the given vector \b X scaled by \a alfa to \a *this.
    vector<T>& add(const std::vector<T>& X, T alfa = T(1));

    //! \brief Dot product between \a *this and another vector.
    //! \param[in] v The vector to dot this vector with
    //! \param[in] off1 Offset for this vector
    //! \param[in] inc1 Increment for this vector
    //! \param[in] off2 Offset for vector \b v
    //! \param[in] inc2 Increment for vector \b v
    T dot(const std::vector<T>& v,
          size_t off1 = 0, int inc1 = 1,
          size_t off2 = 0, int inc2 = 1) const;

    //! \brief Return the Euclidean norm of the vector.
    //! \param[in] off Index offset relative to the first vector component
    //! \param[in] inc Increment in the vector component indices
    T norm2(size_t off = 0, int inc = 1) const;
    //! \brief Return the infinite norm of the vector.
    //! \param off Index offset relative to the first vector component on input,
    //! 1-based index of the largest vector component on output
    //! \param[in] inc Increment in the vector component indices
    T normInf(size_t& off, int inc = 1) const;

    //! \brief Return the sum of the absolute value of the vector elements.
    //! \param[in] off Index offset relative to the first vector component
    //! \param[in] inc Increment in the vector component indices
    T asum(size_t off = 0, int inc = 1) const;

    //! \brief Return the sum of the vector elements.
    //! \param[in] off Index offset relative to the first vector component
    //! \param[in] inc Increment in the vector component indices
    T sum(size_t off = 0, int inc = 1) const
    {
      T xsum = T(0);
      if (inc < 1) return xsum;

      const T* v = this->ptr();
      for (size_t i = off; i < this->size(); i += inc)
        xsum += v[i];
      return xsum;
    }

    //! \brief Resize the vector to length \a n.
    //! \details Will erase the previous content, but only if the size changed,
    //! unless \a forceClear is \e true.
    void resize(size_t n, bool forceClear = false)
    {
      if (forceClear)
      {
        // Erase previous content
        if (n == this->size())
          this->fill(T(0));
        else
          this->clear();
      }

      if (n == this->size()) return; // nothing to do

      if (!forceClear) this->clear();

      this->std::vector<T>::resize(n,T(0));
    }
  };


  /*!
    \brief Two-dimensional rectangular matrix with some algebraic operations.
    \details This is a 2D equivalent to the \a vector class. The matrix elements
    are stored column-wise in a one-dimensional array, such that its pointer
    might be passed to Fortran subroutines requiring 2D arrays as arguments.
  */

  template<class T> class matrix
  {
  public:
    //! \brief Constructor creating an empty matrix.
    matrix() : nrow(0), ncol(0) {}
    //! \brief Constructor creating a matrix of dimension \f$r \times c\f$.
    matrix(size_t r, size_t c) : nrow(r), ncol(c), elem(r*c) {}
    //! \brief Copy constructor, optionally creates the transpose of \a mat.
    matrix(const matrix<T>& mat, bool transposed = false) : elem(mat.size())
    {
      nrow = transposed ? mat.ncol : mat.nrow;
      ncol = transposed ? mat.nrow : mat.ncol;
      if (transposed)
        for (size_t r = 0; r < ncol; r++)
          for (size_t c = 0; c < nrow; c++)
            elem[c+nrow*r] = mat.elem[r+ncol*c];
      else
        elem.fill(mat.elem.ptr());
    }

    //! \brief Clears the matrix and sets its dimension to zero.
    void clear() { nrow = ncol = 0; elem.clear(); }

    //! \brief Resize the matrix to dimension \f$r \times c\f$.
    //! \details Will erase the previous content, but only if both
    //! the total matrix size and the number of rows in the matrix are changed.
    //! It is therefore possible to add or remove a given number of columns to
    //! the matrix without loosing the contents of the remaining columns.
    //! If \a forceClear is \e true, the old matrix content is always erased.
    void resize(size_t r, size_t c, bool forceClear = false)
    {
      if (forceClear)
      {
        // Erase previous content
        if (nrow*ncol == r*c)
          this->fill(T(0));
        else
          this->clear();
      }

      if (r == nrow && c == ncol) return; // nothing to do

      size_t oldNrow = nrow;
      size_t oldSize = this->size();
      nrow = r;
      ncol = c;
      if (this->size() == oldSize) return; // no more to do, size is unchanged

      // If the number of rows is changed the previous content must be cleared
      if (!forceClear && r != oldNrow) elem.clear();

      elem.std::vector<T>::resize(r*c,T(0));
    }

    //! \brief Query number of matrix rows.
    size_t rows() const { return nrow; }
    //! \brief Query number of matrix columns.
    size_t cols() const { return ncol; }
    //! \brief Query total matrix size.
    size_t size() const { return nrow*ncol; }
    //! \brief Check if the matrix is empty.
    bool empty() const { return elem.empty(); }

    //! \brief Access through pointer.
    T* ptr(size_t c = 0) { return &elem[nrow*c]; }
    //! \brief Reference through pointer.
    const T* ptr(size_t c = 0) const { return &elem[nrow*c]; }
    //! \brief Type casting to a one-dimensional vector.
    operator const vector<T>&() const { return elem; }

    //! \brief Index-1 based element access.
    //! \details Assuming column-wise storage as in Fortran.
    T& operator()(size_t r, size_t c)
    {
      CHECK_INDEX("matrix::operator(): Row-index ",r,nrow);
      CHECK_INDEX("matrix::operator(): Column-index ",c,ncol);
      return elem[r-1+nrow*(c-1)];
    }

    //! \brief Index-1 based element reference.
    //! \details Assuming column-wise storage as in Fortran.
    const T& operator()(size_t r, size_t c) const
    {
      CHECK_INDEX("matrix::operator(): Row-index ",r,nrow);
      CHECK_INDEX("matrix::operator(): Column-index ",c,ncol);
      return elem[r-1+nrow*(c-1)];
    }

    //! \brief Extract a row from the matrix.
    vector<T> getRow(size_t r) const
    {
      CHECK_INDEX("matrix::getRow: Row-index ",r,nrow);
      if (nrow < 2) return elem;
      vector<T> row(ncol);
      for (size_t i = 0; i < ncol; i++)
        row[i] = elem[r-1+nrow*i];
      return row;
    }

    //! \brief Extract a column from the matrix.
    vector<T> getColumn(size_t c) const
    {
      CHECK_INDEX("matrix::getColumn: Column-index ",c,ncol);
      vector<T> col(nrow);
      memcpy(col.ptr(),this->ptr(c-1),nrow*sizeof(T));
      return col;
    }

    //! \brief Fill a column of the matrix.
    void fillColumn(size_t c, const std::vector<T>& data)
    {
      CHECK_INDEX("matrix::getColumn: Column-index ",c,ncol);
      size_t ndata = nrow > data.size() ? data.size() : nrow;
      memcpy(this->ptr(c-1),&data.front(),ndata*sizeof(T));
    }

    //! \brief Fill a column of the matrix.
    void fillColumn(size_t c, const T* data)
    {
      CHECK_INDEX("matrix::fillColumn: Column-index ",c,ncol);
      memcpy(this->ptr(c-1),data,nrow*sizeof(T));
    }

    //! \brief Fill the matrix with a scalar value.
    void fill(const T& s) { std::fill(elem.begin(),elem.end(),s); }
    //! \brief Fill the matrix with data from an array.
    void fill(const T* values, size_t n = 0) { elem.fill(values,n); }

    //! \brief Create a diagonal matrix.
    matrix<T>& diag(const T& d, size_t dim = 0)
    {
      if (dim > 0)
        this->resize(dim,dim,true);
      else
        this->resize(nrow,ncol,true);
      for (size_t r = 0; r < nrow && r < ncol; r++)
        elem[r+nrow*r] = d;
      return *this;
    }

    //! \brief Replace the current matrix by its transpose.
    matrix<T>& transpose()
    {
      matrix<T> tmp(*this);
      for (size_t r = 0; r < nrow; r++)
        for (size_t c = 0; c < ncol; c++)
          elem[c+ncol*r] = tmp.elem[r+nrow*c];

      nrow = tmp.ncol;
      ncol = tmp.nrow;
      return *this;
    }

#define THIS(i,j) this->operator()(i,j)

    //! \brief Compute the determinant of a square matrix.
    T det() const
    {
      if (ncol == 1 && nrow >= 1)
        return THIS(1,1);
      else if (ncol == 2 && nrow >= 2)
        return THIS(1,1)*THIS(2,2) - THIS(2,1)*THIS(1,2);
      else if (ncol == 3 && nrow >= 3)
        return THIS(1,1)*(THIS(2,2)*THIS(3,3) - THIS(3,2)*THIS(2,3))
          -    THIS(1,2)*(THIS(2,1)*THIS(3,3) - THIS(3,1)*THIS(2,3))
          +    THIS(1,3)*(THIS(2,1)*THIS(3,2) - THIS(3,1)*THIS(2,2));
      else if (ncol > 0 && nrow > 0) {
        std::cerr <<"matrix::det: Not available for "
                  << nrow <<"x"<< ncol <<" matrices"<< std::endl;
        ABORT_ON_SINGULARITY;
        return T(-999);
      }
      else
        return T(0);
    }

    //! \brief Compute the inverse of a square matrix.
    //! \param[in] tol Division by zero tolerance
    //! \return Determinant of the matrix
    T inverse(const T& tol = T(0))
    {
      T det = this->det();
      if (det == T(-999))
        return det;
      else if (det <= tol && det >= -tol) {
        std::cerr <<"matrix::inverse: Singular matrix |A|="<< det << std::endl;
        ABORT_ON_SINGULARITY;
        return T(0);
      }

      if (ncol == 1)
        THIS(1,1) = T(1) / det;
      else if (ncol == 2) {
        matrix<T> B(2,2);
        B(1,1) =  THIS(2,2) / det;
        B(2,1) = -THIS(2,1) / det;
        B(1,2) = -THIS(1,2) / det;
        B(2,2) =  THIS(1,1) / det;
        *this = B;
      }
      else if (ncol == 3) {
        matrix<T> B(3,3);
        B(1,1) =  (THIS(2,2)*THIS(3,3) - THIS(3,2)*THIS(2,3)) / det;
        B(2,1) = -(THIS(2,1)*THIS(3,3) - THIS(3,1)*THIS(2,3)) / det;
        B(3,1) =  (THIS(2,1)*THIS(3,2) - THIS(3,1)*THIS(2,2)) / det;
        B(1,2) = -(THIS(1,2)*THIS(3,3) - THIS(3,2)*THIS(1,3)) / det;
        B(2,2) =  (THIS(1,1)*THIS(3,3) - THIS(3,1)*THIS(1,3)) / det;
        B(3,2) = -(THIS(1,1)*THIS(3,2) - THIS(3,1)*THIS(1,2)) / det;
        B(1,3) =  (THIS(1,2)*THIS(2,3) - THIS(2,2)*THIS(1,3)) / det;
        B(2,3) = -(THIS(1,1)*THIS(2,3) - THIS(2,1)*THIS(1,3)) / det;
        B(3,3) =  (THIS(1,1)*THIS(2,2) - THIS(2,1)*THIS(1,2)) / det;
        *this = B;
      }

      return det;
    }

    //! \brief Check for symmetry.
    //! \param[in] tol Comparison tolerance
    bool isSymmetric(const T& tol = T(0)) const
    {
      if (nrow != ncol) return false;

      for (size_t r = 0; r < nrow; r++)
        for (size_t c = 0; c < r; c++)
        {
          T diff = elem[r+nrow*c] - elem[c+nrow*r];
          if (diff < -tol || diff > tol) return false;
        }

      return true;
    }

    //! \brief Add the given matrix \b A to \a *this.
    matrix<T>& operator+=(const matrix<T>& A) { return this->add(A); }
    //! \brief Subtract the given matrix \b A from \a *this.
    matrix<T>& operator-=(const matrix<T>& A) { return this->add(A,T(-1)); }
    //! \brief Add the given matrix \b A scaled by \a alfa to \a *this.
    matrix<T>& add(const matrix<T>& A, T alfa = T(1));

    //! \brief Multiplication with a scalar.
    matrix<T>& operator*=(const T& c) { return this->multiply(c); }
    //! \brief Division by a scalar.
    matrix<T>& operator/=(const T& d) { return this->multiply(T(1)/d); }
    //! \brief Multiplication of this matrix by a scalar \a c.
    matrix<T>& multiply(const T& c);

    /*! \brief Matrix-matrix multiplication.
      \details Performs one of the following operations (\b C = \a *this):
      -# \f$ {\bf C} = {\bf A} {\bf B} \f$
      -# \f$ {\bf C} = {\bf A}^T {\bf B} \f$
      -# \f$ {\bf C} = {\bf A} {\bf B}^T \f$
      -# \f$ {\bf C} = {\bf A}^T {\bf B}^T \f$
      -# \f$ {\bf C} = {\bf C} + {\bf A} {\bf B} \f$
      -# \f$ {\bf C} = {\bf C} + {\bf A}^T {\bf B} \f$
      -# \f$ {\bf C} = {\bf C} + {\bf A} {\bf B}^T \f$
      -# \f$ {\bf C} = {\bf C} + {\bf A}^T {\bf B}^T \f$
    */
    matrix<T>& multiply(const matrix<T>& A, const matrix<T>& B,
                        bool transA = false, bool transB = false,
                        bool addTo = false);

    /*! \brief Matrix-matrix multiplication.
      \details Performs one of the following operations (\b C = \a *this):
      -# \f$ {\bf C} = {\bf A} {\bf B} \f$
      -# \f$ {\bf C} = {\bf A}^T {\bf B} \f$
      -# \f$ {\bf C} = {\bf C} + {\bf A} {\bf B} \f$
      -# \f$ {\bf C} = {\bf C} + {\bf A}^T {\bf B} \f$

      The matrix \b B is here represented by a one-dimensional vector,
      and its number of rows is assumed to match the number of columns in
      \b A (or its transpose) and its number of columns is then the total
      vector length divided by the number of rows.
    */
    bool multiplyMat(const matrix<T>& A, const std::vector<T>& B,
                     bool transA = false, bool addTo = false);

    /*! \brief Matrix-matrix multiplication.
      \details Performs one of the following operations (\b C = \a *this):
      -# \f$ {\bf C} = {\bf A} {\bf B} \f$
      -# \f$ {\bf C} = {\bf A} {\bf B}^T \f$
      -# \f$ {\bf C} = {\bf C} + {\bf A} {\bf B} \f$
      -# \f$ {\bf C} = {\bf C} + {\bf A} {\bf B}^T \f$

      The matrix \b A is here represented by a one-dimensional vector,
      and its number of columns is assumed to match the number of rows in
      \b B (or its transpose) and its number of rows is then the total
      vector length divided by the number of columns.
    */
    bool multiplyMat(const std::vector<T>& A, const matrix<T>& B,
                     bool transB = false, bool addTo = false);

    /*! \brief Matrix-vector multiplication.
      \details Performs one of the following operations (\b A = \a *this):
      -# \f$ {\bf Y} = {\bf A} {\bf X} \f$
      -# \f$ {\bf Y} = {\bf A}^T {\bf X} \f$
      -# \f$ {\bf Y} = {\bf Y} + {\bf A} {\bf X} \f$
      -# \f$ {\bf Y} = {\bf Y} + {\bf A}^T {\bf X} \f$
    */
    bool multiply(const std::vector<T>& X, std::vector<T>& Y,
                  bool transA = false, bool addTo = false) const;

    //! \brief Outer product between two vectors.
    bool outer_product(const std::vector<T>& X, const std::vector<T>& Y,
                       bool addTo = false)
    {
      if (!addTo)
        this->resize(X.size(),Y.size());
      else if (X.size() != nrow || Y.size() != ncol)
      {
        std::cerr <<"matrix::outer_product: Incompatible matrix and vectors: A("
                  << nrow <<','<< ncol <<"), X("
                  << X.size() <<"), Y("<< Y.size() <<")\n"
                  <<"                       when computing A += X*Y^t"
                  << std::endl;
        ABORT_ON_INDEX_CHECK;
        return false;
      }

      if (addTo)
        for (size_t j = 0; j < ncol; j++)
          for (size_t i = 0; i < nrow; i++)
            elem[i+nrow*j] += X[i]*Y[j];
      else
        for (size_t j = 0; j < ncol; j++)
          for (size_t i = 0; i < nrow; i++)
            elem[i+nrow*j] = X[i]*Y[j];

      return true;
    }

    //! \brief Return the infinite norm of the matrix.
    T normInf() const
    {
      if (nrow == 0) return T(0);

      // Compute row sums
      vector<T> sums(nrow);
      for (size_t i = 0; i < nrow; i++)
        sums[i] = elem.asum(i,ncol);
      return *std::max_element(sums.begin(),sums.end());
    }

  private:
    //! \brief Check dimension compatibility for matrix-vector multiplication.
    bool compatible(const std::vector<T>& X, bool transA) const
    {
      if (nrow > 0 && ncol > 0)
        if ((transA ? nrow : ncol) == X.size()) return true;

      std::cerr <<"matrix::multiply: Incompatible matrices: A("
                << nrow <<','<< ncol <<"), X("<< X.size() <<")\n"
                <<"                  when computing Y = "
                << (transA ? "A^t":"A") <<" * X"<< std::endl;
      ABORT_ON_INDEX_CHECK;
      return false;
    }

    //! \brief Check dimension compatibility for matrix-matrix multiplication.
    bool compatible(const matrix<T>& A, const matrix<T>& B,
                    bool transA, bool transB, size_t& M, size_t& N, size_t& K)
    {
      M = transA ? A.ncol : A.nrow;
      N = transB ? B.nrow : B.ncol;
      K = transA ? A.nrow : A.ncol;
      if (K == (transB ? B.ncol : B.nrow)) return true;

      std::cerr <<"matrix::multiply: Incompatible matrices: A("
                << A.nrow <<','<< A.ncol <<"), B("
                << B.nrow <<','<< B.ncol <<")\n"
                <<"                  when computing C = "
                << (transA ? "A^t":"A") <<" * "
                << (transB ? "B^t":"B") << std::endl;
      ABORT_ON_INDEX_CHECK;
      return false;
    }

    //! \brief Check dimension compatibility for matrix-matrix multiplication,
    //! when the matrix B is represented by a one-dimensional vector.
    bool compatible(const matrix<T>& A, const std::vector<T>& B,
                    bool transA, size_t& M, size_t& N, size_t& K)
    {
      M = transA ? A.ncol : A.nrow;
      K = transA ? A.nrow : A.ncol;
      N = K > 0 ? B.size()/K : 0;
      if (N*K == B.size() && !B.empty()) return true;

      std::cerr <<"matrix::multiply: Incompatible matrices: A("
                << A.nrow <<','<< A.ncol <<"), B(r*c="<< B.size() <<")\n"
                <<"                  when computing C = "
                << (transA ? "A^t":"A") <<" * B"<< std::endl;
      ABORT_ON_INDEX_CHECK;
      return false;
    }

    //! \brief Check dimension compatibility for matrix-matrix multiplication,
    //! when the matrix A is represented by a one-dimensional vector.
    bool compatible(const std::vector<T>& A, const matrix<T>& B,
                    bool transB, size_t& M, size_t& N, size_t& K)
    {
      N = transB ? B.nrow : B.ncol;
      K = transB ? B.ncol : B.nrow;
      M = K > 0 ? A.size() / K : 0;
      if (M*K == A.size() && !A.empty()) return true;

      std::cerr <<"matrix::multiply: Incompatible matrices: A(r*c="<< A.size()
                <<"), B("<< B.nrow <<","<< B.ncol <<")\n"
                <<"                  when computing C = A * "
                << (transB ? "B^t":"B") << std::endl;
      ABORT_ON_INDEX_CHECK;
      return false;
    }

  private:
    size_t    nrow; //!< Number of matrix rows
    size_t    ncol; //!< Number of matrix columns
    vector<T> elem; //!< Actual matrix elements, stored column by column
  };


  /*!
    \brief Three-dimensional rectangular matrix with some algebraic operations.
    \details This is a 3D equivalent to the \a matrix class. The matrix elements
    are stored column-wise in a one-dimensional array, such that its pointer
    might be passed to Fortran subroutines requiring 3D arrays as arguments.
  */

  template<class T> class matrix3d
  {
  public:
    //! \brief Constructor creating an empty matrix.
    matrix3d() { n[0] = n[1] = n[2] = 0; }
    //! \brief Constructor creating a matrix of dimension
    //! \f$n_1 \times n_2 \times n_3\f$.
    matrix3d(size_t n_1, size_t n_2, size_t n_3) : elem(n_1*n_2*n_3)
    {
      n[0] = n_1;
      n[1] = n_2;
      n[2] = n_3;
    }

    //! \brief Resize the matrix to dimension \f$n_1 \times n_2 \times n_3\f$.
    //! \details Will erase the previous content, but only if both
    //! the total matrix size, and n_1 or n_2 in the matrix are changed.
    //! It is therefore possible to add or remove a given number of elements in
    //! the third dimension of the matrix without loosing the contents of the
    //! first and second dimensions.
    //! If \a forceClear is \e true, the old matrix content is always erased.
    void resize(size_t n_1, size_t n_2, size_t n_3, bool forceClear = false)
    {
      if (forceClear)
      {
        // Erase previous content
        if (n[0]*n[1]*n[2] == n_1*n_2*n_3)
          this->fill(T(0));
        else
        {
          n[0] = n[1] = n[2] = 0;
          elem.clear();
        }
      }

      if (n[0] == n_1 && n[1] == n_2 && n[2] == n_3) return; // nothing to do

      size_t oldn1 = n[0];
      size_t oldn2 = n[1];
      size_t oldSize = this->size();
      n[0] = n_1;
      n[1] = n_2;
      n[2] = n_3;
      if (this->size() == oldSize) return; // no more to do, size is unchanged

      // If the size in any of the first two dimensions are changed
      // the previous content must be cleared
      if (!forceClear && (n[0] != oldn1 || n[1] != oldn2)) elem.clear();

      elem.std::vector<T>::resize(n[0]*n[1]*n[2],T(0));
    }

    //! \brief Query dimensions.
    size_t dim(short int d = 1) const { return d > 0 && d <= 3 ? n[d-1] : 0; }
    //! \brief Query total matrix size.
    size_t size() const { return n[0]*n[1]*n[2]; }
    //! \brief Check if the matrix is empty.
    bool empty() const { return elem.empty(); }

    //! \brief Type casting to a one-dimensional vector.
    operator const vector<T>&() const { return elem; }

    //! \brief Index-1 based element access.
    //! \details Assuming column-wise storage as in Fortran.
    T& operator()(size_t i1, size_t i2, size_t i3)
    {
      CHECK_INDEX("matrix3d::operator(): First index " ,i1,n[0]);
      CHECK_INDEX("matrix3d::operator(): Second index ",i2,n[1]);
      CHECK_INDEX("matrix3d::operator(): Third index " ,i3,n[2]);
      return elem[i1-1 + n[0]*((i2-1) + n[1]*(i3-1))];
    }

    //! \brief Index-1 based element access.
    //! \details Assuming column-wise storage as in Fortran.
    const T& operator()(size_t i1, size_t i2, size_t i3) const
    {
      CHECK_INDEX("matrix3d::operator(): First index " ,i1,n[0]);
      CHECK_INDEX("matrix3d::operator(): Second index ",i2,n[1]);
      CHECK_INDEX("matrix3d::operator(): Third index " ,i3,n[2]);
      return elem[i1-1 + n[0]*((i2-1) + n[1]*(i3-1))];
    }

    //! \brief Access through pointer.
    T* ptr() { return &elem.front(); }
    //! \brief Reference through pointer.
    const T* ptr() const { return &elem.front(); }

    //! \brief Fill the matrix with a scalar value.
    void fill(const T& s) { std::fill(elem.begin(),elem.end(),s); }

    /*! \brief Matrix-matrix multiplication.
      \details Performs one of the following operations (\b C = \a *this):
      -# \f$ {\bf C} = {\bf A} {\bf B} \f$
      -# \f$ {\bf C} = {\bf A}^T {\bf B} \f$
      -# \f$ {\bf C} = {\bf C} + {\bf A} {\bf B} \f$
      -# \f$ {\bf C} = {\bf C} + {\bf A}^T {\bf B} \f$
    */
    bool multiply(const matrix<T>& A, const matrix3d<T>& B,
                  bool transA = false, bool addTo = false);

  private:
    //! \brief Check dimension compatibility for matrix-matrix multiplication.
    bool compatible(const matrix<T>& A, const matrix3d<T>& B,
                    bool transA, size_t& M, size_t& N, size_t& K)
    {
      M = transA ? A.cols() : A.rows();
      N = B.n[1]*B.n[2];
      K = transA ? A.rows() : A.cols();
      if (K == B.n[0]) return true;

      std::cerr <<"matrix3d::multiply: Incompatible matrices: A("
                << A.rows() <<','<< A.cols() <<"), B("
                << B.n[0] <<','<< B.n[1] <<','<< B.n[2] <<")\n"
                <<"                  when computing C = "
                << (transA ? "A^t":"A") <<" * B"<< std::endl;
      ABORT_ON_INDEX_CHECK;
      return false;
    }

  private:
    size_t    n[3]; //!< Dimension of the 3D matrix
    vector<T> elem; //!< Actual matrix elements, stored column by column
  };


#if defined(USE_CBLAS) || defined(USE_MKL)
  //============================================================================
  //===   BLAS-implementation of the matrix/vector multiplication methods   ====
  //============================================================================

  template<> inline
  vector<float>& vector<float>::operator*=(const float& c)
  {
    cblas_sscal(this->size(),c,this->ptr(),1);
    return *this;
  }

  template<> inline
  vector<double>& vector<double>::operator*=(const double& c)
  {
    cblas_dscal(this->size(),c,this->ptr(),1);
    return *this;
  }

  template<> inline
  float vector<float>::dot(const std::vector<float>& v,
                           size_t o1, int i1, size_t o2, int i2) const
  {
    int n1 = i1 > 1 || i1 < -1 ? this->size()/abs(i1) : this->size()-o1;
    int n2 = i2 > 1 || i2 < -1 ? v.size()/abs(i2) : v.size()-o2;
    int n  = n1 < n2 ? n1 : n2;
    return cblas_sdot(n,this->ptr()+o1,i1,&v.front()+o2,i2);
  }

  template<> inline
  double vector<double>::dot(const std::vector<double>& v,
                             size_t o1, int i1, size_t o2, int i2) const
  {
    int n1 = i1 > 1 || i1 < -1 ? this->size()/abs(i1) : this->size()-o1;
    int n2 = i2 > 1 || i2 < -1 ? v.size()/abs(i2) : v.size()-o2;
    int n  = n1 < n2 ? n1 : n2;
    return cblas_ddot(n,this->ptr()+o1,i1,&v.front()+o2,i2);
  }

  template<> inline
  float vector<float>::norm2(size_t off, int inc) const
  {
    int n = inc > 1 || inc < -1 ? this->size()/abs(inc) : this->size()-off;
    return cblas_snrm2(n,this->ptr()+off,inc);
  }

  template<> inline
  double vector<double>::norm2(size_t off, int inc) const
  {
    int n = inc > 1 || inc < -1 ? this->size()/abs(inc) : this->size()-off;
    return cblas_dnrm2(n,this->ptr()+off,inc);
  }

  template<> inline
  float vector<float>::normInf(size_t& off, int inc) const
  {
    if (inc < 1) return 0.0f;

    const float* v = this->ptr() + off;
    off = 1 + cblas_isamax(this->size()/inc,v,inc);
    return fabsf(v[(off-1)*inc]);
  }

  template<> inline
  double vector<double>::normInf(size_t& off, int inc) const
  {
    if (inc < 1) return 0.0;

    const double* v = this->ptr() + off;
    off = 1 + cblas_idamax(this->size()/inc,v,inc);
    return fabs(v[(off-1)*inc]);
  }

  template<> inline
  float vector<float>::asum(size_t off, int inc) const
  {
    int n = inc > 1 || inc < -1 ? this->size()/abs(inc) : this->size()-off;
    return cblas_sasum(n,this->ptr()+off,inc);
  }

  template<> inline
  double vector<double>::asum(size_t off, int inc) const
  {
    int n = inc > 1 || inc < -1 ? this->size()/abs(inc) : this->size()-off;
    return cblas_dasum(n,this->ptr()+off,inc);
  }

  template<> inline
  vector<float>& vector<float>::add(const std::vector<float>& X, float alfa)
  {
    size_t n = this->size() < X.size() ? this->size() : X.size();
    cblas_saxpy(n,alfa,&X.front(),1,this->ptr(),1);
    return *this;
  }

  template<> inline
  vector<double>& vector<double>::add(const std::vector<double>& X, double alfa)
  {
    size_t n = this->size() < X.size() ? this->size() : X.size();
    cblas_daxpy(n,alfa,&X.front(),1,this->ptr(),1);
    return *this;
  }

  template<> inline
  matrix<float>& matrix<float>::add(const matrix<float>& X, float alfa)
  {
    size_t n = this->size() < X.size() ? this->size() : X.size();
    cblas_saxpy(n,alfa,X.ptr(),1,this->ptr(),1);
    return *this;
  }

  template<> inline
  matrix<double>& matrix<double>::add(const matrix<double>& X, double alfa)
  {
    size_t n = this->size() < X.size() ? this->size() : X.size();
    cblas_daxpy(n,alfa,X.ptr(),1,this->ptr(),1);
    return *this;
  }

  template<> inline
  matrix<float>& matrix<float>::multiply(const float& c)
  {
    cblas_sscal(this->size(),c,this->ptr(),1);
    return *this;
  }

  template<> inline
  matrix<double>& matrix<double>::multiply(const double& c)
  {
    cblas_dscal(this->size(),c,this->ptr(),1);
    return *this;
  }

  template<> inline
  bool matrix<float>::multiply(const std::vector<float>& X,
                               std::vector<float>& Y,
                               bool transA, bool addTo) const
  {
    if (!this->compatible(X,transA)) return false;
    if (!addTo) Y.resize(transA ? ncol : nrow);

    cblas_sgemv(CblasColMajor,
                transA ? CblasTrans : CblasNoTrans,
                nrow, ncol, 1.0f,
                this->ptr(), nrow,
                &X.front(), 1, addTo ? 1.0f : 0.0f,
                &Y.front(), 1);

    return true;
  }

  template<> inline
  bool matrix<double>::multiply(const std::vector<double>& X,
                                std::vector<double>& Y,
                                bool transA, bool addTo) const
  {
    if (!this->compatible(X,transA)) return false;
    if (!addTo) Y.resize(transA ? ncol : nrow);

    cblas_dgemv(CblasColMajor,
                transA ? CblasTrans : CblasNoTrans,
                nrow, ncol, 1.0,
                this->ptr(), nrow,
                &X.front(), 1, addTo ? 1.0 : 0.0,
                &Y.front(), 1);

    return true;
  }

  template<> inline
  matrix<float>& matrix<float>::multiply(const matrix<float>& A,
                                         const matrix<float>& B,
                                         bool transA, bool transB, bool addTo)
  {
    size_t M, N, K;
    if (!this->compatible(A,B,transA,transB,M,N,K)) return *this;
    if (!addTo) this->resize(M,N);

    cblas_sgemm(CblasColMajor,
                transA ? CblasTrans : CblasNoTrans,
                transB ? CblasTrans : CblasNoTrans,
                M, N, K, 1.0f,
                A.ptr(), A.nrow,
                B.ptr(), B.nrow,
                addTo ? 1.0f : 0.0f,
                this->ptr(), nrow);

    return *this;
  }

  template<> inline
  matrix<double>& matrix<double>::multiply(const matrix<double>& A,
                                           const matrix<double>& B,
                                           bool transA, bool transB, bool addTo)
  {
    size_t M, N, K;
    if (!this->compatible(A,B,transA,transB,M,N,K)) return *this;
    if (!addTo) this->resize(M,N);

    cblas_dgemm(CblasColMajor,
                transA ? CblasTrans : CblasNoTrans,
                transB ? CblasTrans : CblasNoTrans,
                M, N, K, 1.0,
                A.ptr(), A.nrow,
                B.ptr(), B.nrow,
                addTo ? 1.0 : 0.0,
                this->ptr(), nrow);

    return *this;
  }

  template<> inline
  bool matrix<float>::multiplyMat(const matrix<float>& A,
                                  const std::vector<float>& B,
                                  bool transA, bool addTo)
  {
    size_t M, N, K;
    if (!this->compatible(A,B,transA,M,N,K)) return false;
    if (!addTo) this->resize(M,N);

    cblas_sgemm(CblasColMajor,
                transA ? CblasTrans : CblasNoTrans, CblasNoTrans,
                M, N, K, 1.0f,
                A.ptr(), A.nrow,
                &B.front(), K,
                addTo ? 1.0f : 0.0f,
                this->ptr(), nrow);

    return true;
  }

  template<> inline
  bool matrix<double>::multiplyMat(const matrix<double>& A,
                                   const std::vector<double>& B,
                                   bool transA, bool addTo)
  {
    size_t M, N, K;
    if (!this->compatible(A,B,transA,M,N,K)) return false;
    if (!addTo) this->resize(M,N);

    cblas_dgemm(CblasColMajor,
                transA ? CblasTrans : CblasNoTrans, CblasNoTrans,
                M, N, K, 1.0,
                A.ptr(), A.nrow,
                &B.front(), K,
                addTo ? 1.0 : 0.0,
                this->ptr(), nrow);

    return true;
  }

  template<> inline
  bool matrix<float>::multiplyMat(const std::vector<float>& A,
                                  const matrix<float>& B,
                                  bool transB, bool addTo)
  {
    size_t M, N, K;
    if (!this->compatible(A,B,transB,M,N,K)) return false;
    if (!addTo) this->resize(M,N);

    cblas_sgemm(CblasColMajor,
                CblasNoTrans, transB ? CblasTrans : CblasNoTrans,
                M, N, K, 1.0f,
                &A.front(), M,
                B.ptr(), B.nrow,
                addTo ? 1.0f : 0.0f,
                this->ptr(), nrow);

    return true;
  }

  template<> inline
  bool matrix<double>::multiplyMat(const std::vector<double>& A,
                                   const matrix<double>& B,
                                   bool transB, bool addTo)
  {
    size_t M, N, K;
    if (!this->compatible(A,B,transB,M,N,K)) return false;
    if (!addTo) this->resize(M,N);

    cblas_dgemm(CblasColMajor,
                CblasNoTrans, transB ? CblasTrans : CblasNoTrans,
                M, N, K, 1.0,
                &A.front(), M,
                B.ptr(), B.nrow,
                addTo ? 1.0 : 0.0,
                this->ptr(), nrow);

    return true;
  }

  template<> inline
  bool matrix3d<float>::multiply(const matrix<float>& A,
                                 const matrix3d<float>& B,
                                 bool transA, bool addTo)
  {
    size_t M, N, K;
    if (!this->compatible(A,B,transA,M,N,K)) return false;
    if (!addTo) this->resize(M,B.n[1],B.n[2]);

    cblas_sgemm(CblasColMajor,
                transA ? CblasTrans : CblasNoTrans, CblasNoTrans,
                M, N, K, 1.0f,
                A.ptr(), A.rows(),
                B.ptr(), B.n[0],
                addTo ? 1.0f : 0.0f,
                this->ptr(), n[0]);

    return true;
  }

  template<> inline
  bool matrix3d<double>::multiply(const matrix<double>& A,
                                  const matrix3d<double>& B,
                                  bool transA, bool addTo)
  {
    size_t M, N, K;
    if (!this->compatible(A,B,transA,M,N,K)) return false;
    if (!addTo) this->resize(M,B.n[1],B.n[2]);

    cblas_dgemm(CblasColMajor,
                transA ? CblasTrans : CblasNoTrans, CblasNoTrans,
                M, N, K, 1.0,
                A.ptr(), A.rows(),
                B.ptr(), B.n[0],
                addTo ? 1.0 : 0.0,
                this->ptr(), n[0]);

    return true;
  }

#else
  //============================================================================
  //===   Non-BLAS inlined implementations (slow...)   =========================
  //============================================================================

  template<class T> inline
  vector<T>& vector<T>::operator*=(const T& c)
  {
    for (size_t i = 0; i < this->size(); i++)
      std::vector<T>::operator[](i) *= c;
    return *this;
  }

  template<class T> inline
  T vector<T>::dot(const std::vector<T>& v,
                   size_t o1, int i1, size_t o2, int i2) const
  {
    size_t i, j;
    T dotprod = T(0);
    for (i = o1, j = o2; i < this->size() && j < v.size(); i += i1, j += i2)
      dotprod += this->operator[](i) * v[j];
    return dotprod;
  }

  template<class T> inline
  T vector<T>::norm2(size_t off, int inc) const
  {
    double xsum = 0.0;
    if (inc < 1) return xsum;

    // Warning: This might overflow or underflow for large/small values
    const T* v = this->ptr();
    for (size_t i = off; i < this->size(); i += inc)
      xsum += v[i]*v[i];
    return sqrt(xsum);
  }

  template<class T> inline
  T vector<T>::normInf(size_t& off, int inc) const
  {
    T xmax = T(0);
    if (inc < 1) return xmax;

    const T* v = this->ptr();
    for (size_t i = off; i < this->size(); i += inc)
      if (v[i] > xmax)
      {
        off = 1+i/inc;
        xmax = v[i];
      }
      else if (v[i] < -xmax)
      {
        off = 1+i/inc;
        xmax = -v[i];
      }

    return xmax;
  }

  template<class T> inline
  T vector<T>::asum(size_t off, int inc) const
  {
    T xsum = T(0);
    if (inc < 1) return xsum;

    const T* v = this->ptr();
    for (size_t i = off; i < this->size(); i += inc)
      xsum += v[i] < T(0) ? -v[i] : v[i];
    return xsum;
  }

  template<class T> inline
  vector<T>& vector<T>::add(const std::vector<T>& X, T alfa)
  {
    T* p = this->ptr();
    const T* q = &X.front();
    for (size_t i = 0; i < this->size() && i < X.size(); i++, p++, q++)
      *p += alfa*(*q);
    return *this;
  }

  template<class T> inline
  matrix<T>& matrix<T>::add(const matrix<T>& X, T alfa)
  {
    T* p = this->ptr();
    const T* q = X.ptr();
    for (size_t i = 0; i < this->size() && i < X.size(); i++, p++, q++)
      *p += alfa*(*q);
    return *this;
  }

  template<class T> inline
  matrix<T>& matrix<T>::multiply(const T& c)
  {
    for (size_t i = 0; i < elem.size(); i++)
      elem[i] *= c;
    return *this;
  }

  template<class T> inline
  bool matrix<T>::multiply(const std::vector<T>& X, std::vector<T>& Y,
                           bool transA, bool addTo) const
  {
    if (!this->compatible(X,Y,transA)) return false;
    if (!addTo)
    {
      Y.clear();
      Y.resize(transA ? ncol : nrow, T(0));
    }

    for (size_t i = 0; i < Y.size(); i++)
      for (size_t j = 0; j < X.size(); j++)
        if (transA)
          Y[i] += THIS(j+1,i+1) * X[j];
        else
          Y[i] += THIS(i+1,j+1) * X[j];

    return true;
  }

  template<class T> inline
  matrix<T>& matrix<T>::multiply(const matrix<T>& A,
                                 const matrix<T>& B,
                                 bool transA, bool transB, bool addTo)
  {
    size_t M, N, K;
    if (!this->compatible(A,B,transA,transB,M,N,K)) return *this;
    if (!addTo) this->resize(M,N,true);

    for (size_t i = 1; i <= M; i++)
      for (size_t j = 1; j <= N; j++)
        for (size_t k = 1; k <= K; k++)
          if (transA && transB)
            THIS(i,j) += A(k,i)*B(j,k);
          else if (transA)
            THIS(i,j) += A(k,i)*B(k,j);
          else if (transB)
            THIS(i,j) += A(i,k)*B(j,k);
          else
            THIS(i,j) += A(i,k)*B(k,j);

    return *this;
  }

  template<class T> inline
  bool matrix<T>::multiplyMat(const matrix<T>& A, const std::vector<T>& B,
                              bool transA, bool addTo)
  {
    size_t M, N, K;
    if (!this->compatible(A,B,transA,M,N,K)) return false;
    if (!addTo) this->resize(M,N,true);

    for (size_t i = 1; i <= M; i++)
      for (size_t j = 1; j <= N; j++)
        for (size_t k = 1; k <= K; k++)
          if (transA)
            THIS(i,j) += A(k,i)*B[k-1+K*(j-1)];
          else
            THIS(i,j) += A(i,k)*B[k-1+K*(j-1)];

    return true;
  }

  template<class T> inline
  bool matrix<T>::multiplyMat(const std::vector<T>& A, const matrix<T>& B,
                              bool transB, bool addTo)
  {
    size_t M, N, K;
    if (!this->compatible(A,B,transB,M,N,K)) return false;
    if (!addTo) this->resize(M,N,true);

    for (size_t i = 1; i <= M; i++)
      for (size_t j = 1; j <= N; j++)
        for (size_t k = 1; k <= K; k++)
          if (transB)
            THIS(i,j) += A[i-1+M*(k-1)]*B(j,k);
          else
            THIS(i,j) += A[i-1+M*(k-1)]*B(k,j);

    return true;
  }

  template<class T> inline
  bool matrix3d<T>::multiply(const matrix<T>& A, const matrix3d<T>& B,
                             bool transA, bool addTo)
  {
    size_t M, N, K;
    if (!this->compatible(A,B,transA,M,N,K)) return false;
    if (!addTo) this->resize(M,B.n[1],B.n[2],true);

    for (size_t i = 1; i <= n[0]; i++)
      for (size_t j = 1; j <= n[1]; j++)
        for (size_t k = 1; k <= K; k++)
          for (size_t l = 1; l <= n[2]; l++)
            if (transA)
              this->operator()(i,j,l) += A(k,i)*B(k,j,l);
            else
              this->operator()(i,j,l) += A(i,k)*B(k,j,l);

    return true;
  }
#endif

  //============================================================================
  //===   Global operators   ===================================================
  //============================================================================

  //! \brief Multiplication of a vector and a scalar.
  //! \return \f$ {\bf Y} = c {\bf X} \f$
  template<class T> vector<T> operator*(const vector<T>& X, const T& c)
  {
    vector<T> Y(X.size());
    for (size_t i = 0; i < X.size(); i++)
      Y[i] = X[i]*c;
    return Y;
  }

  //! \brief Division of a vector by a scalar.
  //! \return \f$ {\bf Y} = \frac{1}{d} {\bf X} \f$
  template<class T> vector<T> operator/(const vector<T>& X, const T& d)
  {
    vector<T> Y(X.size());
    T c = T(1) / d;
    for (size_t i = 0; i < X.size(); i++)
      Y[i] = X[i]*c;
    return Y;
  }

  //! \brief Addition of two vectors.
  //! \return \f$ {\bf Z} = {\bf X} + {\bf Y} \f$
  template<class T> vector<T> operator+(const vector<T>& X, const vector<T>& Y)
  {
    vector<T> Z(X.size());
    for (size_t i = 0; i < X.size() && i < Y.size(); i++)
      Z[i] = X[i] + Y[i];
    return Z;
  }

  //! \brief Subtraction of two vectors.
  //! \return \f$ {\bf Z} = {\bf X} - {\bf Y} \f$
  template<class T> vector<T> operator-(const vector<T>& X, const vector<T>& Y)
  {
    vector<T> Z(X.size());
    for (size_t i = 0; i < X.size() && i < Y.size(); i++)
      Z[i] = X[i] - Y[i];
    return Z;
  }

  //! \brief Multiplication of a matrix and a scalar.
  //! \return \f$ {\bf B} = c {\bf A} \f$
  template<class T> matrix<T> operator*(const matrix<T>& A, const T& c)
  {
    matrix<T> B(A);
    return B.multiply(c);
  }

  //! \brief Multiplication of a matrix and a vector.
  //! \return \f$ {\bf Y} = {\bf A} {\bf X} \f$
  template<class T> vector<T> operator*(const matrix<T>& A, const vector<T>& X)
  {
    vector<T> Y;
    A.multiply(X,Y);
    return Y;
  }

  //! \brief Multiplication of two matrices.
  //! \return \f$ {\bf C} = {\bf A} {\bf B} \f$
  template<class T> matrix<T> operator*(const matrix<T>& A, const matrix<T>& B)
  {
    matrix<T> C;
    return C.multiply(A,B);
  }

  extern double zero_print_tol; //!< Zero tolerance for printing numbers
  extern int    nval_per_line;  //!< Number of values to print per line

  //! \brief Truncate a value to zero when it is less than a given threshold.
  //! \details Used when printing matrices for easy comparison with other
  //! matrices when they contain terms that are numerically zero, except for
  //! some round-off noice. The value of the global variable \a zero_print_tol
  //! is used as a tolerance in this method.
  template<class T> inline T trunc(const T& v)
  {
    return v > T(zero_print_tol) || v < T(-zero_print_tol) || std::isnan(v) ?
           v : T(0);
  }

  //! \brief Print the vector \b X to the stream \a s.
  template<class T> std::ostream& operator<<(std::ostream& s,
                                             const vector<T>& X)
  {
    if (X.size() < 1)
      s <<" (empty)";
    else for (size_t i = 0; i < X.size(); i++)
      s << (i%nval_per_line ? ' ':'\n') << trunc(X[i]);

    return s << std::endl;
  }

  //! \brief Print the matrix \b A to the stream \a s.
  //! \details If the matrix is symmetric, only the upper triangular part of
  //! the matrix is printed, with the diagonal elements in the first column.
  //! The global variable \a zero_print_tol is used as a tolerance when
  //! checking whether the matrix is symmetric or not.
  template<class T> std::ostream& operator<<(std::ostream& s,
                                             const matrix<T>& A)
  {
    if (A.rows() < 1 || A.cols() < 1)
      return s <<" (empty)"<< std::endl;

    bool symm = A.isSymmetric(zero_print_tol);
    for (size_t i = 1; i <= A.rows(); i++)
    {
      size_t c1 = symm ? i : 1;
      s <<"\nRow "<< i <<": "<< trunc(A(i,c1));
      for (size_t j = c1+1; j <= A.cols(); j++)
        s <<' '<< trunc(A(i,j));
    }

    return s << std::endl;
  }

  //! \brief Print the 3D matrix \b A to the stream \a s.
  //! \details The matrix is priinted as a set of 2D sub-matrices
  //! based on the first two indices.
  template<class T> std::ostream& operator<<(std::ostream& s,
					     const matrix3d<T>& A)
  {
    if (A.empty())
      return s <<" (empty)"<< std::endl;

    matrix<T> B(A.dim(1),A.dim(2));
    const T* Aptr = A.ptr();
    for (size_t k = 0; k < A.dim(3); k++, Aptr += B.size())
    {
      B.fill(Aptr);
      if (k == 0) s <<"\n";
      s <<"i3="<< k+1 <<":"<< B;
    }

    return s;
  }

  //! \brief Print the vector \b X to the stream \a s in matlab format.
  template<class T> void writeMatlab(const char* label, const vector<T>& X,
                                     std::ostream& s = std::cout)
  {
    if (label)
      s << label <<" = [";
    else
      s <<"[";

    for (size_t i = 1; i <= X.size(); i++)
      s <<' '<< trunc(X(i));
    s <<" ];"<< std::endl;
  }

  //! \brief Print the matrix \b A to the stream \a s in matlab format.
  template<class T> void writeMatlab(const char* label, const matrix<T>& A,
                                     std::ostream& s = std::cout)
  {
    if (label)
      s << label <<" = [";
    else
      s <<"[";

    size_t nsp = 4 + strlen(label);
    for (size_t i = 1; i <= A.rows(); i++)
    {
      if (i > 1)
      {
        s <<";\n";
        for (size_t k = 0; k < nsp; k++) s <<' ';
      }
      for (size_t j = 1; j <= A.cols(); j++)
        s <<' '<< trunc(A(i,j));
    }
    s <<" ];"<< std::endl;
  }
}

#undef THIS
#endif
