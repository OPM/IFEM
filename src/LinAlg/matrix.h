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
//! The multiplication methods are implemented based on the BLAS library if
//! an implementation is provided (HAS_BLAS is defined, see BLAS.h).
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
#include "BLAS.h"

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
    explicit vector<T>(size_t n) { this->resize(n); }
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
    void fill(T s) { std::fill(this->begin(),this->end(),s); }
    //! \brief Fill the vector with data from an array.
    void fill(const T* values, size_t n = 0)
    {
      if (n > this->size()) this->std::vector<T>::resize(n);
      memcpy(&this->front(),values,this->size()*sizeof(T));
    }

    //! \brief Multiplication with a scalar.
    vector<T>& operator*=(T c);
    //! \brief Division by a scalar.
    vector<T>& operator/=(T d) { return this->operator*=(T(1)/d); }

    //! \brief Component-wise multiplication with a vector.
    vector<T>& operator*=(const std::vector<T>& X)
    {
      for (size_t i = 0; i < this->size() && i < X.size(); i++)
        this->operator[](i) *= X[i];
      return *this;
    }
    //! \brief Component-wise division with a vector.
    vector<T>& operator/=(const std::vector<T>& X)
    {
      for (size_t i = 0; i < this->size() && i < X.size(); i++)
        this->operator[](i) *= (X[i] == T(0) ? T(0) : T(1)/X[i]);
      return *this;
    }

    //! \brief Add the given vector \b X to \a *this.
    vector<T>& operator+=(const vector<T>& X) { return this->add(X); }
    //! \brief Subtract the given vector \b X from \a *this.
    vector<T>& operator-=(const vector<T>& X) { return this->add(X,T(-1)); }
    //! \brief Add the given vector \b X scaled by \a alfa to \a *this.
    vector<T>& add(const std::vector<T>& X, const T& alfa = T(1));

    //! \brief Perform \f${\bf Y} = \alpha{\bf Y} + (1-\alpha){\bf X} \f$
    //! where \b Y = \a *this.
    vector<T>& relax(T alfa, const std::vector<T>& X)
    {
      if (alfa != T(1))
      {
        this->operator*=(alfa);
        this->add(X,T(1)-alfa);
      }
      return *this;
    }
    //! \brief Perform \f${\bf Z} = \alpha{\bf Y} + (1-\alpha){\bf X} \f$
    //! where \b Z = \a *this.
    vector<T>& relax(T alfa, const std::vector<T>& X, const std::vector<T>& Y)
    {
      return this->operator=(Y).relax(alfa,X);
    }

    //! \brief Dot product between \a *this and another vector.
    //! \param[in] v The vector to dot this vector with
    //! \param[in] nv Length of the vector \a v
    //! \param[in] off1 Offset for this vector
    //! \param[in] inc1 Increment for this vector
    //! \param[in] off2 Offset for vector \b v
    //! \param[in] inc2 Increment for vector \b v
    T dot(const T* v, size_t nv,
          size_t off1 = 0, int inc1 = 1,
          size_t off2 = 0, int inc2 = 1) const;

    //! \brief Dot product between \a *this and another vector.
    //! \param[in] v The vector to dot this vector with
    //! \param[in] off1 Offset for this vector
    //! \param[in] inc1 Increment for this vector
    //! \param[in] off2 Offset for vector \b v
    //! \param[in] inc2 Increment for vector \b v
    T dot(const std::vector<T>& v,
          size_t off1 = 0, int inc1 = 1,
          size_t off2 = 0, int inc2 = 1) const
    {
      return this->dot(&v.front(),v.size(),off1,inc1,off2,inc2);
    }

    //! \brief Return the Euclidean norm of the vector.
    //! \param[in] off Index offset relative to the first vector component
    //! \param[in] inc Increment in the vector component indices
    T norm2(size_t off = 0, int inc = 1) const;
    //! \brief Return the infinite norm of the vector.
    //! \param off Index offset relative to the first vector component on input,
    //! 1-based index of the largest vector component on output
    //! \param[in] inc Increment in the vector component indices
    T normInf(size_t& off, int inc = 1) const;
    //! \brief Return the infinite norm of the vector (no index offset).
    //! \param[in] inc Increment in the vector component indices
    T normInf(int inc = 1) const { size_t o = 0; return this->normInf(o,inc); }

    //! \brief Return the largest element of the vector.
    T max() const { return *std::max_element(this->begin(),this->end()); }
    //! \brief Return the smallest element of the vector.
    T min() const { return *std::min_element(this->begin(),this->end()); }

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
      if (inc < 1 || this->empty())
        return xsum;

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
    \brief Common base class for multi-dimensional (2D and 3D) matrices.
    \details Contains the methods that the 2D and 3D matrices have in common.
  */

  template<class T> class matrixBase
  {
  protected:
    //! \brief The constructor is protected to allow sub-class instances only.
    matrixBase() { n[0] = n[1] = n[2] = 0; }
    //! \brief Constructor creating a matrix of dimension
    //! \f$n_1 \times n_2 \times n_3\f$.
    matrixBase(size_t n_1, size_t n_2, size_t n_3) : elem(n_1*n_2*n_3)
    {
      n[0] = n_1;
      n[1] = n_2;
      n[2] = n_3;
    }
    //! \brief Copy constructor, only copy the dimension of \a mat.
    matrixBase(const matrixBase<T>& mat) : elem(mat.size())
    {
      memcpy(n,mat.n,sizeof(n));
    }

    //! \brief Resize the matrix to dimension \f$n_1 \times n_2 \times n_3\f$.
    //! \details Will erase the previous content, but only if both the total
    //! matrix size, and any of its dimensions except the last are changed.
    //! It is therefore possible to add or remove a given number of elements in
    //! the last dimension of the matrix without loosing the contents of the
    //! other dimensions.
    //! If \a forceClear is \e true, the old matrix content is always erased.
    void redim(size_t n_1, size_t n_2, size_t n_3, bool forceClear)
    {
      if (forceClear)
      {
        // Erase previous content
        if (this->size() == n_1*n_2*n_3)
          this->fill(T(0));
        else
          this->clear();
      }

      if (n[0] == n_1 && n[1] == n_2 && n[2] == n_3) return; // nothing to do

      size_t oldn1 = n[0];
      size_t oldn2 = n[1];
      size_t oldSize = this->size();
      n[0] = n_1;
      n[1] = n_2;
      n[2] = n_3;
      if (this->size() == oldSize) return; // no more to do, size is unchanged

      // If the size in any of the matrix dimensions, except for the last one,
      // are changed the previous matrix content must be cleared
      if (!forceClear) this->clearIfNrowChanged(oldn1,oldn2);

      elem.std::template vector<T>::resize(n[0]*n[1]*n[2],T(0));
    }

    //! \brief Clears the matrix content if the first dimension(s) changed.
    virtual void clearIfNrowChanged(size_t n1, size_t n2) = 0;

  public:
    //! \brief Query dimensions.
    size_t dim(short int d = 1) const { return d > 0 && d <= 3 ? n[d-1] : 0; }
    //! \brief Query total matrix size.
    size_t size() const { return n[0]*n[1]*n[2]; }
    //! \brief Check if the matrix is empty.
    bool empty() const { return elem.empty(); }

    //! \brief Type casting to a one-dimensional vector.
    operator const vector<T>&() const { return elem; }

    //! \brief Access through pointer.
    T* ptr(size_t c = 0) { return &elem[n[0]*c]; }
    //! \brief Reference through pointer.
    const T* ptr(size_t c = 0) const { return &elem[n[0]*c]; }

    //! \brief Iterator to the start of the matrix elements.
    typename std::vector<T>::iterator begin() { return elem.begin(); }
    //! \brief Iterator to the end of the matrix elements.
    typename std::vector<T>::iterator end() { return elem.end(); }

    //! \brief Clears the matrix and sets its dimension to zero.
    void clear() { n[0] = n[1] = n[2] = 0; elem.clear(); }

    //! \brief Fill the matrix with a scalar value.
    void fill(T s) { std::fill(elem.begin(),elem.end(),s); }
    //! \brief Fill the matrix with data from an array.
    void fill(const T* values, size_t n = 0) { elem.fill(values,n); }

    //! \brief Add the given matrix \b A scaled by \a alfa to \a *this.
    matrixBase<T>& add(const matrixBase<T>& A, const T& alfa);
    //! \brief Multiplication of this matrix by a scalar \a c.
    matrixBase<T>& multiply(const T& c);

    //! \brief Return the Euclidean norm of the matrix.
    //! \param[in] inc Increment in the matrix element vector indices
    T norm2(int inc = 1) const { return elem.norm2(0,inc); }
    //! \brief Return the sum of the absolute value of the matrix elements.
    //! \param[in] inc Increment in the matrix element vector indices
    T asum(int inc = 1) const { return elem.asum(0,inc); }
    //! \brief Return the sum of the matrix elements.
    //! \param[in] inc Increment in the matrix element vector indices
    T sum(int inc = 1) const { return elem.sum(0,inc); }

  protected:
    size_t    n[3]; //!< Dimension of the matrix
    vector<T> elem; //!< Actual matrix elements, stored column by column
  };


  /*!
    \brief Two-dimensional rectangular matrix with some algebraic operations.
    \details This is a 2D equivalent to the \a vector class. The matrix elements
    are stored column-wise in a one-dimensional array, such that its pointer
    might be passed to Fortran subroutines requiring 2D arrays as arguments.
  */

  template<class T> class matrix : public matrixBase<T>
  {
  public:
    //! \brief Constructor creating an empty matrix.
    matrix() : nrow(this->n[0]), ncol(this->n[1]) {}
    //! \brief Constructor creating a matrix of dimension \f$r \times c\f$.
    matrix(size_t r, size_t c)
      : matrixBase<T>(r,c,1), nrow(this->n[0]), ncol(this->n[1]) {}
    //! \brief Copy constructor, optionally creates the transpose of \a mat.
    matrix(const matrix<T>& mat, bool transposed = false)
      : matrixBase<T>(mat), nrow(this->n[0]), ncol(this->n[1])
    {
      nrow = transposed ? mat.ncol : mat.nrow;
      ncol = transposed ? mat.nrow : mat.ncol;
      if (transposed)
        for (size_t r = 0; r < ncol; r++)
          for (size_t c = 0; c < nrow; c++)
            this->elem[c+nrow*r] = mat.elem[r+ncol*c];
      else
        this->elem.fill(mat.elem.ptr());
    }
    //! \brief Empty destructor.
    virtual ~matrix() {}

    //! \brief Resize the matrix to dimension \f$r \times c\f$.
    //! \details Will erase the previous content, but only if both
    //! the total matrix size and the number of rows in the matrix are changed.
    //! It is therefore possible to add or remove a given number of columns to
    //! the matrix without loosing the contents of the remaining columns.
    //! If \a forceClear is \e true, the old matrix content is always erased.
    void resize(size_t r, size_t c, bool forceClear = false)
    {
      this->redim(r,c,1,forceClear);
    }

    //! \brief Increase or decrease the number of rows in the matrix.
    matrix<T>& expandRows(int incRows)
    {
      int newRows = nrow + incRows;
      if (newRows < 1 || ncol < 1)
        // The matrix is empty
        this->clear();
      else if (incRows < 0)
      {
        // The matrix size is reduced
        T* newMat = this->ptr() + newRows;
        for (size_t c = 1; c < ncol; c++, newMat += newRows)
          memmove(newMat,this->ptr(c),newRows*sizeof(T));
        nrow = newRows;
        this->elem.std::template vector<T>::resize(nrow*ncol);
      }
      else if (incRows > 0)
      {
        // The matrix size is increased
        size_t oldRows = nrow;
        nrow = newRows;
        this->elem.std::template vector<T>::resize(nrow*ncol,T(0));
        T* oldMat = this->ptr() + oldRows*(ncol-1);
        for (size_t c = ncol-1; c > 0; c--, oldMat -= oldRows)
        {
          memmove(this->ptr(c),oldMat,oldRows*sizeof(T));
          for (size_t r = nrow-1; r >= oldRows; r--)
            this->elem[r+nrow*(c-1)] = T(0);
        }
      }
      return *this;
    }

    //! \brief Query number of matrix rows.
    size_t rows() const { return nrow; }
    //! \brief Query number of matrix columns.
    size_t cols() const { return ncol; }

    //! \brief Assignment operator.
    matrix<T>& operator=(const matrix<T>& A)
    {
      memcpy(this->n,A.n,sizeof(A.n));
      this->elem = A.elem;
      return *this;
    }

    //! \brief Overloaded assignment operator.
    matrix<T>& operator=(const std::vector<T>& X)
    {
      // Do not use vector<T>::operator= because we don't want to alter size
      size_t nval = X.size() < this->elem.size() ? X.size() : this->elem.size();
      std::copy(X.begin(),X.begin()+nval,this->elem.begin());
      std::fill(this->elem.begin()+nval,this->elem.end(),T(0));
      return *this;
    }

    //! \brief Index-1 based element access.
    //! \details Assuming column-wise storage as in Fortran.
    T& operator()(size_t r, size_t c)
    {
      CHECK_INDEX("matrix::operator(): Row-index ",r,nrow);
      CHECK_INDEX("matrix::operator(): Column-index ",c,ncol);
      return this->elem[r-1+nrow*(c-1)];
    }

    //! \brief Index-1 based element reference.
    //! \details Assuming column-wise storage as in Fortran.
    const T& operator()(size_t r, size_t c) const
    {
      CHECK_INDEX("matrix::operator(): Row-index ",r,nrow);
      CHECK_INDEX("matrix::operator(): Column-index ",c,ncol);
      return this->elem[r-1+nrow*(c-1)];
    }

    //! \brief Extract a row from the matrix.
    vector<T> getRow(size_t r) const
    {
      CHECK_INDEX("matrix::getRow: Row-index ",r,nrow);
      if (nrow < 2) return this->elem;
      vector<T> row(ncol);
      for (size_t i = 0; i < ncol; i++)
        row[i] = this->elem[r-1+nrow*i];
      return row;
    }

    //! \brief Extract a column from the matrix.
    vector<T> getColumn(size_t c) const
    {
      CHECK_INDEX("matrix::getColumn: Column-index ",c,ncol);
      if (ncol < 2) return this->elem;
      vector<T> col(nrow);
      memcpy(col.ptr(),this->ptr(c-1),nrow*sizeof(T));
      return col;
    }

    //! \brief Fill a column of the matrix.
    void fillColumn(size_t c, const std::vector<T>& data)
    {
      CHECK_INDEX("matrix::fillColumn: Column-index ",c,ncol);
      size_t ndata = nrow > data.size() ? data.size() : nrow;
      memcpy(this->ptr(c-1),&data.front(),ndata*sizeof(T));
    }

    //! \brief Fill a column of the matrix.
    void fillColumn(size_t c, const T* data)
    {
      CHECK_INDEX("matrix::fillColumn: Column-index ",c,ncol);
      memcpy(this->ptr(c-1),data,nrow*sizeof(T));
    }

    //! \brief Fill a row of the matrix.
    void fillRow(size_t r, const T* data)
    {
      CHECK_INDEX("matrix::fillRow: Row-index ",r,nrow);
      if (nrow < 2)
        this->elem.fill(data);
      else for (size_t i = 0; i < ncol; i++)
        this->elem[r-1+nrow*i] = data[i];
    }

    //! \brief Fill a block of the matrix with another matrix.
    void fillBlock(const matrix<T>& block, size_t r, size_t c,
                   bool transposed = false)
    {
      size_t nr = transposed ? block.cols() : block.rows();
      size_t nc = transposed ? block.rows() : block.cols();
      for (size_t i = 1; i <= nr && i+r-1 <= nrow; i++)
      {
        size_t ip = i+r-2 + nrow*(c-1);
        for (size_t j = 1; j <= nc && j+c-1 <= ncol; j++, ip += nrow)
          this->elem[ip] = transposed ? block(j,i) : block(i,j);
      }
    }

    //! \brief Add a scalar multiple of another matrix to a block of the matrix.
    void addBlock(const matrix<T>& block, T s, size_t r, size_t c,
                  bool transposed = false)
    {
      size_t nr = transposed ? block.cols() : block.rows();
      size_t nc = transposed ? block.rows() : block.cols();
      for (size_t i = 1; i <= nr && i+r-1 <= nrow; i++)
      {
        size_t ip = i+r-2 + nrow*(c-1);
        for (size_t j = 1; j <= nc && j+c-1 <= ncol; j++, ip += nrow)
          this->elem[ip] += s*(transposed ? block(j,i) : block(i,j));
      }
    }

    //! \brief Create a diagonal matrix.
    matrix<T>& diag(T d, size_t dim = 0)
    {
      if (dim > 0)
        this->resize(dim,dim,true);
      else
        this->resize(nrow,ncol,true);
      for (size_t r = 0; r < nrow && r < ncol; r++)
        this->elem[r+nrow*r] = d;
      return *this;
    }

    //! \brief Replace the current matrix by its transpose.
    matrix<T>& transpose()
    {
      matrix<T> tmp(*this);
      for (size_t r = 0; r < nrow; r++)
        for (size_t c = 0; c < ncol; c++)
          this->elem[c+ncol*r] = tmp.elem[r+nrow*c];

      nrow = tmp.ncol;
      ncol = tmp.nrow;
      return *this;
    }

    //! \brief Return the trace of the matrix (sum of its diagonal elements).
    T trace() const { return this->elem.sum(0,nrow+1); }

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
    T inverse(T tol = T(0))
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
    bool isSymmetric(T tol = T(0)) const
    {
      if (nrow != ncol) return false;

      for (size_t r = 0; r < nrow; r++)
        for (size_t c = 0; c < r; c++)
        {
          T diff = this->elem[r+nrow*c] - this->elem[c+nrow*r];
          if (diff < -tol || diff > tol) return false;
        }

      return true;
    }

    //! \brief Add the given matrix \b A to \a *this.
    matrix<T>& operator+=(const matrix<T>& A) { return this->add(A); }
    //! \brief Subtract the given matrix \b A from \a *this.
    matrix<T>& operator-=(const matrix<T>& A) { return this->add(A,T(-1)); }
    //! \brief Add the given matrix \b A scaled by \a alfa to \a *this.
    matrix<T>& add(const matrix<T>& A, T alfa = T(1))
    {
      return static_cast<matrix<T>&>(this->matrixBase<T>::add(A,alfa));
    }

    //! \brief Multiplication with a scalar.
    matrix<T>& operator*=(T c) { return this->multiply(c); }
    //! \brief Division by a scalar.
    matrix<T>& operator/=(T d) { return this->multiply(T(1)/d); }
    //! \brief Multiplication of this matrix by a scalar \a c.
    matrix<T>& multiply(T c)
    {
      return static_cast<matrix<T>&>(this->matrixBase<T>::multiply(c));
    }

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
                        bool addTo = false, const T& alpha = T(1));

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
      -# \f$ {\bf Y} = {\bf Y} - {\bf A} {\bf X} \f$
      -# \f$ {\bf Y} = {\bf Y} - {\bf A}^T {\bf X} \f$
    */
    bool multiply(const std::vector<T>& X, std::vector<T>& Y,
                  bool transA = false, char addTo = 0) const;

    /*! \brief Matrix-vector multiplication.
      \details Performs the following operations (\b A = \a *this):
      -# \f$ {\bf Y} = {\alpha}{\bf A} {\bf X} + {\beta}{\bf Y}\f$
      -# \f$ {\bf Y} = {\alpha}{\bf A}^T {\bf X} + {\beta}{\bf Y}\f$
    */
    bool multiply(const std::vector<T>& X, std::vector<T>& Y,
                  const T& alpha, const T& beta = T(0),
                  bool transA = false, int stridex = 1, int stridey = 1,
                  int ofsx = 0, int ofsy = 0) const;

    //! \brief Outer product between two vectors.
    bool outer_product(const std::vector<T>& X, const std::vector<T>& Y,
                       bool addTo = false, T alpha = T(1));

    //! \brief Return the infinite norm of the matrix.
    T normInf() const
    {
      if (nrow == 0) return T(0);

      // Compute row sums
      vector<T> sums(nrow);
      for (size_t i = 0; i < nrow; i++)
        sums[i] = this->elem.asum(i,ncol);
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

    //! \brief Check dimension compatibility for outer product multiplication.
    bool compatible(const std::vector<T>& X, const std::vector<T>& Y)
    {
      if (X.size() == nrow && Y.size() == ncol) return true;

      std::cerr <<"matrix::outer_product: Incompatible matrix and vectors: A("
                << nrow <<','<< ncol <<"), X("
                << X.size() <<"), Y("<< Y.size() <<")\n"
                <<"                       when computing A += X*Y^t"
                << std::endl;
      ABORT_ON_INDEX_CHECK;
      return false;
    }

  protected:
    //! \brief Clears the content if the number of rows changed.
    virtual void clearIfNrowChanged(size_t n1, size_t)
    {
      if (n1 != nrow) this->elem.clear();
    }

  private:
    size_t& nrow; //!< Number of matrix rows
    size_t& ncol; //!< Number of matrix columns
  };


#ifdef HAS_BLAS
  //============================================================================
  //===   BLAS-implementation of the matrix/vector multiplication methods   ====
  //============================================================================

  template<> inline
  vector<float>& vector<float>::operator*=(float c)
  {
    cblas_sscal(this->size(),c,this->ptr(),1);
    return *this;
  }

  template<> inline
  vector<double>& vector<double>::operator*=(double c)
  {
    cblas_dscal(this->size(),c,this->ptr(),1);
    return *this;
  }

  template<> inline
  float vector<float>::dot(const float* v, size_t nv,
                           size_t o1, int i1, size_t o2, int i2) const
  {
    int n1 = i1 > 1 || i1 < -1 ? this->size()/abs(i1) : this->size()-o1;
    int n2 = i2 > 1 || i2 < -1 ? nv/abs(i2) : nv-o2;
    int n  = n1 < n2 ? n1 : n2;
    return cblas_sdot(n,this->ptr()+o1,i1,v+o2,i2);
  }

  template<> inline
  double vector<double>::dot(const double* v, size_t nv,
                             size_t o1, int i1, size_t o2, int i2) const
  {
    int n1 = i1 > 1 || i1 < -1 ? this->size()/abs(i1) : this->size()-o1;
    int n2 = i2 > 1 || i2 < -1 ? nv/abs(i2) : nv-o2;
    int n  = n1 < n2 ? n1 : n2;
    return cblas_ddot(n,this->ptr()+o1,i1,v+o2,i2);
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
    if (inc < 1 || this->empty())
      return 0.0f;

    const float* v = this->ptr() + off;
    off = 1 + cblas_isamax(this->size()/inc,v,inc);
    return fabsf(v[(off-1)*inc]);
  }

  template<> inline
  double vector<double>::normInf(size_t& off, int inc) const
  {
    if (inc < 1 || this->empty())
      return 0.0;

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
  vector<float>& vector<float>::add(const std::vector<float>& X,
                                    const float& alfa)
  {
    size_t n = this->size() < X.size() ? this->size() : X.size();
    cblas_saxpy(n,alfa,&X.front(),1,this->ptr(),1);
    return *this;
  }

  template<> inline
  vector<double>& vector<double>::add(const std::vector<double>& X,
                                      const double& alfa)
  {
    size_t n = this->size() < X.size() ? this->size() : X.size();
    cblas_daxpy(n,alfa,&X.front(),1,this->ptr(),1);
    return *this;
  }

  template<> inline
  matrixBase<float>& matrixBase<float>::add(const matrixBase<float>& A,
                                            const float& alfa)
  {
    size_t n = this->size() < A.size() ? this->size() : A.size();
    cblas_saxpy(n,alfa,A.ptr(),1,this->ptr(),1);
    return *this;
  }

  template<> inline
  matrixBase<double>& matrixBase<double>::add(const matrixBase<double>& A,
                                              const double& alfa)
  {
    size_t n = this->size() < A.size() ? this->size() : A.size();
    cblas_daxpy(n,alfa,A.ptr(),1,this->ptr(),1);
    return *this;
  }

  template<> inline
  matrixBase<float>& matrixBase<float>::multiply(const float& c)
  {
    cblas_sscal(this->size(),c,this->ptr(),1);
    return *this;
  }

  template<> inline
  matrixBase<double>& matrixBase<double>::multiply(const double& c)
  {
    cblas_dscal(this->size(),c,this->ptr(),1);
    return *this;
  }

  template<> inline
  bool matrix<float>::multiply(const std::vector<float>& X,
                               std::vector<float>& Y,
                               bool transA, char addTo) const
  {
    if (!this->compatible(X,transA)) return false;
    if (!addTo) Y.resize(transA ? ncol : nrow);

    cblas_sgemv(CblasColMajor,
                transA ? CblasTrans : CblasNoTrans,
                nrow, ncol, addTo < 0 ? -1.0f : 1.0f,
                this->ptr(), nrow,
                &X.front(), 1, addTo ? 1.0f : 0.0f,
                &Y.front(), 1);

    return true;
  }

  template<> inline
  bool matrix<double>::multiply(const std::vector<double>& X,
                                std::vector<double>& Y,
                                bool transA, char addTo) const
  {
    if (!this->compatible(X,transA)) return false;
    if (!addTo) Y.resize(transA ? ncol : nrow);

    cblas_dgemv(CblasColMajor,
                transA ? CblasTrans : CblasNoTrans,
                nrow, ncol, addTo < 0 ? -1.0 : 1.0,
                this->ptr(), nrow,
                &X.front(), 1, addTo ? 1.0 : 0.0,
                &Y.front(), 1);

    return true;
  }

  template<> inline
  bool matrix<float>::multiply(const std::vector<float>& X,
                               std::vector<float>& Y,
                               const float& alpha, const float& beta,
                               bool transA, int stridex, int stridey,
                               int ofsx, int ofsy) const
  {
    if (ofsx == 0 && stridex == 1)
      if (!this->compatible(X,transA)) return false;

    if (beta == 0.0f)
    {
      Y.resize(ofsy + 1 + ((transA ? ncol : nrow)-1)*stridey);
      if (stridey > 0) std::fill(Y.begin(),Y.end(),0.0f);
    }

    cblas_sgemv(CblasColMajor,
                transA ? CblasTrans : CblasNoTrans,
                nrow, ncol, alpha,
                this->ptr(), nrow,
                &X[ofsx], stridex, beta,
                &Y[ofsy], stridey);

    return true;
  }

  template<> inline
  bool matrix<double>::multiply(const std::vector<double>& X,
                                std::vector<double>& Y,
                                const double& alpha, const double& beta,
                                bool transA, int stridex, int stridey,
                                int ofsx, int ofsy) const
  {
    if (ofsx == 0 && stridex == 1)
      if (!this->compatible(X,transA)) return false;

    if (beta == 0.0)
    {
      Y.resize(ofsy + 1 + ((transA ? ncol : nrow)-1)*stridey);
      if (stridey > 0) std::fill(Y.begin(),Y.end(),0.0);
    }

    cblas_dgemv(CblasColMajor,
                transA ? CblasTrans : CblasNoTrans,
                nrow, ncol, alpha,
                this->ptr(), nrow,
                &X[ofsx], stridex, beta,
                &Y[ofsy], stridey);

    return true;
  }

  template<> inline
  matrix<float>& matrix<float>::multiply(const matrix<float>& A,
                                         const matrix<float>& B,
                                         bool transA, bool transB,
                                         bool addTo, const float& alpha)
  {
    size_t M, N, K;
    if (!this->compatible(A,B,transA,transB,M,N,K))
      this->clear();
    else if (!addTo)
      this->resize(M,N);

    if (!this->empty())
      cblas_sgemm(CblasColMajor,
                  transA ? CblasTrans : CblasNoTrans,
                  transB ? CblasTrans : CblasNoTrans,
                  M, N, K, alpha,
                  A.ptr(), A.nrow,
                  B.ptr(), B.nrow,
                  addTo ? 1.0f : 0.0f,
                  this->ptr(), nrow);

    return *this;
  }

  template<> inline
  matrix<double>& matrix<double>::multiply(const matrix<double>& A,
                                           const matrix<double>& B,
                                           bool transA, bool transB,
                                           bool addTo, const double& alpha)
  {
    size_t M, N, K;
    if (!this->compatible(A,B,transA,transB,M,N,K))
      this->clear();
    else if (!addTo)
      this->resize(M,N);

    if (!this->empty())
      cblas_dgemm(CblasColMajor,
                  transA ? CblasTrans : CblasNoTrans,
                  transB ? CblasTrans : CblasNoTrans,
                  M, N, K, alpha,
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
  bool matrix<float>::outer_product(const std::vector<float>& X,
                                    const std::vector<float>& Y,
                                    bool addTo, float alpha)
  {
    if (!addTo)
      this->resize(X.size(),Y.size());
    else if (!this->compatible(X,Y))
      return false;

    cblas_sgemm(CblasColMajor,
                CblasNoTrans, CblasTrans,
                nrow, ncol, 1, alpha,
                X.data(), nrow,
                Y.data(), ncol,
                addTo ? 1.0f : 0.0f,
                this->ptr(), nrow);

    return true;
  }

  template<> inline
  bool matrix<double>::outer_product(const std::vector<double>& X,
                                     const std::vector<double>& Y,
                                     bool addTo, double alpha)
  {
    if (!addTo)
      this->resize(X.size(),Y.size());
    else if (!this->compatible(X,Y))
      return false;

    cblas_dgemm(CblasColMajor,
                CblasNoTrans, CblasTrans,
                nrow, ncol, 1, alpha,
                X.data(), nrow,
                Y.data(), ncol,
                addTo ? 1.0 : 0.0,
                this->ptr(), nrow);

    return true;
  }

#else
  //============================================================================
  //===   Non-BLAS inlined implementations (slow...)   =========================
  //============================================================================

  template<class T> inline
  vector<T>& vector<T>::operator*=(T c)
  {
    for (size_t i = 0; i < this->size(); i++)
      std::vector<T>::operator[](i) *= c;
    return *this;
  }

  template<class T> inline
  T vector<T>::dot(const T* v, size_t nv,
                   size_t o1, int i1, size_t o2, int i2) const
  {
    size_t i, j;
    T dotprod = T(0);
    for (i = o1, j = o2; i < this->size() && j < nv; i += i1, j += i2)
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
  vector<T>& vector<T>::add(const std::vector<T>& X, const T& alfa)
  {
    T* p = this->ptr();
    const T* q = &X.front();
    for (size_t i = 0; i < this->size() && i < X.size(); i++, p++, q++)
      *p += alfa*(*q);
    return *this;
  }

  template<class T> inline
  matrixBase<T>& matrixBase<T>::add(const matrixBase<T>& A, const T& alfa)
  {
    T* p = this->ptr();
    const T* q = A.ptr();
    for (size_t i = 0; i < this->size() && i < A.size(); i++, p++, q++)
      *p += alfa*(*q);
    return *this;
  }

  template<class T> inline
  matrixBase<T>& matrixBase<T>::multiply(const T& c)
  {
    for (size_t i = 0; i < this->elem.size(); i++)
      this->elem[i] *= c;
    return *this;
  }

  template<class T> inline
  bool matrix<T>::multiply(const std::vector<T>& X, std::vector<T>& Y,
                           bool transA, char addTo) const
  {
    if (!this->compatible(X,transA)) return false;
    if (!addTo)
    {
      Y.clear();
      Y.resize(transA ? ncol : nrow, T(0));
    }

    for (size_t i = 0; i < Y.size(); i++)
      for (size_t j = 0; j < X.size(); j++)
        if (transA)
          Y[i] += THIS(j+1,i+1) * (addTo < 0 ? -X[j] : X[j]);
        else
          Y[i] += THIS(i+1,j+1) * (addTo < 0 ? -X[j] : X[j]);

    return true;
  }

  template<class T> inline
  bool matrix<T>::multiply(const std::vector<T>& X, std::vector<T>& Y,
                           const T& alpha, const T& beta,
                           bool transA, int stridex, int stridey,
                           int ofsx, int ofsy) const
  {
    if (ofsx == 0 && stridex == 1)
      if (!this->compatible(X,transA)) return false;

    if (beta == T(0))
    {
      Y.resize(ofsy + 1 + ((transA ? ncol : nrow)-1)*stridey);
      std::fill(Y.begin(),Y.end(),T(0));
    }
    else for (size_t i = ofsy; i < Y.size(); i += stridey)
      Y[i] *= beta;

    size_t a, b, i, j;
    for (a = 1, i = ofsy; i < Y.size(); a++, i += stridey)
      for (b = 1, j = ofsx; j < X.size(); b++, j += stridex)
        if (transA)
          Y[i] += alpha * THIS(b,a) * X[j];
        else
          Y[i] += alpha * THIS(a,b) * X[j];

    return true;
  }

  template<class T> inline
  matrix<T>& matrix<T>::multiply(const matrix<T>& A,
                                 const matrix<T>& B,
                                 bool transA, bool transB, bool addTo,
                                 const T& alpha)
  {
    size_t M, N, K;
    if (!this->compatible(A,B,transA,transB,M,N,K))
    {
      this->clear();
      M = 0;
    }
    else if (!addTo)
      this->resize(M,N,true);

    for (size_t i = 1; i <= M; i++)
      for (size_t j = 1; j <= N; j++)
        for (size_t k = 1; k <= K; k++)
          if (transA && transB)
            THIS(i,j) += alpha*A(k,i)*B(j,k);
          else if (transA)
            THIS(i,j) += alpha*A(k,i)*B(k,j);
          else if (transB)
            THIS(i,j) += alpha*A(i,k)*B(j,k);
          else
            THIS(i,j) += alpha*A(i,k)*B(k,j);

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
  bool matrix<T>::outer_product(const std::vector<T>& X,
                                const std::vector<T>& Y,
                                bool addTo, T alpha)
  {
    if (!addTo)
      this->resize(X.size(),Y.size());
    else if (!this->compatible(X,Y))
      return false;

    if (addTo)
      for (size_t j = 0; j < ncol; j++)
        for (size_t i = 0; i < nrow; i++)
          this->elem[i+nrow*j] += alpha*X[i]*Y[j];
    else
      for (size_t j = 0; j < ncol; j++)
        for (size_t i = 0; i < nrow; i++)
          this->elem[i+nrow*j] = alpha*X[i]*Y[j];

    return true;
  }

#endif

  //============================================================================
  //===   Global operators   ===================================================
  //============================================================================

  extern double zero_print_tol; //!< Zero tolerance for printing numbers
  extern int    nval_per_line;  //!< Number of values to print per line

  //! \brief Truncate a value to zero when it is less than a given threshold.
  //! \details Used when printing matrices for easy comparison with other
  //! matrices when they contain terms that are numerically zero, except for
  //! some round-off noise. The value of the global variable \a zero_print_tol
  //! is used as a tolerance in this method.
  template<class T> inline T trunc(T v)
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
      s << ((i%nval_per_line) ? ' ':'\n') << trunc(X[i]);

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

    size_t nsp = label ? 4 + strlen(label) : 1;
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
