// $Id$
//==============================================================================
//!
//! \file matrixnd.h
//!
//! \date Oct 1 2007
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Simple template classes for dense rectangular 3d and 4d matrices.
//!
//==============================================================================

#ifndef UTL_MATRIXND_H
#define UTL_MATRIXND_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cmath>
#include "matrix.h"


namespace utl //! General utility classes and functions.
{
  /*!
    \brief Three-dimensional rectangular matrix with some algebraic operations.
    \details This is a 3D equivalent to the \a matrix class. The matrix elements
    are stored column-wise in a one-dimensional array, such that its pointer
    might be passed to Fortran subroutines requiring 3D arrays as arguments.
  */

  template<class T> class matrix3d : public matrixBase<T>
  {
  public:
    //! \brief Constructor creating an empty matrix.
    matrix3d() {}
    //! \brief Constructor creating a matrix of dimension
    //! \f$n_1 \times n_2 \times n_3\f$.
    matrix3d(size_t n_1, size_t n_2, size_t n_3) : matrixBase<T>(n_1,n_2,n_3,1) {}
    //! \brief Constructor to read a matrix from a stream.
    matrix3d(std::istream& is, std::streamsize max = 10)
    {
      // Read size
      size_t n0 = 0, n1 = 0, n2 = 0;
      is.ignore(10, ':');
      is.ignore(1);
      is >> n0 >> n1 >> n2;
      is.ignore(1, '\n');

      // Read contents
      this->resize(n0,n1,n2,1);
      for (size_t k = 0; k < n2; k++) {
        is.ignore(max, '\n');
        for (size_t i = 0; i < n0; i++) {
          is.ignore(max, ':');
          for (size_t j = 0; j < n1; j++)
            is >> this->elem[n0*(n1*k + j) + i];
          is.ignore(1, '\n');
        }
      }
    }

    //! \brief Empty destructor.
    virtual ~matrix3d() {}

    //! \brief Resize the matrix to dimension \f$n_1 \times n_2 \times n_3\f$.
    //! \details Will erase the previous content, but only if both
    //! the total matrix size, and n_1 or n_2 in the matrix are changed.
    //! It is therefore possible to add or remove a given number of elements in
    //! the third dimension of the matrix without loosing the contents of the
    //! first and second dimensions.
    //! If \a forceClear is \e true, the old matrix content is always erased.
    void resize(size_t n_1, size_t n_2, size_t n_3, bool forceClear = false)
    {
      this->redim(n_1,n_2,n_3,1,forceClear);
    }

    //! \brief Index-1 based element access.
    //! \details Assuming column-wise storage as in Fortran.
    T& operator()(size_t i1, size_t i2, size_t i3)
    {
      CHECK_INDEX("matrix3d::operator(): First index " ,i1,this->n[0]);
      CHECK_INDEX("matrix3d::operator(): Second index ",i2,this->n[1]);
      CHECK_INDEX("matrix3d::operator(): Third index " ,i3,this->n[2]);
      return this->elem[i1-1 + this->n[0]*((i2-1) + this->n[1]*(i3-1))];
    }

    //! \brief Index-1 based element access.
    //! \details Assuming column-wise storage as in Fortran.
    const T& operator()(size_t i1, size_t i2, size_t i3) const
    {
      CHECK_INDEX("matrix3d::operator(): First index " ,i1,this->n[0]);
      CHECK_INDEX("matrix3d::operator(): Second index ",i2,this->n[1]);
      CHECK_INDEX("matrix3d::operator(): Third index " ,i3,this->n[2]);
      return this->elem[i1-1 + this->n[0]*((i2-1) + this->n[1]*(i3-1))];
    }

    //! \brief Extract a column from the matrix.
    vector<T> getColumn(size_t i2, size_t i3) const
    {
      CHECK_INDEX("matrix3d::getColumn: Second index ",i2,this->n[1]);
      CHECK_INDEX("matrix3d::getColumn: Third index " ,i3,this->n[2]);
      if (this->n[1] < 2 && this->n[2] < 2) return this->elem;
      vector<T> col(this->n[0]);
      memcpy(col.ptr(),this->ptr(i2-1+this->n[1]*(i3-1)),col.size()*sizeof(T));
      return col;
    }

    //! \brief Fill a column of the matrix.
    void fillColumn(size_t i2, size_t i3, const std::vector<T>& data)
    {
      CHECK_INDEX("matrix3d::fillColumn: Second index ",i2,this->n[1]);
      CHECK_INDEX("matrix3d::fillColumn: Third index " ,i3,this->n[2]);
      size_t ndata = this->n[0] > data.size() ? data.size() : this->n[0];
      memcpy(this->ptr(i2-1+this->n[1]*(i3-1)),&data.front(),ndata*sizeof(T));
    }

    //! \brief Return the trace of the \a i1'th sub-matrix.
    T trace(size_t i1) const
    {
      CHECK_INDEX("matrix3d::trace(): Index ",i1,this->n[0]);
      return this->elem.sum(i1-1,this->n[0]*(this->n[1]+1));
    }

    //! \brief Add the given matrix \b X to \a *this.
    matrix3d<T>& operator+=(const matrix3d<T>& X) { return this->add(X); }
    //! \brief Subtract the given matrix \b X from \a *this.
    matrix3d<T>& operator-=(const matrix3d<T>& X) { return this->add(X,T(-1)); }
    //! \brief Add the given matrix \b X scaled by \a alfa to \a *this.
    matrix3d<T>& add(const matrix3d<T>& X, T alfa = T(1))
    {
      return static_cast<matrix3d<T>&>(this->matrixBase<T>::add(X,alfa));
    }

    //! \brief Multiplication of this matrix by a scalar \a c.
    matrix3d<T>& multiply(T c)
    {
      return static_cast<matrix3d<T>&>(this->matrixBase<T>::multiply(c));
    }

    /*! \brief Matrix-matrix multiplication.
      \details Performs one of the following operations (\b C = \a *this):
      -# \f$ {\bf C} = {\bf A} {\bf B} \f$
      -# \f$ {\bf C} = {\bf A}^T {\bf B} \f$
      -# \f$ {\bf C} = {\bf C} + {\bf A} {\bf B} \f$
      -# \f$ {\bf C} = {\bf C} + {\bf A}^T {\bf B} \f$
    */
    bool multiply(const matrix<T>& A, const matrix3d<T>& B,
                  bool transA = false, bool addTo = false);

    /*! \brief Matrix-matrix multiplication.
      \details Performs one of the following operations (\b C = \a *this):
      -# \f$ {\bf C} = {\bf A} {\bf B} \f$
      -# \f$ {\bf C} = {\bf C} + {\bf A} {\bf B} \f$

      The matrix \b A is here represented by a one-dimensional vector, and its
      number of columns is assumed to match the first dimension of \b B and its
      number of rows is then the total vector length divided by
      the number of columns.
    */
    bool multiplyMat(const std::vector<T>& A, const matrix3d<T>& B,
                     bool addTo = false);

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

    //! \brief Check dimension compatibility for matrix-matrix multiplication,
    //! when the matrix A is represented by a one-dimensional vector.
    bool compatible(const std::vector<T>& A, const matrix3d<T>& B,
                    size_t& M, size_t& N, size_t& K)
    {
      N = B.n[1]*B.n[2];
      K = B.n[0];
      M = K > 0 ? A.size() / K : 0;
      if (M*K == A.size() && !A.empty()) return true;

      std::cerr <<"matrix3d::multiply: Incompatible matrices: A(r*c="<< A.size()
                <<"), B("<< B.n[0] <<","<< B.n[1] <<","<< B.n[2] <<")\n"
                <<"                  when computing C = A * B"<< std::endl;
      ABORT_ON_INDEX_CHECK;
      return false;
    }

  protected:
    //! \brief Clears the content if any of the first two dimensions changed.
    virtual void clearIfNrowChanged(size_t n1, size_t n2, size_t)
    {
      if (n1 != this->n[0] || n2 != this->n[1]) this->elem.clear();
    }
  };


  /*!
    \brief Four-dimensional rectangular matrix with some algebraic operations.
    \details This is a 4D equivalent to the \a matrix class. The matrix elements
    are stored column-wise in a one-dimensional array, such that its pointer
    might be passed to Fortran subroutines requiring 4D arrays as arguments.
  */

  template<class T> class matrix4d : public matrixBase<T>
  {
  public:
    //! \brief Constructor creating an empty matrix.
    matrix4d() {}
    //! \brief Constructor creating a matrix of dimension
    //! \f$n_1 \times n_2 \times n_3 \times n_4\f$.
    matrix4d(size_t n_1, size_t n_2, size_t n_3, size_t n_4) :
      matrixBase<T>(n_1,n_2,n_3,n_4) {}

    //! \brief Empty destructor.
    virtual ~matrix4d() {}

    //! \brief Resize the matrix to dimension \f$n_1 \times n_2 \times n_3 \times n_4\f$.
    //! \details Will erase the previous content, but only if both
    //! the total matrix size, and n_1 or n_2 in the matrix are changed.
    //! It is therefore possible to add or remove a given number of elements in
    //! the third dimension of the matrix without loosing the contents of the
    //! first and second dimensions.
    //! If \a forceClear is \e true, the old matrix content is always erased.
    void resize(size_t n_1, size_t n_2, size_t n_3,
                size_t n_4, bool forceClear = false)
    {
      this->redim(n_1,n_2,n_3,n_4,forceClear);
    }

    //! \brief Index-1 based element access.
    //! \details Assuming column-wise storage as in Fortran.
    T& operator()(size_t i1, size_t i2, size_t i3, size_t i4)
    {
      CHECK_INDEX("matrix4d::operator(): First index " ,i1,this->n[0]);
      CHECK_INDEX("matrix4d::operator(): Second index ",i2,this->n[1]);
      CHECK_INDEX("matrix4d::operator(): Third index " ,i3,this->n[2]);
      CHECK_INDEX("matrix4d::operator(): Third index " ,i4,this->n[3]);
      return this->elem[i1-1 + this->n[0]*((i2-1) + this->n[1]*((i3-1) + this->n[2]*(i4-1)))];
    }

    //! \brief Index-1 based element access.
    //! \details Assuming column-wise storage as in Fortran.
    const T& operator()(size_t i1, size_t i2, size_t i3, size_t i4) const
    {
      CHECK_INDEX("matrix3d::operator(): First index " ,i1,this->n[0]);
      CHECK_INDEX("matrix3d::operator(): Second index ",i2,this->n[1]);
      CHECK_INDEX("matrix3d::operator(): Third index " ,i3,this->n[2]);
      CHECK_INDEX("matrix3d::operator(): Third index " ,i3,this->n[3]);
      return this->elem[i1-1 + this->n[0]*((i2-1) + this->n[1]*((i3-1) + this->n[2]*(i4-1)))];
    }

    //! \brief Add the given matrix \b X to \a *this.
    matrix4d<T>& operator+=(const matrix4d<T>& X) { return this->add(X); }
    //! \brief Subtract the given matrix \b X from \a *this.
    matrix4d<T>& operator-=(const matrix4d<T>& X) { return this->add(X,T(-1)); }
    //! \brief Add the given matrix \b X scaled by \a alfa to \a *this.
    matrix4d<T>& add(const matrix4d<T>& X, T alfa = T(1))
    {
      return static_cast<matrix4d<T>&>(this->matrixBase<T>::add(X,alfa));
    }

    //! \brief Multiplication of this matrix by a scalar \a c.
    matrix4d<T>& multiply(T c)
    {
      return static_cast<matrix4d<T>&>(this->matrixBase<T>::multiply(c));
    }

    //! \brief Extract a column from the matrix.
    vector<T> getColumn(size_t i2, size_t i3, size_t i4) const
    {
      CHECK_INDEX("matrix3d::getColumn: Second index ",i2,this->n[1]);
      CHECK_INDEX("matrix3d::getColumn: Third index " ,i3,this->n[2]);
      CHECK_INDEX("matrix3d::getColumn: Fourth index ",i4,this->n[3]);
      if (this->n[1] < 2 && this->n[2] < 2 && this->n[3] < 2) return this->elem;
      vector<T> col(this->n[0]);
      memcpy(col.ptr(),this->ptr(i2-1+this->n[1]*((i3-1) + (i4-1)*this->n[2])),col.size()*sizeof(T));
      return col;
    }


    //! \brief Fill a column of the matrix.
    void fillColumn(size_t i2, size_t i3, size_t i4, const std::vector<T>& data)
    {
      CHECK_INDEX("matrix3d::fillColumn: Second index ",i2,this->n[1]);
      CHECK_INDEX("matrix3d::fillColumn: Third index " ,i3,this->n[2]);
      CHECK_INDEX("matrix3d::fillColumn: Third index " ,i4,this->n[3]);
      size_t ndata = this->n[0] > data.size() ? data.size() : this->n[0];
      memcpy(this->ptr(i2-1+this->n[1]*((i3-1))+(i4-1)*this->n[2]),&data.front(),ndata*sizeof(T));
    }

  protected:
    //! \brief Clears the content if any of the first two dimensions changed.
    virtual void clearIfNrowChanged(size_t n1, size_t n2, size_t n3)
    {
      if (n1 != this->n[0] || n2 != this->n[1] || n3 != this->n[2])
        this->elem.clear();
    }
  };


#ifdef HAS_BLAS
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

  template<> inline
  bool matrix3d<float>::multiplyMat(const std::vector<float>& A,
                                    const matrix3d<float>& B, bool addTo)
  {
    size_t M, N, K;
    if (!this->compatible(A,B,M,N,K)) return false;
    if (!addTo) this->resize(M,B.n[1],B.n[2]);

    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                M, N, K, 1.0f,
                A.data(), M,
                B.ptr(), B.n[0],
                addTo ? 1.0f : 0.0f,
                this->ptr(), n[0]);

    return true;
  }

  template<> inline
  bool matrix3d<double>::multiplyMat(const std::vector<double>& A,
                                     const matrix3d<double>& B, bool addTo)
  {
    size_t M, N, K;
    if (!this->compatible(A,B,M,N,K)) return false;
    if (!addTo) this->resize(M,B.n[1],B.n[2]);

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                M, N, K, 1.0,
                A.data(), M,
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
  bool matrix3d<T>::multiply(const matrix<T>& A, const matrix3d<T>& B,
                             bool transA, bool addTo)
  {
    size_t M, N, K;
    if (!this->compatible(A,B,transA,M,N,K)) return false;
    if (!addTo) this->resize(M,B.n[1],B.n[2],true);

    for (size_t i = 1; i <= this->n[0]; i++)
      for (size_t j = 1; j <= this->n[1]; j++)
        for (size_t k = 1; k <= K; k++)
          for (size_t l = 1; l <= this->n[2]; l++)
            if (transA)
              this->operator()(i,j,l) += A(k,i)*B(k,j,l);
            else
              this->operator()(i,j,l) += A(i,k)*B(k,j,l);

    return true;
  }

  template<class T> inline
  bool matrix3d<T>::multiplyMat(const std::vector<T>& A, const matrix3d<T>& B,
                                bool addTo)
  {
    size_t M, N, K;
    if (!this->compatible(A,B,M,N,K)) return false;
    if (!addTo) this->resize(M,B.n[1],B.n[2],true);

    for (size_t i = 1; i <= this->n[0]; i++)
      for (size_t j = 1; j <= this->n[1]; j++)
        for (size_t k = 1; k <= K; k++)
          for (size_t l = 1; l <= this->n[2]; l++)
            this->operator()(i,j,l) += A[i-1+M*(k-1)]*B(k,j,l);

    return true;
  }
#endif

  //! \brief Print the 3D matrix \b A to the stream \a s.
  //! \details The matrix is printed as a set of 2D sub-matrices
  //! based on the first two indices.
  template<class T> std::ostream& operator<<(std::ostream& s,
                                             const matrix3d<T>& A)
  {
    if (A.empty())
      return s <<" (empty)"<< std::endl;

    s <<"Dimension: "<< A.dim(1) <<" "<< A.dim(2) <<" "<< A.dim(3);
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
}
#endif
