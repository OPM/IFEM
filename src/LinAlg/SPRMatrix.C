// $Id$
//==============================================================================
//!
//! \file SPRMatrix.C
//!
//! \date Jan 4 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of the system matrix on the SPR-format.
//!
//==============================================================================

#include "SPRMatrix.h"
#include "SAM.h"
#include <numeric>
#include <fstream>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#if defined(_WIN32) && !defined(__MINGW32__) && !defined(__MINGW64__)
#define sprprp_ SPRPRP
#define sprsas_ SPRSAS
#define sprrnm_ SPRRNM
#define sprtrs_ SPRTRS
#define sprsmb_ SPRSMB
#define sprpmp_ SPRPMP
#define spradm_ SPRADM
#define sprdad_ SPRDAD
#define sprprm_ SPRPRM
#define sprsol_ SPRSOL
#define sprlax_ SPRLAX
#define sprprn_ SPRPRN
#define sprcnv_ SPRCNV
#define openftnfile_ OPENFTNFILE
#define closeftnfile_ CLOSEFTNFILE
#elif defined(_AIX)
#define sprprp_ sprprp
#define sprsas_ sprsas
#define sprrnm_ sprrnm
#define sprtrs_ sprtrs
#define sprsmb_ sprsmb
#define sprpmp_ sprpmp
#define spradm_ spradm
#define sprdad_ sprdad
#define sprprm_ sprprm
#define sprsol_ sprsol
#define sprlax_ sprlax
#define sprprn_ sprprn
#define sprcnv_ sprcnv
#define openftnfile_ openftnfile
#define closeftnfile_ closeftnfile
#endif

extern "C" {
//! \brief Prepares the control information for the sparse assembly process.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprprp_(const Int_* madof, const Int_* minex,  const Int_* mpmnpc,
             const Int_* mmnpc, const Int_* mpmceq, const Int_* mmceq,
             const Int_* msc,   const Int_& nspar,  const Int_& lpu,
             const Int_* mpar,  Int_* mspar,
             const Int_* meqn,  Int_* iWork, Int_& ierr);
//! \brief Computes the SPR-version of the connectivity array \a mmnpc.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprsas_(const Int_* mpar,  const Int_* mpmnpc, const Int_* mmnpc,
             const Int_* madof, const Int_* msc,    const Int_* mpmceq,
             const Int_* mmceq, const Int_* meqn,
             Int_* mspar, Int_* msica, Int_* iWork,
             const Int_& nspar, const Int_& lpu, Int_& ierr);
//! \brief Reorders the equations by means of the SPR-node partition.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprrnm_(Int_* mspar, Int_* msica, Int_* iWork,
             const Int_& lpu, Int_& ierr);
//! \brief Analyses the sparsity pattern of the system matrix.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprtrs_(Int_* mspar, Int_* msica, Int_* mtrees, const Int_* meqn,
             const Int_& nspar, const Int_& nmsica, const Int_& ntrees,
             const Int_& ndof, const Int_& niwork, Int_* iWork, Real* rinfo,
             const Int_& lpu, Int_& ierr);
//! \brief Performs the symbolic factorization of the system matrix.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprsmb_(Int_* mspar, Int_* msica, Int_* mtrees, Int_* msifa,
             const Int_* meqn, Int_* iWork, const Int_& lpu, Int_& ierr);
//! \brief Finalizes the SPR control arrays \a msica, \a mtrees and \a msifa.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprpmp_(const Int_* msica, const Int_* mtrees, const Int_* meqn,
             const Int_* mpar, Int_* mspar, Int_* mvarnc);
//! \brief Assembles an element matrix \a EM into the system matrix \a SM.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void spradm_(const Real* eK, const Real* ttcc, const Int_* mpar,
             const Int_* mspar, const Int_* madof, const Int_* meqn,
             const Int_* mpmnpc, const Int_* mmnpc, const Int_* mpmceq,
             const Int_* mmceq, const Int_* msica, const Int_* mtrees,
             const Int_* msifa, const Int_* mvarnc, Real* values, Real* sysRHS,
             Int_* work, const Int_& iel, const Int_& nedof, const Int_& lpu,
             const Int_& nrhs, Int_& ierr);
//! \brief Adds a scalar value into the diagonal of the system matrix.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprdad_(const Int_* mpar, const Int_* mtrees, const Int_* msifa,
             Real* values, const Real& sigma, const Int_& lpu, Int_& ierr);
//! \brief Performs a matrix-matrix multiplication.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprprm_(const Real* A, const Real* B, Real* C, Real* rWork,
             const Int_* mspar, const Int_* mtrees, const Int_* msifa,
             const Int_& m, const Int_& n, const Int_& ksa, const Int_& iflag,
             const Int_& lpu, Int_& ierr);
//! \brief Solves the linear equation system \a Ax=b.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprsol_(const Int_& iop, const Int_* mspar, const Int_* mtrees,
             const Int_* msifa, Real* value, Real* B,
             const Int_& ldB, const Int_& nrhs, Real* tol,
             Int_* iWork, Real* rWork, const Int_& lpu, Int_& ierr);
//! \brief Solves the generalized eigenproblem \a Ax=&lambda;Bx.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprlax_(Real* A, Real* B, const Real* tol,
             const Int_* mparA, const Int_* mtreeA, const Int_* msifA,
             const Int_* mparB, const Int_* mtreeB, const Int_* msifB,
             const Int_& iop, Real* val, Real* vec, Real* rWork, Int_* iWork,
             const Real& shift, const Int_& n, const Int_& nv, const Int_& maxl,
             const Int_& lpu, const Int_& ipsw, Int_& ierr);
//! \brief Prints the matrix content to the specified Fortran unit number.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprprn_(const Real* L,
             const Int_* mspar, const Int_* mtrees, const Int_* msifa,
             const Int_& lpu, Int_& ierr);
//! \brief Converts the matrix to a full matrix.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprcnv_(const Real* L, Real* A,
             const Int_* mspar, const Int_* mtrees, const Int_* msifa,
             const Int_& m, const Int_& ksa, const Int_& lpu, Int_& ierr);
//! \brief Opens a temporary file for Fortran print.
int openftnfile_(const char* fname, const ssize_t nchar);
//! \brief Closes a Fortran file.
void closeftnfile_(const int& iunit);
}

using Int_Vec = std::vector<Int_>; //!< Convenience type alias


/*!
  \brief Generic function for copying a C-array.
*/

template<class T> T* copyArr(const T* array, int n)
{
  T* newarr = new T[n];
  memcpy(newarr,array,n*sizeof(T));
  return newarr;
}


SPRMatrix::SPRMatrix () : rWork(1)
{
  ierr = 0;
  memset(mpar,0,sizeof(mpar));
  msica = msifa = mtrees = mvarnc = nullptr;
  values = nullptr;

#ifdef USE_OPENMP
  if (int nt = omp_get_max_threads(); nt > 0)
    jWork = new Int_Vec[nt];
#else
  jWork = &iWork;
#endif
  mySam = nullptr;
}


SPRMatrix::SPRMatrix (const SPRMatrix& A) : SystemMatrix(A), rWork(1)
{
  ierr   = 0;
  memcpy(mpar,A.mpar,sizeof(A.mpar));
  msica  = copyArr(A.msica,A.mpar[1]);
  msifa  = copyArr(A.msifa,A.mpar[2]);
  mtrees = copyArr(A.mtrees,A.mpar[35]);
  mvarnc = copyArr(A.mvarnc,A.mpar[7]*2);
  values = copyArr(A.values,A.mpar[7]+A.mpar[15]);

#ifdef USE_OPENMP
  if (int nt = omp_get_max_threads(); nt > 0)
    jWork = new Int_Vec[nt];
#else
  jWork = &iWork;
#endif
  mySam = A.mySam;
#ifdef USE_INT64
  sharedSam = true;
#endif
}


SPRMatrix::~SPRMatrix ()
{
  delete[] msica;
  delete[] msifa;
  delete[] mtrees;
  delete[] mvarnc;
  delete[] values;
#ifdef USE_OPENMP
  delete[] jWork;
#endif
#ifdef USE_INT64
  if (!sharedSam)
    delete mySam;
#endif
}


#ifdef USE_INT64
SPRMatrix::SAM64::SAM64 (const SAM& sam)
{
  // Lambda function making a 64-bit version of a (32-bit) int array.
  auto&& copyTo64 = [](const int* array, int n)
  {
    Int_* array64 = new Int_[n];
    std::copy(array,array+n,array64);
    return array64;
  };

  mpar   = copyTo64(sam.mpar,sizeof(sam.mpar)/sizeof(int));
  mpmnpc = copyTo64(sam.mpmnpc,sam.nel+1);
  mmnpc  = copyTo64(sam.mmnpc,sam.nmmnpc);
  madof  = copyTo64(sam.madof,sam.nnod+1);
  mpmceq = copyTo64(sam.mpmceq,sam.nceq+1);
  mmceq  = copyTo64(sam.mmceq,sam.nmmceq);
  meqn   = copyTo64(sam.meqn,sam.ndof);
  minex  = sam.minex;
}


SPRMatrix::SAM64::~SAM64 ()
{
  delete[] mpar;
  delete[] mpmnpc;
  delete[] mmnpc;
  delete[] madof;
  delete[] mpmceq;
  delete[] mmceq;
  delete[] meqn;
}


std::pair<int,int> SPRMatrix::SAM64::getNodeAndLocalDof (int ieq, bool) const
{
  const Int_ nnod = mpar[0];
  const Int_ ndof = mpar[2];
  if (Int_* idof = std::find(meqn,meqn+ndof,ieq); idof != meqn+ndof)
    ieq = std::distance(meqn,idof);
  else
    return { 0, 0 };

  for (int n = 1; n <= nnod; n++)
    if (madof[n] > static_cast<Int_>(ieq))
      return { minex ? minex[n-1] : n, ieq-madof[n-1]+1 };

  return { 0, 0 };
}
#endif


size_t SPRMatrix::dim (int idim) const
{
  switch (idim) {
  case 1:
  case 2:
    return mpar[7];
  case 3:
    return mpar[7] * mpar[7];
  default:
    return mpar[7] + mpar[15];
  }
}


/*!
  Must be called once before the element assembly loop.
  The SPR data structures are initialized and the all symbolic operations
  that are need before the actual assembly can start are performed.
*/

void SPRMatrix::initAssembly (const SAM& sam, char)
{
  memset(mpar,0,sizeof(mpar));
  msica = msifa = mtrees = mvarnc = nullptr;
  values = nullptr;

#ifdef HAS_SPR
  if (!sam.madof || !sam.mpmnpc)
  {
    std::cerr <<"SPRMatrix: SAM object is not properly initialized"<< std::endl;
    return;
  }

  iWork.resize(2*sam.nnod + 2*sam.ndof);
#ifdef USE_INT64
  Int_* minex = nullptr;
#else
  int* minex = sam.minex;
#endif
  if (!minex)
  {
    minex = iWork.data();
    if (sam.minex)
      std::copy(sam.minex,sam.minex+sam.nnod,minex);
    else // External node numbers are not provided - generate an identity array
      std::iota(minex,minex+sam.nnod,1);
  }

#ifdef USE_INT64
  // Make 64-bit versions of the SAM arrays
  mySam = new SAM64(sam);
  Int_* msc = new Int_[sam.ndof];
  std::copy(sam.msc,sam.msc+sam.ndof,msc);
#else
  mySam = const_cast<SAM*>(&sam);
  int* msc = sam.msc;
#endif
  for (int idof = 0; idof < sam.ndof; idof++)
    if (msc[idof] < 0)
    {
      // Make a temporary copy of msc with non-negative values only
#ifndef USE_INT64
      if (msc == sam.msc)
      {
        msc = new Int_[sam.ndof];
        std::copy(sam.msc,sam.msc+sam.ndof,msc);
      }
#endif
      msc[idof] = 0;
    }

  bool finished = true;
  sprprp_(mySam->madof, minex, mySam->mpmnpc, mySam->mmnpc,
          mySam->mpmceq, mySam->mmceq, msc, NS, 6, mySam->mpar,
          mpar, mySam->meqn, iWork.data()+sam.nnod, ierr);
  if (ierr < 0)
    std::cerr <<"SAM::SPRPRM: Failure "<< ierr << std::endl;
  else if (sam.neq != mpar[7])
    std::cerr <<"SAM::SPRPRM: Internal error, NEQ = "<< sam.neq
              <<" != "<< mpar[7] << std::endl;
  else if (mpar[7] > 0)
    finished = false;

#ifdef USE_INT64
  std::copy(mySam->meqn,mySam->meqn+sam.ndof,sam.meqn);
#else
  if (msc != sam.msc) delete[] msc;
#endif
  if (finished) return;

  // Restore the original msc
#ifdef USE_INT64
  std::copy(sam.msc,sam.msc+sam.ndof,msc);
#else
  msc = sam.msc;
#endif

  // Perform symbolic assembly
  Int_Vec itemp(mpar[34]);
  sprsas_(mySam->mpar, mySam->mpmnpc, mySam->mmnpc, mySam->madof, msc,
          mySam->mpmceq, mySam->mmceq, mySam->meqn, mpar,
          itemp.data(), iWork.data(), NS, 6, ierr);
  if (ierr < 0)
  {
    std::cerr <<"SAM::SPRSAS: Failure "<< ierr << std::endl;
    return;
  }

  // Reallocation of msica with the correct size
  msica = new Int_[mpar[1]];
  memcpy(msica,itemp.data(),mpar[1]*sizeof(Int_));

  // Perform nodal reordering
  iWork.resize(mpar[36]);
  sprrnm_(mpar, msica, iWork.data(), 6, ierr);
  if (ierr < 0)
  {
    std::cerr <<"SAM::SPRRNM: Failure "<< ierr << std::endl;
    return;
  }

  // Analyze the sparsity pattern
  mtrees = new Int_[mpar[35]];
  iWork.resize(mpar[37]);
  Real rinfo[4];
  sprtrs_(mpar, msica, mtrees, mySam->meqn, NS, mpar[1], mpar[35], sam.ndof,
          mpar[37], iWork.data(), rinfo, 6, ierr);
  if (ierr < 0)
  {
    std::cerr <<"SAM::SPRTRS: Failure "<< ierr << std::endl;
    return;
  }
#ifdef USE_INT64
  std::copy(mySam->meqn,mySam->meqn+sam.ndof,sam.meqn);
#endif

  // Perform symbolic factorization
  itemp.resize(mpar[2]);
  iWork.resize(3*mpar[5] + 2*mpar[7] + 1);
  sprsmb_(mpar, msica, mtrees, itemp.data(), mySam->meqn,
          iWork.data(), 6, ierr);
  if (ierr < 0)
  {
    std::cerr <<"SAM::SPRSMB: Failure "<< ierr << std::endl;
    return;
  }

  // Reallocation of msifa with the correct size
  msifa = new Int_[mpar[2]];
  memcpy(msifa,itemp.data(),mpar[2]*sizeof(Int_));

  // Finalize the SPR datastructure
  mvarnc = new Int_[2*mpar[7]];
  sprpmp_(msica, mtrees, mySam->meqn, mySam->mpar, mpar, mvarnc);

  // Allocate space for the matrix itself
  values = new Real[mpar[7] + mpar[15]];
  memset(values,0,(mpar[7] + mpar[15])*sizeof(Real));
#else
  std::cerr <<"SPRMatrix: SPR solver not available"<< std::endl;
#endif
}


void SPRMatrix::init ()
{
  if (mpar[0] < 4) return;

  mpar[0] = 4;
  memset(values,0,(mpar[7] + mpar[15])*sizeof(Real));
}


void SPRMatrix::dump (std::ostream& os, LinAlg::StorageFormat fmt, const char*)
{
  if (fmt == LinAlg::BINARY)
  {
    os.write("#binary:",8);
    os.write((const char*)mpar,sizeof(mpar));
    os.write((const char*)msica,mpar[1]*sizeof(Int_));
    os.write((const char*)msifa,mpar[2]*sizeof(Int_));
    os.write((const char*)mtrees,mpar[35]*sizeof(Int_));
    os.write((const char*)mvarnc,2*mpar[7]*sizeof(Int_));
    os.write((const char*)values,(mpar[7]+mpar[15])*sizeof(Real));
  }
  else
  {
    os << mpar[0];
    for (int i = 1; i < NS; i++)
      os <<' '<< mpar[i];

    if (msica && mpar[1] > 0)
    {
      os <<'\n'<< msica[0];
      for (int i = 1; i < mpar[1]; i++)
        os <<' '<< msica[i];
    }

    if (msifa && mpar[2] > 0)
    {
      os <<'\n'<< msifa[0];
      for (int i = 1; i < mpar[2]; i++)
        os <<' '<< msifa[i];
    }

    if (mtrees && mpar[35] > 0)
    {
      os <<'\n'<< mtrees[0];
      for (int i = 1; i < mpar[35]; i++)
        os <<' '<< mtrees[i];
    }

    if (mvarnc && mpar[7] > 0)
    {
      os <<'\n'<< mvarnc[0];
      for (int i = 1; i < 2*mpar[7]; i++)
        os <<' '<< mvarnc[i];
    }

    if (values && mpar[7]+mpar[15] > 0)
      for (int i = 0; i < mpar[7]+mpar[15]; i++)
        os << (i%100 ? ' ' : '\n') << values[i];

    os << std::endl;
  }
}


bool SPRMatrix::load (const char* fileName, bool binary)
{
  std::ifstream is(fileName,
                   binary ? std::ifstream::binary : std::ifstream::in);
  if (!is)
  {
    std::cerr <<" *** SPRMatrix::load: Failed to read file "<< fileName
              << std::endl;
    return false;
  }

  if (binary)
  {
    if (char header[8]; !is.read(header,8) || strncmp(header,"#binary:",8))
    {
      std::cerr <<" *** SPRMatrix::load: Invalid binary matrix file "
                << fileName <<" ("<< std::string(header,8)
                <<")."<< std::endl;
      return false;
    }
    is.read((char*)mpar,sizeof(mpar));
    if (msica)
      is.read((char*)msica,mpar[1]*sizeof(Int_));
    if (msifa)
      is.read((char*)msifa,mpar[2]*sizeof(Int_));
    if (mtrees)
      is.read((char*)mtrees,mpar[35]*sizeof(Int_));
    if (mvarnc)
      is.read((char*)mvarnc,2*mpar[7]*sizeof(Int_));
    if (values)
      is.read((char*)values,(mpar[7]+mpar[15])*sizeof(Real));
  }
  else
  {
    for (int i = 0; i < NS && is; i++)
      is >> mpar[i];

    if (msica)
      for (int i = 0; i < mpar[1] && is; i++)
        is >> msica[i];

    if (msifa)
      for (int i = 0; i < mpar[2] && is; i++)
        is >> msifa[i];

    if (mtrees)
      for (int i = 0; i < mpar[35] && is; i++)
        is >> mtrees[i];

    if (mvarnc)
      for (int i = 0; i < 2*mpar[7] && is; i++)
        is >> mvarnc[i];

    if (values)
      for (int i = 0; i < mpar[7]+mpar[15] && is; i++)
        is >> values[i];
  }

  return is.good() && this->flagNonZeroEqs();
}


bool SPRMatrix::assemble (const Matrix& eM, const SAM& sam, int e)
{
  return this->assemble(e,eM,sam);
}


bool SPRMatrix::assemble (const Matrix& eM, const SAM& sam,
                          SystemVector& B, int e)
{
  if (B.getType() != LinAlg::DENSE) return false;

  return this->assemble(e,eM,sam,B.getPtr());
}


bool SPRMatrix::assemble (const Matrix&, const SAM&,
                          SystemVector&, const IntVec&)
{
  std::cerr <<"SPRMatrix::assemble(const Matrix&,const SAM&,"
            <<"SystemVector&,const IntVec&): Not implemented."
            << std::endl;
  return false;
}


bool SPRMatrix::assemble (const Matrix&, const IntVec&)
{
  std::cerr <<"SPRMatrix::assemble(const Matrix&,const IntVec&):"
            <<" Not implemented."<< std::endl;
  return false;
}


bool SPRMatrix::assemble (int e, const Matrix& eM, const SAM& sam, Real* B)
{
#ifdef HAS_SPR
#ifdef USE_OPENMP
  Int_Vec& IWORK = jWork[omp_in_parallel() ? omp_get_thread_num() : 0];
#else
  Int_Vec& IWORK = *jWork;
#endif
  IWORK.resize(mpar[38]);
  spradm_(eM.ptr(), sam.ttcc, mySam->mpar, mpar, mySam->madof, mySam->meqn,
          mySam->mpmnpc, mySam->mmnpc, mySam->mpmceq, mySam->mmceq,
          msica, mtrees, msifa, mvarnc, values, B ? B : rWork.data(),
          IWORK.data(), e, eM.rows(), 6, B ? 1 : 0, ierr);
  if (ierr == 0)
    return this->flagNonZeroEqs({IWORK.begin(),IWORK.begin()+eM.rows()});

  std::cerr <<"SAM::SPRADM: Failure "<< ierr << std::endl;
#endif
  return false;
}


void SPRMatrix::mult (Real alpha)
{
  cblas_dscal(mpar[7]+mpar[15], alpha, values, 1);
}


bool SPRMatrix::add (const SystemMatrix& B, Real alpha)
{
  const SPRMatrix* Bptr = dynamic_cast<const SPRMatrix*>(&B);
  if (!Bptr) return false;

  if (mpar[7] != Bptr->mpar[7] || mpar[15] != Bptr->mpar[15]) return false;

  if (B.isZero()) return true;

  cblas_daxpy(mpar[7]+mpar[15], alpha, Bptr->values, 1, values, 1);

  return this->flagNonZeroEqs(B);
}


bool SPRMatrix::add (Real sigma, int ieq)
{
#ifdef HAS_SPR
  ierr = ieq > 0 ? ieq : 0;
  sprdad_(mpar, mtrees, msifa, values, sigma, 6, ierr);
  if (ierr < 0)
    std::cerr <<"SAM::SPRDAD: Failure "<< ierr << std::endl;
  else if (ieq > 0)
    return this->flagNonZeroEqs({ieq});
  else
    return this->flagNonZeroEqs();
#endif

  return false;
}


bool SPRMatrix::multiply (size_t n, const Real* b, Real* c)
{
#ifdef HAS_SPR
  if (rWork.size() < n) rWork.resize(n);
  sprprm_(values+n, b, c, rWork.data(), mpar, mtrees, msifa,
          n, 1, 1, 1, 6, ierr);
  if (!ierr) return true;

  std::cerr <<"SAM::SPRPRM: Failure "<< ierr << std::endl;
#endif
  return false;
}


bool SPRMatrix::multiply (const SystemVector& B, SystemVector& C) const
{
  if (mpar[7] <= 0) return true;

  if (B.getType() != LinAlg::DENSE) return false;
  if (C.getType() != LinAlg::DENSE) return false;

  C.redim(mpar[7]);
  return const_cast<SPRMatrix*>(this)->multiply(C.dim(),B.getRef(),C.getPtr());
}


bool SPRMatrix::solve (SystemVector& B, Real*)
{
  if (mpar[7] < 1) return true; // No equations to solve

  if (B.getType() != LinAlg::DENSE) return false;

#ifdef HAS_SPR
  Real tol[3] = { Real(1.0e-12), Real(0), Real(0) };
  iWork.resize(std::max(mpar[12],mpar[13]+1));
  rWork.resize(std::max(mpar[16],mpar[13]));
  Int_ iop = mpar[0] < 5 ? 3 : 4;
  sprsol_(iop, mpar, mtrees, msifa, values, B.getPtr(),
          B.dim(), 1, tol, iWork.data(), rWork.data(), 6, ierr);
  if (ierr < 0)
    std::cerr <<"SAM::SPRSOL: Failure ("<< ierr <<")";
  else if (ierr == 3)
    std::cerr <<"SAM::SPRSOL: Negative pivot";
  else
    return true;

  std::pair<int,int> dof = mySam->getNodeAndLocalDof(mpar[26],true);
  if (dof.first > 0)
    std::cerr <<" detected in node "<< dof.first <<", local DOF "<< dof.second;
  else if (mpar[26] > 0)
    std::cerr <<", for equation "<< mpar[26] <<" (unknown) ";
  std::cerr << std::endl;
#else
  ierr = -99;
#endif
  return ierr >= 0;
}


/*!
  The eigenproblem is assumed to be on the form
  \b A \b x = &lambda; \b B \b x where \b A ( = \a *this ) and \b B
  both are assumed to be symmetric and \b B also to be positive definite.
  The eigenproblem is solved by the SAM library subroutine \a SPRLAN.
  \sa SAM library documentation.
*/

bool SPRMatrix::solveEig (SPRMatrix& B, RealArray& val, Matrix& vec, int nev,
                          Real shift, int iop)
{
  const Int_ n = mpar[7];
  if (n < 1 || nev == 0) return true; // No equations to solve
  if (n != B.mpar[7]) return false; // Incompatible matrices

#ifdef USE_F77SAM
#ifdef SP_DEBUG
  const int ipsw = SP_DEBUG;
#else
  const int ipsw = 0;
#endif
  const Real tol[3] = { Real(1.0e-8), Real(1.0e-12), Real(1.0e-8) };
  Int_ maxlan = 3*nev+12;
  if (maxlan > n) maxlan = n;
  Int_ nrWork = std::max(maxlan*(maxlan+7), std::max(mpar[13],mpar[16]));
  val.resize(n);
  vec.resize(n,maxlan);
  iWork.resize(std::max(2*maxlan+mpar[12],n));
  rWork.resize(maxlan + 2*n + nrWork);
  sprlax_(values, B.values, tol, mpar, mtrees, msifa, B.mpar, B.mtrees, B.msifa,
          iop, val.data(), vec.ptr(), rWork.data(), iWork.data(),
          shift, n, nev, maxlan, 6, ipsw, ierr);
  val.resize(nev);
  vec.resize(n,nev);
  if (!ierr) return true;

  std::cerr <<"SAM::SPRLAX: Failure "<< ierr << std::endl;
#else
  std::cerr <<"SAM::SPRLAX eigensolver not available"<< std::endl;
#endif
  return false;
}


Real SPRMatrix::Linfnorm () const
{
  RealArray B(mpar[7],Real(1)), C(mpar[7],Real(0));
  if (const_cast<SPRMatrix*>(this)->multiply(B.size(),B.data(),C.data()))
    return *std::max_element(C.begin(),C.end());

  return Real(0);
}


//! \cond DO_NOT_DOCUMENT
namespace SPR { bool printAsFull = false; }
//! \endcond

std::ostream& SPRMatrix::write (std::ostream& os) const
{
  if (SPR::printAsFull)
  {
    Matrix fullMat;
    if (this->convert(fullMat))
      return os << fullMat;
  }

#ifdef HAS_SPR
  // Open a temporary file for the Fortran output
  const char* tmpfile = "/tmp/spr.tmp";
  Int_ lpu = openftnfile_(tmpfile,strlen(tmpfile)), lerr;
  if (lpu <= 0) return os;
  // Write to the temporary file
  const Int_ n = mpar[7];
  sprprn_(values+n, mpar, mtrees, msifa, lpu, lerr);
  closeftnfile_(static_cast<int>(lpu));
  // Copy temporary file to output stream
  char c;
  std::ifstream is(tmpfile);
  while (is.get(c))
    os.put(c);
#endif
  return os;
}


bool SPRMatrix::convert (Matrix& fullMat) const
{
  const Int_ n = mpar[7];
  Int_ lerr = -99;
  fullMat.resize(n,n);
#ifdef HAS_SPR
  sprcnv_(values+n, fullMat.ptr(), mpar, mtrees, msifa, n, 1, 6, lerr);
#endif
  return lerr == 0;
}
