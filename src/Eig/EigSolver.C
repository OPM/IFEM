// $Id$
//==============================================================================
//!
//! \file EigSolver.C
//!
//! \date Apr 20 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Interface to LAPack and ARPack eigenvalue solvers.
//!
//==============================================================================

#include "EigSolver.h"
#include "DenseMatrix.h"
#include "SPRMatrix.h"
#ifdef HAS_SLEPC
#include "PETScMatrix.h"
#endif

#if defined(_WIN32) && !defined(__MINGW32__) && !defined(__MINGW64__)
#define eig_av_   EIG_AV
#define eig_mv_   EIG_MV
#define eig_sol_  EIG_SOL
#define eig_drv1_ EIG_DRV1
#define eig_drv2_ EIG_DRV2
#define eig_drv3_ EIG_DRV3
#define eig_drv4_ EIG_DRV4
#define eig_drv5_ EIG_DRV5
#define eig_drv6_ EIG_DRV6
#elif defined(_AIX)
#define eig_av_   eig_av
#define eig_mv_   eig_mv
#define eig_sol_  eig_sol
#define eig_drv1_ eig_drv1
#define eig_drv2_ eig_drv2
#define eig_drv3_ eig_drv3
#define eig_drv4_ eig_drv4
#define eig_drv5_ eig_drv5
#define eig_drv6_ eig_drv6
#endif

extern "C" {
//! \brief Driver to solve a standard eigenvalue problem.
//! \details The regular mode is used.
//! This is a FORTRAN 77 subroutine based on the example driver DSDRV1 from
//! the ARPACK distribution.
//! \sa ARPACK documentation.
void eig_drv1_(const int& n, const int& nev, const int& ncv,
	       double* d, double* v, double* work, int& ierr);
//! \brief Driver to solve a standard eigenvalue problem.
//! \details The shift-and-invert mode is used.
//! This is a FORTRAN 77 subroutine based on the example driver DSDRV2 from
//! the ARPACK distribution.
//! \sa ARPACK documentation.
void eig_drv2_(const int& n, const int& nev, const int& ncv, const double& sig,
	       double* d, double* v, double* work, int& ierr);
//! \brief Driver to solve a generalized eigenvalue problem.
//! \details The inverse mode is used.
//! This is a FORTRAN 77 subroutine based on the example driver DSDRV3 from
//! the ARPACK distribution.
//! \sa ARPACK documentation.
void eig_drv3_(const int& n, const int& nev, const int& ncv,
	       double* d, double* v, double* work, int& ierr);
//! \brief Driver to solve a generalized eigenvalue problem.
//! \details The shift-and-invert mode is used.
//! This is a FORTRAN 77 subroutine based on the example driver DSDRV4 from
//! the ARPACK distribution.
//! \sa ARPACK documentation.
void eig_drv4_(const int& n, const int& nev, const int& ncv, const double& sig,
	       double* d, double* v, double* work, int& ierr);
//! \brief Driver to solve a generalized eigenvalue problem.
//! \details The buckling mode is used.
//! This is a FORTRAN 77 subroutine based on the example driver DSDRV5 from
//! the ARPACK distribution.
//! \sa ARPACK documentation.
void eig_drv5_(const int& n, const int& nev, const int& ncv, const double& sig,
	       double* d, double* v, double* work, int& ierr);
//! \brief Driver to solve a generalized eigenvalue problem.
//! \details The Cayley mode is used.
//! This is a FORTRAN 77 subroutine based on the example driver DSDRV6 from
//! the ARPACK distribution.
//! \sa ARPACK documentation.
void eig_drv6_(const int& n, const int& nev, const int& ncv, const double& sig,
	       double* d, double* v, double* work, int& ierr);
}

static SystemMatrix* K  = 0; //!< Pointer to coefficient matrix A
static SystemMatrix* M  = 0; //!< Pointer to coefficient matrix B
static SystemMatrix* AM = 0; //!< Pointer to the matrix to invert


bool eig::solve (SystemMatrix* A, SystemMatrix* B,
		 Vector& eigVal, Matrix& eigVec, int nev)
{
  DenseMatrix* dK = dynamic_cast<DenseMatrix*>(A);
  DenseMatrix* dM = dynamic_cast<DenseMatrix*>(B);

  if (dK)
    if (dM) // Try the LAPack solver DSYGVX
      return dK->solveEig(*dM,eigVal,eigVec,nev);
    else    // Try the LAPack solver DSYEVX
      return dK->solveEig(eigVal,eigVec,nev);

  SPRMatrix* sK = dynamic_cast<SPRMatrix*>(A);
  SPRMatrix* sM = dynamic_cast<SPRMatrix*>(B);

  if (sK && sM) // Try the Lanczos solver SPRLAN
    return sK->solveEig(*sM,eigVal,eigVec,nev);

#ifdef HAS_SLEPC
  PETScMatrix *pK = dynamic_cast<PETScMatrix*>(A);
  PETScMatrix *pM = dynamic_cast<PETScMatrix*>(B);

  if (pK && pM) // Try the SLEPc solver EPSSolve
    return pK->solveEig(*pM,eigVal,eigVec,nev);
#endif

  return false;
}


bool eig::solve (SystemMatrix* A, SystemMatrix* B,
		 Vector& eigVal, Matrix& eigVec, int nev, int ncv,
		 int mode, double shift)
{
  K = A;
  M = B;
  int ierr = 1;
  int n = K->dim();
  int nwork = 4*n+ncv*(ncv+9);
  if (mode == 6) nwork += n;
  eigVal.resize(2*ncv);
  eigVec.resize(n,ncv);
  double* work = new double[nwork];

  switch (mode) {
  case 1:
    eig_drv1_(n,nev,ncv,eigVal.ptr(),eigVec.ptr(),work,ierr);
    break;
  case 2:
    AM = K->copy();
    if (shift != 0.0 && !AM->add(-shift))
      ierr = 123;
    else
      eig_drv2_(n,nev,ncv,shift,eigVal.ptr(),eigVec.ptr(),work,ierr);
    break;
  case 3:
    AM = M->copy();
    eig_drv3_(n,nev,ncv,eigVal.ptr(),eigVec.ptr(),work,ierr);
    break;
  case 4:
    AM = K->copy();
    if (shift != 0.0 && !AM->add(*M,-shift))
      ierr = 123;
    else
      eig_drv4_(n,nev,ncv,shift,eigVal.ptr(),eigVec.ptr(),work,ierr);
    break;
  case 5:
    AM = K->copy();
    if (shift != 0.0 && !AM->add(*M,+shift)) // Notice the +sign on the shift!
      ierr = 123;
    else
      eig_drv5_(n,nev,ncv,shift,eigVal.ptr(),eigVec.ptr(),work,ierr);
    break;
  case 6:
    AM = K->copy();
    if (shift != 0.0 && AM->add(*M,-shift))
      ierr = 123;
    else
      eig_drv6_(n,nev,ncv,shift,eigVal.ptr(),eigVec.ptr(),work,ierr);
    break;
  }

  if (ierr == 123)
    std::cerr <<" *** eig::solve: Failed to add system matrices.\n"
	      <<"                 Check matrix type or dimensions."<< std::endl;

  delete[] work;
  if (AM) delete AM;
  AM = 0;
  K = 0;
  M = 0;
  return ierr == 0;
}


//! \brief Performs the matrix-vector multiplication \b y = \b K * \b x.
extern "C" void eig_av_(const double* x, double* y)
{
  if (!K) return;

  StdVector Y;
  K->multiply(StdVector(x,K->dim()),Y);
  memcpy(y,Y.ptr(),Y.dim()*sizeof(double));
}


//! \brief Performs the matrix-vector multiplication \b y = \b M * \b x.
extern "C" void eig_mv_(const double* x, double* y)
{
  if (!M) return;

  StdVector Y;
  M->multiply(StdVector(x,M->dim()),Y);
  memcpy(y,Y.ptr(),Y.dim()*sizeof(double));
}


//! \brief Solves the linear system of equations \b AM \b y = \b x.
//! \details The matrix \b AM depends on the eigensolver method.
extern "C" void eig_sol_(const double* x, double* y, int& ierr)
{
  ierr = -99;
  if (!AM) return;

  StdVector RHS(x,AM->dim());
  ierr = AM->solve(RHS) ? 0 : -1;
  memcpy(y,RHS.ptr(),RHS.dim()*sizeof(double));
}
