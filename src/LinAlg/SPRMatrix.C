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
#endif

extern "C" {
//! \brief Prepares the control information for the sparse assembly process.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprprp_(const int* madof, const int* minex,  const int* mpmnpc,
	     const int* mmnpc, const int* mpmceq, const int* mmceq,
	     const int* msc,   const int& nspar,  const int& lpu,
	     const int* mpar,  int* mspar,
	     const int* meqn,  int* iWork, int& ierr);
//! \brief Computes the SPR-version of the connectivity array \a mmnpc.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprsas_(const int* mpar,  const int* mpmnpc, const int* mmnpc,
	     const int* madof, const int* msc,    const int* mpmceq,
	     const int* mmceq, const int* meqn,
	     int* mspar, int* msica, int* iWork, const int& lpu, int& ierr);
//! \brief Reorders the equations by means of the SPR-node partition.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprrnm_(int* mspar, int* msica, int* iWork, const int& lpu, int& ierr);
//! \brief Analyses the sparsity pattern of the system matrix.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprtrs_(int* mspar, int* msica, int* mtrees, const int* meqn,
	     const int& ndof, int* iWork, Real* rinfo,
	     const int& lpu, int& ierr);
//! \brief Performs the symbolic factorization of the system matrix.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprsmb_(int* mspar, int* msica, int* mtrees, int* msifa,
	     const int* meqn, int* iWork, const int& lpu, int& ierr);
//! \brief Finalizes the SPR control arrays \a msica, \a mtrees and \a msifa.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprpmp_(const int* msica, const int* mtrees, const int* meqn,
	     const int* mpar, int* mspar, int* mvarnc);
//! \brief Assembles an element matrix \a EM into the system matrix \a SM.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void spradm_(const Real* eK, const Real* ttcc, const int* mpar,
	     const int* mspar, const int* madof, const int* meqn,
	     const int* mpmnpc, const int* mmnpc, const int* mpmceq,
	     const int* mmceq, const int* msica, const int* mtrees,
	     const int* msifa, const int* mvarnc, Real* values, Real* sysRHS,
	     int* work, const int& iel, const int& nedof, const int& lpu,
	     const int& nrhs, int& ierr);
//! \brief Adds a scalar value into the diagonal of the system matrix.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprdad_(const int* mpar, const int* mtrees, const int* msifa,
	     Real* values, const Real& sigma, const int& lpu, int& ierr);
//! \brief Performs a matrix-matrix multiplication.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprprm_(const Real* A, const Real* B, Real* C, Real* rWork,
	     const int* mspar, const int* mtrees, const int* msifa,
	     const int& m, const int& n, const int& ksa, const int& iflag,
	     const int& lpu, int& ierr);
//! \brief Solves the linear equation system \a Ax=b.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprsol_(const int& iop, const int* mspar, const int* mtrees,
	     const int* msifa, Real* value, Real* B,
	     const int& ldB, const int& nrhs, Real* tol,
	     int* iWork, Real* rWork, const int& lpu, int& ierr);
//! \brief Solves the generalized eigenproblem \a A*x=(lambda)*B*x.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprlax_(Real* A, Real* B, const Real* tol,
	     const int* mparA, const int* mtreeA, const int* msifA,
	     const int* mparB, const int* mtreeB, const int* msifB,
	     const int& iop, Real* val, Real* vec, Real* rWork, int* iWork,
	     const Real& shift, const int& n, const int& nv, const int& maxlan,
	     const int& lpu, const int& ipsw, int& ierr);
}


SPRMatrix::SPRMatrix (const SPRMatrix& A)
{
  msica  = new int[A.mpar[1]];
  msifa  = new int[A.mpar[2]];
  mtrees = new int[A.mpar[35]];
  mvarnc = new int[2*A.mpar[7]];
  values = new Real[A.mpar[7] + A.mpar[15]];

  memcpy(mpar  ,A.mpar  ,NS*sizeof(int));
  memcpy(msica ,A.msica ,A.mpar[1]*sizeof(int));
  memcpy(msifa ,A.msifa ,A.mpar[2]*sizeof(int));
  memcpy(mtrees,A.mtrees,A.mpar[35]*sizeof(int));
  memcpy(mvarnc,A.mvarnc,2*A.mpar[7]*sizeof(int));
  memcpy(values,A.values,(A.mpar[7]+A.mpar[15])*sizeof(Real));
  iWork.clear();
  rWork.clear();
}


SPRMatrix::~SPRMatrix ()
{
  if (msica)  delete[] msica;
  if (msifa)  delete[] msifa;
  if (mtrees) delete[] mtrees;
  if (mvarnc) delete[] mvarnc;
  if (values) delete[] values;
}


void SPRMatrix::initAssembly (const SAM& sam, bool)
{
  memset(mpar,0,NS*sizeof(int));
  msica = 0;
  msifa = 0;
  mtrees = 0;
  mvarnc = 0;
  values = 0;

#ifdef HAS_SPR
  if (!sam.madof || !sam.mpmnpc)
  {
    std::cerr <<"SPRMatrix: SAM object is not properly initialized"<< std::endl;
    return;
  }

  int ierr = 0;
  iWork.resize(sam.nnod + 2*sam.ndof);
  int* minex = sam.minex;
  if (!minex)
  {
    // External node numbers are not provided - generate an identity array
    minex = &iWork.front();
    for (int n = 0; n < sam.nnod; n++)
      minex[n] = n+1;
  }
  sprprp_(sam.madof, minex, sam.mpmnpc, sam.mmnpc,
	  sam.mpmceq, sam.mmceq, sam.msc, NS, 6, sam.mpar,
	  mpar, sam.meqn, &iWork.front(), ierr);
  if (ierr < 0)
  {
    std::cerr <<"SAM::SPRPRM: Failure "<< ierr << std::endl;
    return;
  }
  else if (sam.neq != mpar[7])
  {
    std::cerr <<"SAM::SPRPRM: Internal error, NEQ = "<< sam.neq
	      <<" != "<< mpar[7] << std::endl;
    return;
  }
  else if (mpar[7] < 1)
    return; // No equations to solve

  // Perform symbolic assembly
  std::vector<int> itemp(mpar[34]);
  sprsas_(sam.mpar, sam.mpmnpc, sam.mmnpc, sam.madof, sam.msc, sam.mpmceq,
	  sam.mmceq, sam.meqn, mpar, &itemp.front(), &iWork.front(), 6, ierr);
  if (ierr < 0)
  {
    std::cerr <<"SAM::SPRSAS: Failure "<< ierr << std::endl;
    return;
  }

  // Reallocation of msica with the correct size
  msica = new int[mpar[1]];
  memcpy(msica,&itemp.front(),mpar[1]*sizeof(int));

  // Perform nodal reordering
  iWork.resize(mpar[36]);
  sprrnm_(mpar, msica, &iWork.front(), 6, ierr);
  if (ierr < 0)
  {
    std::cerr <<"SAM::SPRRNM: Failure "<< ierr << std::endl;
    return;
  }

  // Analyze the sparsity pattern
  mtrees = new int[mpar[35]];
  iWork.resize(mpar[37]);
  Real rinfo[4];
  sprtrs_(mpar, msica, mtrees, sam.meqn, sam.ndof,
	  &iWork.front(), rinfo, 6, ierr);
  if (ierr < 0)
  {
    std::cerr <<"SAM::SPRTRS: Failure "<< ierr << std::endl;
    return;
  }

  // Perform symbolic factorization
  itemp.resize(mpar[2]);
  iWork.resize(3*mpar[5] + 2*mpar[7] + 1);
  sprsmb_(mpar, msica, mtrees, &itemp.front(), sam.meqn,
	  &iWork.front(), 6, ierr);
  if (ierr < 0)
  {
    std::cerr <<"SAM::SPRSMB: Failure "<< ierr << std::endl;
    return;
  }

  // Reallocation of msifa with the correct size
  msifa = new int[mpar[2]];
  memcpy(msifa,&itemp.front(),mpar[2]*sizeof(int));

  // Finalize the SPR datastructure
  mvarnc = new int[2*mpar[7]];
  sprpmp_(msica, mtrees, sam.meqn, sam.mpar, mpar, mvarnc);

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


bool SPRMatrix::assemble (const Matrix& eM, const SAM& sam, int e)
{
#ifdef HAS_SPR
  iWork.resize(mpar[38]);
  int ierr;
  Real Bdummy;
  spradm_(eM.ptr(), sam.ttcc, sam.mpar, mpar,
	  sam.madof, sam.meqn, sam.mpmnpc, sam.mmnpc, sam.mpmceq, sam.mmceq,
	  msica, mtrees, msifa, mvarnc, values, &Bdummy, &iWork.front(),
	  e, eM.rows(), 6, 0, ierr);
  if (!ierr) return true;

  std::cerr <<"SAM::SPRADM: Failure "<< ierr << std::endl;
#endif
  return false;
}


bool SPRMatrix::assemble (const Matrix& eM, const SAM& sam,
			  SystemVector& B, int e)
{
  if (B.getType() != LinAlg::DENSE) return false;

#ifdef HAS_SPR
  iWork.resize(mpar[38]);
  int ierr;
  spradm_(eM.ptr(), sam.ttcc, sam.mpar, mpar,
	  sam.madof, sam.meqn, sam.mpmnpc, sam.mmnpc, sam.mpmceq, sam.mmceq,
	  msica, mtrees, msifa, mvarnc, values, B.getPtr(), &iWork.front(),
	  e, eM.rows(), 6, 1, ierr);
  if (!ierr) return true;

  std::cerr <<"SAM::SPRADM: Failure "<< ierr << std::endl;
#endif
  return false;
}


bool SPRMatrix::add (const SystemMatrix& B, Real alpha)
{
  const SPRMatrix* Bptr = dynamic_cast<const SPRMatrix*>(&B);
  if (!Bptr) return false;

  if (mpar[7] != Bptr->mpar[7] || mpar[15] != Bptr->mpar[15]) return false;

  for (int i = 0; i < mpar[7]+mpar[15]; i++)
    values[i] += alpha*Bptr->values[i];

  return true;
}


bool SPRMatrix::add (Real sigma)
{
#ifdef HAS_SPR
  int ierr;
  sprdad_(mpar, mtrees, msifa, values, sigma, 6, ierr);
  if (!ierr) return true;

  std::cerr <<"SAM::SPRDAD: Failure "<< ierr << std::endl;
#endif
  return false;
}


bool SPRMatrix::multiply (const SystemVector& B, SystemVector& C) const
{
  const int n = mpar[7];
  if (n < 1) return true;

  if (B.getType() != LinAlg::DENSE) return false;
  if (C.getType() != LinAlg::DENSE) return false;

  C.redim(n);
#ifdef HAS_SPR
  if (n > 0 && rWork.size() < (size_t)n) rWork.resize(n);
  int ierr;
  sprprm_(values+n, B.getRef(), C.getPtr(), &rWork.front(),
          mpar, mtrees, msifa, n, 1, 1, 1, 6, ierr);
  if (!ierr) return true;

  std::cerr <<"SPRMatrix::SPRPRM: Failure "<< ierr << std::endl;
#endif
  return false;
}


//! \brief Convenience macro.
#define MAX(a,b) (a) > (b) ? (a) : (b)

bool SPRMatrix::solve (SystemVector& B, bool, Real*)
{
  if (mpar[7] < 1) return true; // No equations to solve

  if (B.getType() != LinAlg::DENSE) return false;

#ifdef HAS_SPR
  Real tol[3] = { Real(1.0e-12), Real(0), Real(0) };
  iWork.resize(MAX(mpar[12],mpar[13]+1));
  rWork.resize(MAX(mpar[16],mpar[13]));
  int iop = mpar[0] < 5 ? 3 : 4;
  int ierr;
  sprsol_(iop, mpar, mtrees, msifa, values, B.getPtr(),
	  B.dim(), 1, tol, &iWork.front(), &rWork.front(), 6, ierr);
  if (!ierr) return true;

  std::cerr <<"SPRMatrix::SPRSOL: Failure "<< ierr << std::endl;
#endif
  return false;
}


bool SPRMatrix::solveEig (SPRMatrix& B, RealArray& val, Matrix& vec, int nv,
			  Real shift, int iop)
{
  const int n = mpar[7];
  if (n < 1 || nv == 0) return true; // No equations to solve
  if (n != B.mpar[7]) return false; // Incompatible matrices

#if HAS_SPR > 1
  std::cout <<"  Solving sparse eigenproblem using SAM::SPRLAX"<< std::endl;
  Real tol[3] = { Real(1.0e-8), Real(1.0e-12), Real(1.0e-8) };
  int maxlan = 3*nv+12;
  if (maxlan > n) maxlan = n;
  int ierr = 0;
  val.resize(n);
  vec.resize(n,maxlan);
  iWork.resize(MAX(2*maxlan+mpar[12],n));
  rWork.resize(MAX(maxlan*(maxlan+7),MAX(mpar[13],mpar[16])) + maxlan + 2*n);
  sprlax_(values, B.values, tol, mpar, mtrees, msifa, B.mpar, B.mtrees, B.msifa,
	  iop, &val.front(), vec.ptr(), &rWork.front(), &iWork.front(),
  	  shift, n, nv, maxlan, 6, 0, ierr);
  val.resize(nv);
  vec.resize(n,nv);
  if (!ierr) return true;

  std::cerr <<"SAM::SPRLAX: Failure "<< ierr << std::endl;
#else
  std::cerr <<"SPRMatrix: SPRLAX eigensolver not available"<< std::endl;
#endif
  return false;
}


Real SPRMatrix::Linfnorm () const
{
  std::cerr <<"SPRMatrix::Linfnorm not yet implemented"<< std::endl;
  return Real(0);
}
