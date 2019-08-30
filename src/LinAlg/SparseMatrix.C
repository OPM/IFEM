// $Id$
//==============================================================================
//!
//! \file SparseMatrix.C
//!
//! \date Jan 8 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of the system matrix on an unstructured sparse format.
//!
//==============================================================================

#include "SparseMatrix.h"
#include "IFEM.h"
#include "SAM.h"
#if defined(HAS_SUPERLU_MT)
#include "slu_mt_ddefs.h"
#elif defined(HAS_SUPERLU)
#include "slu_ddefs.h"
#endif
#ifdef HAS_SAMG
#include "samg.h"
#endif
#ifdef HAS_UMFPACK
#include <umfpack.h>
#endif
#ifdef USE_OPENMP
#include <omp.h>
#endif
#include <algorithm>

#if defined(HAS_SUPERLU_MT)
#define sluop_t superlumt_options_t
#elif defined(HAS_SUPERLU)
#define sluop_t superlu_options_t
#endif


/*!
  \brief Data structures for the SuperLU equation solver.
*/

struct SuperLUdata
{
#if defined(HAS_SUPERLU) || defined(HAS_SUPERLU_MT)
  SuperMatrix A; //!< The unfactored coefficient matrix
  SuperMatrix L; //!< The lower triangle factor
  SuperMatrix U; //!< The upper triangle factor
  Real*       R; //!< The row scale factors for \a A
  Real*       C; //!< The column scale factors for \a A
  int*   perm_r; //!< Row permutation vector
  int*   perm_c; //!< Column permutation vector
  int*    etree; //!< The elimination tree
  sluop_t* opts; //!< Input options for the SuperLU driver routine
#ifdef HAS_SUPERLU_MT
  equed_t equed; //!< Specifies the form of equilibration that was done
#else
  char equed[1]; //!< Specifies the form of equilibration that was done
#endif
  Real    rcond; //!< Reciprocal condition number
  Real      rpg; //!< Reciprocal pivot growth

  //! \brief The constructor initializes the default input options.
  explicit SuperLUdata(int numThreads = 0) :
    A{}, L{}, U{}
  {
    equed[0] = 0;
    R = C = 0;
    perm_r = perm_c = etree = 0;
    rcond = rpg = 0.0;
    if (numThreads > 0)
    {
      opts = new sluop_t;
#ifdef HAS_SUPERLU_MT
      opts->nprocs = numThreads;
      opts->fact = DOFACT;
      opts->trans = NOTRANS;
      opts->refact = NO;
      opts->panel_size = sp_ienv(1);
      opts->relax = sp_ienv(2);
      opts->diag_pivot_thresh = 1.0;
      opts->drop_tol = 0.0;
      opts->ColPerm = MMD_ATA;
      opts->usepr = NO;
      opts->SymmetricMode = NO;
      opts->PrintStat = NO;
      opts->perm_c = 0;
      opts->perm_r = 0;
      opts->work = 0;
      opts->lwork = 0;
      opts->etree = 0;
      opts->colcnt_h = 0;
      opts->part_super_h = 0;
#else
      set_default_options(opts);
      opts->SymmetricMode = YES;
      opts->ColPerm = MMD_AT_PLUS_A;
      opts->DiagPivotThresh = 0.001;
#endif
    }
    else
      opts = 0;
  }

  //! \brief No copying of this class.
  SuperLUdata(const SuperLUdata&) = delete;

  //! \brief The destructor frees the dynamically allocated data members.
  ~SuperLUdata()
  {
    Destroy_SuperMatrix_Store(&A);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    if (R)      delete[] R;
    if (C)      delete[] C;
    if (perm_r) delete[] perm_r;
    if (perm_c) delete[] perm_c;
    if (etree)  delete[] etree;
#ifdef HAS_SUPERLU_MT
    delete[] opts->etree;
    delete[] opts->colcnt_h;
    delete[] opts->part_super_h;
#endif
    if (opts)   delete   opts;
  }
#endif
};


bool SparseMatrix::printSLUstat = false;


SparseMatrix::SparseMatrix (SparseSolver eqSolver, int nt)
{
  editable = 'P';
  factored = false;
  nrow = ncol = 0;
  solver = eqSolver;
  numThreads = nt;
#ifdef HAS_UMFPACK
  umfSymbolic = nullptr;
#endif
  slu = 0;
}


SparseMatrix::SparseMatrix (size_t m, size_t n)
{
  editable = 'P';
  factored = false;
  nrow = m;
  ncol = n > 0 ? n : m;
  solver = NONE;
  numThreads = 0;
  slu = 0;
#ifdef HAS_UMFPACK
  umfSymbolic = nullptr;
#endif
}


SparseMatrix::SparseMatrix (const SparseMatrix& B) :
  elem(B.elem), IA(B.IA), JA(B.JA), A(B.A)
{
  editable = B.editable;
  factored = false;
  nrow = B.nrow;
  ncol = B.ncol;
  solver = B.solver;
  numThreads = B.numThreads;
  slu = 0; // The SuperLU data (if any) is not copied
#ifdef HAS_UMFPACK
  umfSymbolic = nullptr;
#endif
}


SparseMatrix::~SparseMatrix ()
{
  if (slu) delete slu;
#ifdef HAS_UMFPACK
  if (umfSymbolic)
    umfpack_di_free_symbolic(&umfSymbolic);
#endif
}


LinAlg::MatrixType SparseMatrix::getType () const
{
  return solver == S_A_M_G ? LinAlg::SAMG : LinAlg::SPARSE;
}


bool SparseMatrix::lockPattern (bool doLock)
{
  bool wasLocked = editable != 'P';
  if (editable) editable = doLock ? 'V' : 'P';
  return wasLocked;
}


void SparseMatrix::resize (size_t r, size_t c, bool forceEditable)
{
  factored = false;
  if (r == nrow && c == ncol && !forceEditable)
  {
    // Clear the matrix content but retain its sparsity pattern
    for (ValueMap::iterator it = elem.begin(); it != elem.end(); ++it)
      it->second = Real(0);
    std::fill(A.begin(),A.end(),Real(0));
    return;
  }

  // Clear the matrix completely, including its sparsity pattern
  editable = 'P';
  elem.clear();
  IA.clear();
  JA.clear();
  A.clear();

  nrow = r;
  ncol = c > 0 ? c : r;

  if (slu) delete slu;
  slu = 0;
#ifdef HAS_UMFPACK
  if (umfSymbolic) {
    umfpack_di_free_symbolic(&umfSymbolic);
    umfSymbolic = nullptr;
  }
#endif
}


bool SparseMatrix::redim (size_t r, size_t c)
{
  if (editable != 'P') return false;

  nrow = r;
  ncol = c > 0 ? c : r;
  for (ValueIter it = elem.begin(); it != elem.end();)
    if (it->first.first > nrow || it->first.second > ncol)
    {
      ValueIter jt = it++;
      elem.erase(jt->first);
    }
    else
      ++it;

  return true;
}


size_t SparseMatrix::dim (int idim) const
{
  if (idim == 1)
    return nrow;
  else if (idim == 2)
    return ncol;
  else if (idim == 3)
    return nrow*ncol;
  else
    return this->size();
}


Real& SparseMatrix::operator () (size_t r, size_t c)
{
  if (r < 1 || r > nrow || c < 1 || c > ncol)
    std::cerr <<"SparseMatrix::operator(): Indices ("
              << r <<","<< c <<") out of range "
              << nrow <<"x"<< ncol << std::endl;
  else if (editable) {
    IJPair index(r,c);
    ValueMap::iterator vit = elem.find(index);
    if (vit != elem.end())
      return vit->second; // This non-zero term already exists
    else if (editable == 'P') {
      // Editable pattern, insert a new non-zero entry
      Real& value = elem[index] = Real(0);
      return value;
    }
  }
  else if (solver == SUPERLU || solver == UMFPACK) {
    // Column-oriented format with 0-based indices
    IntVec::const_iterator begin = JA.begin() + IA[c-1];
    IntVec::const_iterator end = JA.begin() + IA[c];
    IntVec::const_iterator it = std::find(begin, end, r-1);
    if (it != end) return A[it - JA.begin()];
  }
  else {
    // Row-oriented format with 1-based indices
    IntVec::const_iterator begin = JA.begin() + (IA[r-1]-1);
    IntVec::const_iterator end = JA.begin() + (IA[r]-1);
    IntVec::const_iterator it = std::find(begin, end, c);
    if (it != end) return A[it - JA.begin()];
  }

  // If we arrive here, we have tried to update the sparsity pattern when it is
  // locked. The behavior then is unpredictable, especially when multithreading.
  std::cerr <<" *** Non-existing SparseMatrix entry (r,c)="<< r <<","<< c
            << std::endl;
  static Real anyValue = Real(0);
#if INDEX_CHECK > 1
  abort();
#endif
  return anyValue;
}


const Real& SparseMatrix::operator () (size_t r, size_t c) const
{
  if (r < 1 || r > nrow || c < 1 || c > ncol)
    std::cerr <<"SparseMatrix::operator(): Indices ("
              << r <<","<< c <<") out of range "
              << nrow <<"x"<< ncol << std::endl;
  else if (editable) {
    ValueIter vit = elem.find(IJPair(r,c));
    if (vit != elem.end()) return vit->second;
  }
  else if (solver == SUPERLU || solver == UMFPACK) {
    // Column-oriented format with 0-based indices
    IntVec::const_iterator begin = JA.begin() + IA[c-1];
    IntVec::const_iterator end = JA.begin() + IA[c];
    IntVec::const_iterator it = std::find(begin, end, r-1);
    if (it != end) return A[it - JA.begin()];
  }
  else {
    // Row-oriented format with 1-based indices
    IntVec::const_iterator begin = JA.begin() + (IA[r-1]-1);
    IntVec::const_iterator end = JA.begin() + (IA[r]-1);
    IntVec::const_iterator it = std::find(begin, end, c);
    if (it != end) return A[it - JA.begin()];
  }

  // Return zero for any non-existing non-zero term
  static const Real zero = Real(0);
  return zero;
}


void SparseMatrix::dump (std::ostream& os, char format, const char* label)
{
  if (label) os << label <<" = [\n";
  switch (format)
    {
    case 'M':
    case 'm':
      if (editable)
        for (const auto& it : elem)
          os << it.first.first <<' '<< it.first.second <<" "<< it.second
             <<";\n";
      else if (solver == SUPERLU || solver == UMFPACK) {
        // Column-oriented format with 0-based indices
        os << JA[0]+1 <<" 1 "<< A[0];
        for (size_t j = 1; j <= ncol; j++)
          for (int i = IA[j-1]; i < IA[j]; i++) {
            if(j==1 && i==IA[0]) continue;
            os << ";\n" << JA[i]+1 <<' '<< j <<' '<< A[i] ;
        }
      } else {
        // Row-oriented format with 1-based indices
        os << "1 " << JA[0] <<' '<< A[0] ;
        for (size_t i = 1; i <= nrow; i++)
          for (int j = IA[i-1]; j < IA[i]; j++) {
            if(i==1 && j==IA[0]) continue;
            os << ";\n" << i <<' '<< JA[j-1] <<' '<< A[i-1] ;
        }
      }
      os << "];\n";
      break;

    default:
      this->write(os);
    }
}


std::ostream& SparseMatrix::write (std::ostream& os) const
{
  os << nrow <<' '<< ncol <<' '<< this->size();
  if (editable)
    for (const auto& it : elem)
      os <<'\n'<< it.first.first <<' '<< it.first.second <<" : "<< it.second;
  else {
    size_t i;
    os <<'\n';
    for (i = 0; i < A.size(); i++) os << A[i] <<' ';

    os <<'\n'<< IA.size() <<'\n';
    for (i = 0; i < IA.size(); i++) os << IA[i] <<' ';

    os <<'\n'<< JA.size() <<'\n';
    for (i = 0; i < JA.size(); i++) os << JA[i] <<' ';
  }
  return os << std::endl;
}


void SparseMatrix::printSparsity (std::ostream& os) const
{
  if (nrow < 1 || ncol < 1) return;

  size_t r, c;
  os <<'\t';
  for (c = 1; c <= ncol; c++)
    (c%10) ? os << c%10 : os << ' ';
  os <<'\n';

  for (r = 1; r <= nrow; r++) {
    os << r <<'\t';
    for (c = 1; c <= ncol; c++)
      if (editable)
        os << (elem.find(IJPair(r,c)) == elem.end() ? '.' : 'X');
      else if (solver == SUPERLU || solver == UMFPACK) {
        // Column-oriented format with 0-based indices
        IntVec::const_iterator begin = JA.begin() + IA[c-1];
        IntVec::const_iterator end = JA.begin() + IA[c];
        os << (std::find(begin,end,r-1) == end ? '.' : 'X');
      }
      else {
        // Row-oriented format with 1-based indices
        IntVec::const_iterator begin = JA.begin() + (IA[r-1]-1);
        IntVec::const_iterator end = JA.begin() + (IA[r]-1);
        os << (std::find(begin,end,c) == end ? '.' : 'X');
      }
    os <<'\n';
  }
  os << std::endl;
}


void SparseMatrix::printFull (std::ostream& os) const
{
  for (size_t r = 1; r <= nrow; r++)
    for (size_t c = 1; c <= ncol; c++)
      os << this->operator()(r,c) << (c < ncol ? '\t' : '\n');

  os << std::endl;
}


bool SparseMatrix::augment (const SystemMatrix& B, size_t r0, size_t c0)
{
  if (editable != 'P') return false;

  const SparseMatrix* Bptr = dynamic_cast<const SparseMatrix*>(&B);
  if (!Bptr) return false;
  if (!Bptr->editable) return false;

  if (r0+Bptr->nrow > nrow) nrow = r0 + Bptr->nrow;
  if (r0+Bptr->nrow > ncol) ncol = r0 + Bptr->nrow;
  if (c0+Bptr->ncol > ncol) ncol = c0 + Bptr->ncol;
  if (c0+Bptr->ncol > nrow) nrow = c0 + Bptr->ncol;

  for (const auto& it : Bptr->elem)
  {
    elem[std::make_pair(r0+it.first.first,c0+it.first.second)] += it.second;
    elem[std::make_pair(c0+it.first.second,r0+it.first.first)] += it.second;
  }

  return true;
}


bool SparseMatrix::truncate (Real threshold)
{
  if (!editable) return false; // Not implemented for non-editable matrices yet

  Real tol = Real(0);
  for (const auto& it : elem)
    if (it.first.first == it.first.second)
      if (it.second > tol)
        tol = it.second;
      else if (it.second < -tol)
        tol = -it.second;

  tol *= threshold;
  size_t nnz = elem.size();
  for (ValueMap::iterator it = elem.begin(); it != elem.end();)
    if (it->second <= -tol || it->second >= tol)
      ++it;
    else if (editable == 'P')
    {
      ValueIter jt = it++;
      elem.erase(jt->first);
    }
    else
      (it++)->second = Real(0);

  if (nnz > elem.size())
    IFEM::cout <<"SparseMatrix: Truncated "<< nnz-elem.size()
               <<" elements smaller than "<< tol <<" to zero"<< std::endl;
  return true;
}


bool SparseMatrix::add (const SystemMatrix& B, Real alpha)
{
  const SparseMatrix* Bptr = dynamic_cast<const SparseMatrix*>(&B);
  if (!Bptr) return false;

  if (Bptr->nrow > nrow || Bptr->ncol > ncol) return false;

  if (editable == 'P' && Bptr->editable)
    for (const auto& it : Bptr->elem)
      elem[it.first] += alpha*it.second;

  else if (!editable && !Bptr->editable)
  {
    // For non-editable matrices the sparsity patterns must match
    if (A.size() == Bptr->A.size() && IA == Bptr->IA && JA == Bptr->JA)
      A.add(Bptr->A,alpha);
    else
      return false;
  }
  else if (editable == 'P')
  {
    if (solver == SUPERLU || solver == UMFPACK)
      // Column-oriented format with 0-based indices
      for (size_t j = 1; j <= Bptr->ncol; j++)
        for (int i = Bptr->IA[j-1]; i < Bptr->IA[j]; i++)
          elem[IJPair(Bptr->JA[i]+1,j)] += alpha*Bptr->A[i];
    else
      // Row-oriented format with 1-based indices
      for (size_t i = 1; i <= Bptr->nrow; i++)
        for (int j = Bptr->IA[i-1]; j < Bptr->IA[i]; j++)
          elem[IJPair(i,Bptr->JA[j-1])] += alpha*Bptr->A[j-1];
  }
  else
    return false;

  return true;
}


bool SparseMatrix::add (Real sigma)
{
  for (size_t i = 1; i <= nrow && i <= ncol; i++)
    this->operator()(i,i) += sigma;

  return true;
}


bool SparseMatrix::multiply (const SystemVector& B, SystemVector& C) const
{
  C.resize(nrow,true);
  if (B.dim() < ncol) return false;

  const StdVector* Bptr = dynamic_cast<const StdVector*>(&B);
  if (!Bptr) return false;
  StdVector*       Cptr = dynamic_cast<StdVector*>(&C);
  if (!Cptr) return false;

  if (editable)
    for (const auto& it : elem)
      (*Cptr)(it.first.first) += it.second*(*Bptr)(it.first.second);
  else if (solver == SUPERLU) {
#ifdef notyet_USE_OPENMP // TODO: akva needs to fix this, gives wrong result!
    if (omp_get_max_threads() > 1) {
      std::vector<Vector> V(omp_get_max_threads());
#pragma omp parallel for schedule(static)
      for (size_t j = 1; j <= ncol; j++) {
        Vector& myV = V[omp_get_thread_num()];
        if (myV.empty())
          myV.resize(nrow);
        for (int i = IA[j-1]; i < IA[j]; i++)
          myV(JA[i]+1) += A[i]*(*Bptr)(j);
      }
      for (int j = 0; j < omp_get_max_threads(); j++)
        for (size_t i = 1; i <= V[j].size(); i++)
          (*Cptr)(i) += V[j](i);
    } else
#endif
    // Column-oriented format with 0-based indices
    for (size_t j = 1; j <= ncol; j++)
      for (int i = IA[j-1]; i < IA[j]; i++)
        (*Cptr)(JA[i]+1) += A[i]*(*Bptr)(j);
  }
  else // Row-oriented format with 1-based indices
    for (size_t i = 1; i <= nrow; i++)
      for (int j = IA[i-1]; j < IA[i]; j++)
        (*Cptr)(i) += A[j-1]*(*Bptr)(JA[j-1]);

  return true;
}


/*!
  \brief This is a C++ version of the F77 subroutine ADDEM2 (SAM library).
  \details It performs exactly the same tasks, except that \a NRHS always is 1,
  and that the system matrix \a SM here is an object of the SparseMatrix class.
*/

static void assemSparse (const Matrix& eM, SparseMatrix& SM, Vector& SV,
                         const IntVec& meen, const int* meqn,
                         const int* mpmceq, const int* mmceq, const Real* ttcc)
{
  // Add elements corresponding to free dofs in eM into SM
  int i, j, ip, nedof = meen.size();
  for (j = 1; j <= nedof; j++)
  {
    int jeq = meen[j-1];
    if (jeq < 1) continue;

    SM(jeq,jeq) += eM(j,j);

    for (i = 1; i < j; i++)
    {
      int ieq = meen[i-1];
      if (ieq < 1) continue;

      SM(ieq,jeq) += eM(i,j);
      SM(jeq,ieq) += eM(j,i);
    }
  }

  // Add (appropriately weighted) elements corresponding to constrained
  // (dependent and prescribed) dofs in eM into SM and/or SV
  for (j = 1; j <= nedof; j++)
  {
    int jceq = -meen[j-1];
    if (jceq < 1) continue;

    int jp = mpmceq[jceq-1];
    Real c0 = ttcc[jp-1];

    // Add contributions to SV (right-hand-side)
    if (!SV.empty())
      for (i = 1; i <= nedof; i++)
      {
        int ieq = meen[i-1];
        int iceq = -ieq;
        if (ieq > 0)
          SV(ieq) -= c0*eM(i,j);
        else if (iceq > 0)
          for (ip = mpmceq[iceq-1]; ip < mpmceq[iceq]-1; ip++)
            if (mmceq[ip] > 0)
            {
              ieq = meqn[mmceq[ip]-1];
              SV(ieq) -= c0*ttcc[ip]*eM(i,j);
            }
      }

    // Add contributions to SM
    for (jp = mpmceq[jceq-1]; jp < mpmceq[jceq]-1; jp++)
      if (mmceq[jp] > 0)
      {
        int jeq = meqn[mmceq[jp]-1];
        for (i = 1; i <= nedof; i++)
        {
          int ieq = meen[i-1];
          int iceq = -ieq;
          if (ieq > 0)
          {
            SM(ieq,jeq) += ttcc[jp]*eM(i,j);
            SM(jeq,ieq) += ttcc[jp]*eM(j,i);
          }
          else if (iceq > 0)
            for (ip = mpmceq[iceq-1]; ip < mpmceq[iceq]-1; ip++)
              if (mmceq[ip] > 0)
              {
                ieq = meqn[mmceq[ip]-1];
                SM(ieq,jeq) += ttcc[ip]*ttcc[jp]*eM(i,j);
              }
        }
      }
  }
}


/*!
  \brief Adds a nodal vector into a non-symmetric rectangular sparse matrix.
  \details The nodal values are added into the columns \a col to \a col+2.
*/

static void assemSparse (const RealArray& V, SparseMatrix& SM, size_t col,
                         const IntVec& mnen, const int* meqn,
                         const int* mpmceq, const int* mmceq, const Real* ttcc)
{
  for (size_t d = 0; d < mnen.size(); d++, col++)
  {
    Real vd = d < V.size() ? V[d] : V.back();
    int ieq = mnen[d];
    int ceq = -ieq;
    if (ieq > 0)
      SM(ieq,col) += vd;
    else if (ceq > 0)
      for (int ip = mpmceq[ceq-1]; ip < mpmceq[ceq]-1; ip++)
      {
        ieq = meqn[mmceq[ip]-1];
        SM(ieq,col) += vd;
      }
  }
}


void SparseMatrix::initAssembly (const SAM& sam, bool delayLocking)
{
  this->resize(sam.neq,sam.neq);
#ifdef USE_OPENMP
  if (omp_get_max_threads() > 1)
    this->preAssemble(sam,delayLocking);
#endif
}


void SparseMatrix::preAssemble (const SAM& sam, bool delayLocking)
{
  if (editable != 'P')
    return;

  // Compute the sparsity pattern
  std::vector<IntSet> dofc;
  if (!sam.getDofCouplings(dofc))
    return;

  // If we are not locking the sparsity pattern yet, the index pair map over
  // the non-zero matrix elements needs to be initialized before the assembly.
  // This is used when SAM::getDofCouplings does not return all connectivities
  if (delayLocking) // that will exist in the final matrix.
    for (size_t i = 0; i < dofc.size(); i++)
      for (const int& it : dofc[i])
        (*this)(i+1,it) = 0.0;

  editable = 'V'; // Temporarily lock the sparsity pattern
  if (delayLocking)
    return; // The final sparsity pattern is not fixed yet

  IFEM::cout <<"\nPre-computing sparsity pattern for system matrix ("
             << nrow <<"x"<< ncol <<"): "<< std::flush;

  switch (solver) {
  case UMFPACK:
  case SUPERLU: this->optimiseSLU(dofc); break;
  case S_A_M_G: this->optimiseSAMG(); break;
  default: break;
  }

  // The sparsity pattern is now permanently locked (until resize is invoked)
  IFEM::cout <<"nNZ = "<< this->size() << std::endl;
}


void SparseMatrix::preAssemble (const std::vector<IntVec>& MMNPC, size_t nel)
{
#ifdef USE_OPENMP
  if (omp_get_max_threads() < 2)
    return;

  // Compute the nodal sparsity pattern
  int inod, jnod;
  for (size_t iel = 0; iel < nel; iel++)
    for (size_t j = 0; j < MMNPC[iel].size(); j++)
      if ((jnod = MMNPC[iel][j]+1) > 0)
      {
        (*this)(jnod,jnod) = 0.0;
        for (size_t i = 0; i < j; i++)
          if ((inod = MMNPC[iel][i]+1) > 0)
            (*this)(inod,jnod) = (*this)(jnod,inod) = 0.0;
      }

  switch (solver) {
  case UMFPACK:
  case SUPERLU: this->optimiseSLU(); break;
  case S_A_M_G: this->optimiseSAMG(); break;
  default: break;
  }
#endif
}


void SparseMatrix::init ()
{
  this->resize(nrow,ncol);
}


bool SparseMatrix::assemble (const Matrix& eM, const SAM& sam, int e)
{
  IntVec meen;
  if (!sam.getElmEqns(meen,e,eM.rows()))
    return false;

  Vector dummyB;
  assemSparse(eM,*this,dummyB,meen,sam.meqn,sam.mpmceq,sam.mmceq,sam.ttcc);
  return true;
}


bool SparseMatrix::assemble (const Matrix& eM, const SAM& sam,
                             SystemVector& B, int e)
{
  StdVector* Bptr = dynamic_cast<StdVector*>(&B);
  if (!Bptr) return false;

  IntVec meen;
  if (!sam.getElmEqns(meen,e,eM.rows()))
    return false;

  assemSparse(eM,*this,*Bptr,meen,sam.meqn,sam.mpmceq,sam.mmceq,sam.ttcc);
  return true;
}


bool SparseMatrix::assemble (const Matrix& eM, const SAM& sam,
                             SystemVector& B, const IntVec& meen)
{
  StdVector* Bptr = dynamic_cast<StdVector*>(&B);
  if (!Bptr) return false;

  if (eM.rows() < meen.size() || eM.cols() < meen.size())
    return false;

  assemSparse(eM,*this,*Bptr,meen,sam.meqn,sam.mpmceq,sam.mmceq,sam.ttcc);
  return true;
}


bool SparseMatrix::assembleCol (const RealArray& V, const SAM& sam,
                                int n, size_t col)
{
  if (V.empty() || col > ncol) return false;

  IntVec mnen;
  if (!sam.getNodeEqns(mnen,n)) return false;

  assemSparse(V,*this,col,mnen,sam.meqn,sam.mpmceq,sam.mmceq,sam.ttcc);
  return true;
}


bool SparseMatrix::optimiseSAMG (bool transposed)
{
  if (!editable) return false;

  ValueMap trans; // only computed if needed
  ValueIter begin, end;

  if (transposed) {
    for (const auto& it : elem)
      trans[IJPair(it.first.second,it.first.first)] = it.second;
    begin = trans.begin();
    end = trans.end();
    std::swap(nrow,ncol);
  }
  else {
    begin = elem.begin();
    end = elem.end();
  }

  size_t nnz = this->size();
  A.resize(nnz);
  JA.resize(nnz);
  IA.resize(nrow+1,nnz+1);

  IA[0] = 1; // first row start at index 1
  size_t cur_row = 1, ix = 0;
  for (ValueIter it = begin; it != end; ++it, ix++) {
    A[ix] = it->second; // storing element value
    JA[ix] = it->first.second;
    while (it->first.first > cur_row)
      IA[cur_row++] = ix+1;
  }

  editable = false;
  elem.clear(); // Erase the editable matrix elements

  // convert to row storage format required by SAMG (diagonal term first)
  for (size_t r = 0; r < nrow; r++) {
    int rstart = IA[r]-1;
    int rend = IA[r+1]-1;
    // looking for diagonal element
    for (int diag_ix = rstart; diag_ix < rend; diag_ix++)
      if (JA[diag_ix] == (int)(1+r)) {
        // swapping (if necessary) with first element on this row
        if (diag_ix > rstart) {
          std::swap(A[rstart],A[diag_ix]);
          std::swap(JA[rstart],JA[diag_ix]);
        }
        break;
      }
  }

  return true;
}


/*!
  This method is based on the function dreadtriple() from the SuperLU package.
*/

bool SparseMatrix::optimiseSLU ()
{
  if (!editable) return false;

  size_t nnz = this->size();
  A.resize(nnz);
  JA.resize(nnz);
  IA.resize(ncol+1,0);

  // Initialize the array of column pointers
  ValueIter it;
  for (const auto& it : elem)
    if (it.first.first <= nrow && it.first.second <= ncol)
      IA[it.first.second-1]++;
    else
      return false;

  size_t j, nz;
  int k, jsize = IA[0];
  for (j = 1, k = IA[0] = 0; j < ncol; j++) {
    k += jsize;
    jsize = IA[j];
    IA[j] = k;
  }

  // Copy the triplets into the column-oriented storage
  for (nz = 0, it = elem.begin(); nz < nnz; nz++, ++it) {
    k = IA[it->first.second-1]++;
    JA[k] = it->first.first-1;
    A[k]  = it->second;
  }

  // Reset the column pointers to the beginning of each column
  for (j = ncol; j > 0; j--)
    IA[j] = IA[j-1];
  IA[0] = 0;

  editable = false;
  elem.clear(); // Erase the editable matrix elements

  return true;
}


/*!
  This method does not use the internal index-pair to value map \a elem.
*/

bool SparseMatrix::optimiseSLU (const std::vector<IntSet>& dofc)
{
  if (!editable) return false;

  // Initialize the array of column pointers
  size_t i, j, nnz = 0;
  IA.resize(ncol+1,0);
  for (i = 0; i < dofc.size(); i++)
  {
    nnz += dofc[i].size();
    for (const int& it : dofc[i])
      if (i < nrow && it <= (int)ncol)
        IA[it-1]++;
      else
        return false;
  }

  int k, jsize = IA[0];
  for (j = 1, k = IA[0] = 0; j < ncol; j++) {
    k += jsize;
    jsize = IA[j];
    IA[j] = k;
  }

  // Initialize the array of row indices
  JA.resize(nnz);
  for (i = 0; i < dofc.size(); i++)
    for (const int& it : dofc[i])
      JA[IA[it-1]++] = i;

  // Reset the column pointers to the beginning of each column
  for (j = ncol; j > 0; j--)
    IA[j] = IA[j-1];
  IA[0] = 0;

  editable = false;
  A.resize(nnz); // Allocate the non-zero matrix element storage

  return true;
}


bool SparseMatrix::solve (SystemVector& B, bool, Real* rc)
{
  if (this->size() < 1) return true; // No equations to solve

  StdVector* Bptr = dynamic_cast<StdVector*>(&B);
  if (!Bptr) return false;

  switch (solver)
    {
    case SUPERLU: return this->solveSLUx(*Bptr,rc);
    case S_A_M_G: return this->solveSAMG(*Bptr);
    case UMFPACK: return this->solveUMF(*Bptr,rc);
    default: std::cerr <<"SparseMatrix::solve: No equation solver"<< std::endl;
    }

  return false;
}


bool SparseMatrix::solveSLU (Vector& B)
{
  int ierr = ncol+1;
  if (!factored) this->optimiseSLU();

#ifdef HAS_SUPERLU_MT
  if (!slu) {
    // Create a new SuperLU matrix
    slu = new SuperLUdata;
    slu->perm_c = new int[ncol];
    slu->perm_r = new int[nrow];
    dCreate_CompCol_Matrix(&slu->A, nrow, ncol, this->size(),
                           &A.front(), &JA.front(), &IA.front(),
                           SLU_NC, SLU_D, SLU_GE);
  }
  else {
    Destroy_SuperMatrix_Store(&slu->A);
    Destroy_SuperNode_Matrix(&slu->L);
    Destroy_CompCol_Matrix(&slu->U);
    dCreate_CompCol_Matrix(&slu->A, nrow, ncol, this->size(),
                           &A.front(), &JA.front(), &IA.front(),
                           SLU_NC, SLU_D, SLU_GE);
  }

  // Get column permutation vector perm_c[], according to permc_spec:
  //   permc_spec = 0: natural ordering
  //   permc_spec = 1: minimum degree ordering on structure of A'*A
  //   permc_spec = 2: minimum degree ordering on structure of A'+A
  //   permc_spec = 3: approximate minimum degree for unsymmetric matrices
  int permc_spec = 1;
  get_perm_c(permc_spec, &slu->A, slu->perm_c);

  // Create right-hand-side/solution vector(s)
  size_t nrhs = B.size() / nrow;
  SuperMatrix Bmat;
  dCreate_Dense_Matrix(&Bmat, nrow, nrhs, B.ptr(), nrow,
                       SLU_DN, SLU_D, SLU_GE);

  // Invoke the simple driver
  pdgssv(numThreads, &slu->A, slu->perm_c, slu->perm_r,
         &slu->L, &slu->U, &Bmat, &ierr);

  if (ierr > 0)
    std::cerr <<"SuperLU_MT Failure "<< ierr << std::endl;

  Destroy_SuperMatrix_Store(&Bmat);

#elif defined(HAS_SUPERLU)
  if (!slu) {
    // Create a new SuperLU matrix
    slu = new SuperLUdata(1);
    slu->perm_c = new int[ncol];
    slu->perm_r = new int[nrow];
    dCreate_CompCol_Matrix(&slu->A, nrow, ncol, this->size(),
                           &A.front(), &JA.front(), &IA.front(),
                           SLU_NC, SLU_D, SLU_GE);
  }
  else if (factored)
    slu->opts->Fact = FACTORED; // Re-use previous factorization
  else {
    Destroy_SuperMatrix_Store(&slu->A);
    Destroy_SuperNode_Matrix(&slu->L);
    Destroy_CompCol_Matrix(&slu->U);
    dCreate_CompCol_Matrix(&slu->A, nrow, ncol, this->size(),
                           &A.front(), &JA.front(), &IA.front(),
                           SLU_NC, SLU_D, SLU_GE);
  }

  // Create right-hand-side/solution vector(s)
  size_t nrhs = B.size() / nrow;
  SuperMatrix Bmat;
  dCreate_Dense_Matrix(&Bmat, nrow, nrhs, B.ptr(), nrow,
                       SLU_DN, SLU_D, SLU_GE);

  SuperLUStat_t stat;
  StatInit(&stat);

  // Invoke the simple driver
  dgssv(slu->opts, &slu->A, slu->perm_c, slu->perm_r,
        &slu->L, &slu->U, &Bmat, &stat, &ierr);

  if (ierr > 0)
    std::cerr <<"SuperLU Failure "<< ierr << std::endl;
  else
    factored = true;

  if (printSLUstat)
    StatPrint(&stat);
  StatFree(&stat);

  Destroy_SuperMatrix_Store(&Bmat);
#else
  std::cerr <<"SparseMatrix::solve: SuperLU solver not available"<< std::endl;
#endif
  return ierr == 0;
}


bool SparseMatrix::solveSLUx (Vector& B, Real* rcond)
{
  int ierr = ncol+1;
  if (!factored) this->optimiseSLU();

#ifdef HAS_SUPERLU_MT
  if (!slu) {
    // Create a new SuperLU matrix
    slu = new SuperLUdata(numThreads);
    slu->equed = NOEQUIL;
    slu->perm_c = new int[ncol];
    slu->perm_r = new int[nrow];
    slu->C = new Real[ncol];
    slu->R = new Real[nrow];
    slu->opts->etree = new int[ncol];
    slu->opts->colcnt_h = new int[ncol];
    slu->opts->part_super_h = new int[ncol];
    memset(slu->opts->colcnt_h, 0, ncol*sizeof(int));
    memset(slu->opts->part_super_h, 0, ncol*sizeof(int));
    memset(slu->opts->etree, 0, ncol*sizeof(int));
    dCreate_CompCol_Matrix(&slu->A, nrow, ncol, this->size(),
                           &A.front(), &JA.front(), &IA.front(),
                           SLU_NC, SLU_D, SLU_GE);

    // Get column permutation vector perm_c[], according to permc_spec:
    //   permc_spec = 0: natural ordering
    //   permc_spec = 1: minimum degree ordering on structure of A'*A
    //   permc_spec = 2: minimum degree ordering on structure of A'+A
    //   permc_spec = 3: approximate minimum degree for unsymmetric matrices
    int permc_spec = 1;
    get_perm_c(permc_spec, &slu->A, slu->perm_c);
  }
  else if (factored)
    slu->opts->fact = FACTORED; // Re-use previous factorization
  else
    slu->opts->refact = YES; // Re-use previous ordering

  // Create right-hand-side and solution vector(s)
  Vector      X(B.size());
  SuperMatrix Bmat, Xmat;
  const size_t nrhs = B.size() / nrow;
  dCreate_Dense_Matrix(&Bmat, nrow, nrhs, B.ptr(), nrow,
                       SLU_DN, SLU_D, SLU_GE);
  dCreate_Dense_Matrix(&Xmat, nrow, nrhs, X.ptr(), nrow,
                       SLU_DN, SLU_D, SLU_GE);

  Real ferr[nrhs], berr[nrhs];
  superlu_memusage_t mem_usage;

  // Invoke the expert driver
  pdgssvx(numThreads, slu->opts, &slu->A, slu->perm_c, slu->perm_r,
          &slu->equed, slu->R, slu->C, &slu->L, &slu->U, &Bmat, &Xmat,
          &slu->rpg, &slu->rcond, ferr, berr, &mem_usage, &ierr);

  B.swap(X);

  if (ierr > 0)
    std::cerr <<"SuperLU_MT Failure "<< ierr << std::endl;
  else if (!factored)
  {
    factored = true;
    if (rcond)
      *rcond = slu->rcond;
  }

  Destroy_SuperMatrix_Store(&Bmat);
  Destroy_SuperMatrix_Store(&Xmat);

#elif defined(HAS_SUPERLU)
  if (!slu) {
    // Create a new SuperLU matrix
    slu = new SuperLUdata(1);
    slu->perm_c = new int[ncol];
    slu->perm_r = new int[nrow];
    slu->etree = new int[ncol];
    slu->C = new Real[ncol];
    slu->R = new Real[nrow];
    dCreate_CompCol_Matrix(&slu->A, nrow, ncol, this->size(),
                           &A.front(), &JA.front(), &IA.front(),
                           SLU_NC, SLU_D, SLU_GE);
  }
  else if (factored)
    slu->opts->Fact = FACTORED; // Re-use previous factorization
  else {
    Destroy_SuperMatrix_Store(&slu->A);
    Destroy_SuperNode_Matrix(&slu->L);
    Destroy_CompCol_Matrix(&slu->U);
    dCreate_CompCol_Matrix(&slu->A, nrow, ncol, this->size(),
                           &A.front(), &JA.front(), &IA.front(),
                           SLU_NC, SLU_D, SLU_GE);
  }

  // Create right-hand-side vector and solution vector
  Vector      X(B.size());
  SuperMatrix Bmat, Xmat;
  const  size_t nrhs = B.size() / nrow;
  dCreate_Dense_Matrix(&Bmat, nrow, nrhs, B.ptr(), nrow,
                       SLU_DN, SLU_D, SLU_GE);
  dCreate_Dense_Matrix(&Xmat, nrow, nrhs, X.ptr(), nrow,
                       SLU_DN, SLU_D, SLU_GE);

  slu->opts->ConditionNumber = printSLUstat || rcond ? YES : NO;
  slu->opts->PivotGrowth = printSLUstat ? YES : NO;

  void* work = 0;
  int  lwork = 0;
  Real ferr[nrhs], berr[nrhs];
  mem_usage_t mem_usage;

  SuperLUStat_t stat;
  StatInit(&stat);

  // Invoke the expert driver
#if SUPERLU_VERSION == 5
  GlobalLU_t Glu;
  dgssvx(slu->opts, &slu->A, slu->perm_c, slu->perm_r, slu->etree, slu->equed,
         slu->R, slu->C, &slu->L, &slu->U, work, lwork, &Bmat, &Xmat,
         &slu->rpg, &slu->rcond, ferr, berr, &Glu, &mem_usage, &stat, &ierr);
#else
  dgssvx(slu->opts, &slu->A, slu->perm_c, slu->perm_r, slu->etree, slu->equed,
         slu->R, slu->C, &slu->L, &slu->U, work, lwork, &Bmat, &Xmat,
         &slu->rpg, &slu->rcond, ferr, berr, &mem_usage, &stat, &ierr);
#endif

  B.swap(X);

  if (ierr > 0)
    std::cerr <<"SuperLU Failure "<< ierr << std::endl;
  else if (!factored)
  {
    factored = true;
    if (rcond)
      *rcond = slu->rcond;
  }

  if (printSLUstat)
  {
    StatPrint(&stat);
    IFEM::cout <<"Reciprocal condition number = "<< slu->rcond
               <<"\nReciprocal pivot growth = "<< slu->rpg << std::endl;
  }
  StatFree(&stat);

  Destroy_SuperMatrix_Store(&Bmat);
  Destroy_SuperMatrix_Store(&Xmat);
#else
  std::cerr <<"SparseMatrix::solve: SuperLU solver not available"<< std::endl;
#endif
  return ierr == 0;
}


bool SparseMatrix::solveUMF (Vector& B, Real* rcond)
{
  if (!factored) this->optimiseSLU();

#ifdef HAS_UMFPACK
  double info[UMFPACK_INFO];
  Vector X(B.size());
  if (!umfSymbolic) {
    umfpack_di_symbolic(nrow, ncol, IA.data(), JA.data(),
                        A.data(), &umfSymbolic, nullptr, info);
    if (info[UMFPACK_STATUS] != UMFPACK_OK)
      return false;
  }

  void* numeric;
  umfpack_di_numeric(IA.data(), JA.data(), A.data(), umfSymbolic,
                     &numeric, nullptr, info);
  if (info[UMFPACK_STATUS] != UMFPACK_OK)
    return false;
  if (rcond)
    *rcond = info[UMFPACK_RCOND];
  umfpack_di_solve(UMFPACK_A,
                   IA.data(), JA.data(), A.data(),
                   &X[0], &B[0], numeric, nullptr, info);
  if (info[UMFPACK_STATUS] == UMFPACK_OK)
    B = X;
  umfpack_di_free_numeric(&numeric);
  return info[UMFPACK_STATUS] == UMFPACK_OK;
#else
  std::cerr <<"SparseMatrix::solve: UMFPACK solver not available"<< std::endl;
  return false;
#endif
}


bool SparseMatrix::solveSAMG (Vector& B)
{
  int ierr = 1;
  if (!factored) this->optimiseSAMG();

#ifdef HAS_SAMG
  // Setting up additional parameters
  int nnu = this->rows(); // number of solution components (variables)
  int nna = this->size(); // number of nonzero entries in system matrix
  int nsys = 1; // one unknown - scalar system
  int iu = 1; // dummy - since nsys = 1;
  int ndiu = 1; // dummy - since nsys = 1;
  int ip = 1; // dummy - since nsys = 1 and no point-based approach used
  int ndip = 1; // dummy - since nsys = 1 and no point-based approach used
  int matrix = 120; // symmetric, NO zero rowsum, NO modification of matrix
  int iscale = 1; // dummy - since nsys = 1;
  Real res_in, res_out; // output parameters
  int ncyc_done; // output parameters
  int nsolve = 2; // default solution strategy
                  // (this parameter can quickly become complicated to set)
  int ifirst = 1; // no initial solution guess
  Real eps = -1.0e-16; // stopping criterion
  int ncyc = 110200; // define the cycling and accelleration strategy
  int iswitch = !factored ? 4140 : 1140; // No memory release in first run,
  // assuming identical coefficent matrix in all subsequent runs
  Real a_cmplx = 0.0; // l
  Real g_cmplx = 0.0; //  l  specifies preallocation of memory, etc.
  Real p_cmplx = 0.0; //  /
  Real w_avrge = 0.0; // /
  Real chktol = 0.0; // standard checking for logical correctness
                     // @@ set this to negative value for production run
  int iout = -2; // minimal output on results and timings
  int idump = -2; // printout of coarsening history

  int mode_mess = -2;
  SAMG_SET_MODE_MESS(&mode_mess);

  Vector X(B.size());

  SAMG(&nnu, &nna, &nsys,
       &IA.front(), &JA.front(), &A.front(), B.ptr(), X.ptr(),
       &iu, &ndiu, &ip, &ndip, &matrix, &iscale,
       &res_in, &res_out, &ncyc_done, &ierr,
       &nsolve, &ifirst, &eps, &ncyc, &iswitch,
       &a_cmplx, &g_cmplx, &p_cmplx, &w_avrge,
       &chktol, &idump, &iout);

  B.swap(X);

  if (ierr > 0)
    std::cerr <<"SAMG Failure "<< ierr << std::endl;
  else
  {
    factored = true;
    if (ierr < 0)
      std::cerr <<"SAMG warning "<< -ierr << std::endl;
  }
#else
  std::cerr <<"SparseMatrix::solve: SAMG solver not available"<< std::endl;
#endif
  return ierr <= 0;
}


Real SparseMatrix::Linfnorm () const
{
  RealArray sums(nrow,Real(0));

  if (editable)
    for (ValueIter it = elem.begin(); it != elem.end(); ++it)
      sums[it->first.first-1] += fabs(it->second);
  else if (solver == SUPERLU || solver == UMFPACK)
    // Column-oriented format with 0-based row-indices
    for (size_t j = 1; j <= ncol; j++)
      for (int i = IA[j-1]; i < IA[j]; i++)
        sums[JA[i]] += fabs(A[i]);
  else
    // Row-oriented format with 1-based row-indices
    for (size_t i = 1; i <= nrow; i++)
      for (int j = IA[i-1]; j < IA[i]; j++)
        sums[i-1] += fabs(A[j-1]);

  return *std::max_element(sums.begin(),sums.end());
}


void SparseMatrix::calcCSR(IntVec& IA, IntVec& JA,
                           size_t nrow, const ValueMap& elem)
{
  size_t nnz = elem.size();
  IA.resize(nrow+1,nnz);
  JA.resize(nnz);

  IA[0] = 0;
  size_t cur_row = 1, ix = 0;
  for (auto it = elem.begin(); it != elem.end(); ++it, ix++) {
    JA[ix] = it->first.second-1;
    while (it->first.first > cur_row)
      IA[cur_row++] = ix;
  }
}
