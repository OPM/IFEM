// $Id$
//==============================================================================
//!
//! \file SystemMatrix.C
//!
//! \date Jan 4 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief General representation of the system matrix on different forms.
//!
//==============================================================================

#include "SystemMatrix.h"
#include "DenseMatrix.h"
#include "SPRMatrix.h"
#include "SparseMatrix.h"
#include "DiagMatrix.h"
#ifdef HAS_PETSC
#include "PETScMatrix.h"
#endif
#ifdef HAS_ISTL
#include "ISTLMatrix.h"
#endif
#include "LinSolParams.h"
#include <fstream>


SystemVector* SystemVector::create (const ProcessAdm* adm,
                                    LinAlg::MatrixType vType)
{
  switch (vType)
    {
    case LinAlg::ISTL:
#ifdef HAS_ISTL
      if (adm) return new ISTLVector(*adm);
#endif
      break;

    case LinAlg::PETSC:
#ifdef HAS_PETSC
      if (adm) return new PETScVector(*adm);
#endif
      break;

    default:
      return new StdVector();
    }

  std::cerr <<"SystemVector::create: Unsupported vector type "
            << vType << std::endl;
  return nullptr;
}


SystemVector& SystemVector::copy (const SystemVector& x)
{
  this->redim(x.dim());
  memcpy(this->getPtr(),x.getRef(),x.dim()*sizeof(Real));
  return *this;
}


void StdVector::assemble (const Vectors& vecs,
                          const IntVec& meqn,
                          const int neq)
{
  size_t ofs = 0;
  for (const Vector& v : vecs) {
    for (size_t i = 0; i < v.size(); ++i)
      (*this)(ofs+meqn[i]+1) += v[i];
    ofs += neq;
  }
}


void StdVector::dump (const utl::vector<Real>& x, const char* label,
                      LinAlg::StorageFormat format, std::ostream& os)
{
  switch (format)
    {
    case LinAlg::MATLAB:
      utl::writeMatlab(label,x,os);
      break;

    case LinAlg::MATRIX_MARKET:
      os <<"%%MatrixMarket matrix array real general\n"<< x.size() <<" 1";
      for (Real v : x) os <<"\n"<< v;
      os << std::endl;
      break;

    case LinAlg::FLAT:
      if (label)
        os << label <<" =";
      os << x;
      break;

    default:
      break;
    }
}


SystemMatrix* SystemMatrix::create (const ProcessAdm* adm,
                                    LinAlg::MatrixType mType,
                                    const LinSolParams& spar)
{
#ifdef HAS_PETSC
  if (mType == LinAlg::PETSC && adm)
    return new PETScMatrix(*adm,spar);
#endif
#ifdef HAS_ISTL
  if (mType == LinAlg::ISTL && adm)
    return new ISTLMatrix(*adm,spar);
#endif

  return SystemMatrix::create(adm,mType);
}


SystemMatrix* SystemMatrix::create (const ProcessAdm* adm,
                                    LinAlg::MatrixType mType,
                                    int num_thread_SLU)
{
#if defined(HAS_PETSC) || defined(HAS_ISTL)
  // Use default settings when no parameters are provided by user
  static LinSolParams defaultPar;
#endif

  switch (mType)
    {
    case LinAlg::DENSE:
      return new DenseMatrix();

    case LinAlg::SPR:
      return new SPRMatrix();

    case LinAlg::SPARSE:
      return new SparseMatrix(SparseMatrix::SUPERLU,num_thread_SLU);

    case LinAlg::SAMG:
      return new SparseMatrix(SparseMatrix::S_A_M_G);

    case LinAlg::PETSC:
#ifdef HAS_PETSC
      if (adm) return new PETScMatrix(*adm,defaultPar);
#else
      std::cerr <<"SystemMatrix::create: PETSc not compiled in."<< std::endl;
#endif
      break;

    case LinAlg::ISTL:
#ifdef HAS_ISTL
      if (adm) return new ISTLMatrix(*adm,defaultPar);
#else
      std::cerr <<"SystemMatrix::create: ISTL not compiled in."<< std::endl;
#endif
      break;

    case LinAlg::UMFPACK:
      return new SparseMatrix(SparseMatrix::UMFPACK);

    case LinAlg::DIAG:
      return new DiagMatrix();

    default:
      break;
    }

  std::cerr <<"SystemMatrix::create: Unsupported matrix type "
            << mType << std::endl;
  return nullptr;
}


void SystemMatrix::initNonZeroEqs ()
{
  nonZeroEqs.clear();
  nonZeroEqs.resize(this->dim(1),false);
}


bool SystemMatrix::flagNonZeroEqs (const IntVec& meq)
{
  if (meq.empty() || nonZeroEqs.size() == 1)
    std::fill(nonZeroEqs.begin(),nonZeroEqs.end(),true);
  else
  {
    for (int ieq : meq)
      if (ieq < 0 || ieq > static_cast<int>(nonZeroEqs.size()))
        return false; //TODO: for ieq < 0, find the actual master dofs
      else if (ieq > 0)
        nonZeroEqs[ieq-1] = true;
  }

  return true;
}


bool SystemMatrix::flagNonZeroEqs (const SystemMatrix& B)
{
  if (nonZeroEqs.size() == 1)
    nonZeroEqs.front() = true;
  else
  {
    for (size_t i = 0; i < nonZeroEqs.size(); i++)
      if (i >= B.nonZeroEqs.size())
        break;
      else if (B.nonZeroEqs[i])
        nonZeroEqs[i] = true;
  }

  return true;
}


bool SystemMatrix::isZero () const
{
  if (this->dim(1) < 1)
    return false;

  return std::none_of(nonZeroEqs.begin(),nonZeroEqs.end(),
                      [](bool nz){ return nz; });
}


/*!
  This method will insert 1.0e9 on the diagonal, for all equations without
  contributions, such that the remaining equations can be solved for,
  if the method initNonZeroEqs() was invoked before the assembly started.
*/

bool SystemMatrix::endAssembly ()
{
  if (nonZeroEqs.size() > 1 && !this->isZero())
    for (size_t ieq = 1; ieq <= nonZeroEqs.size(); ieq++)
      if (!nonZeroEqs[ieq-1] && !this->add(1.0e9,ieq))
        return false;

  return true;
}


StdVector SystemMatrix::operator* (const SystemVector& b) const
{
  StdVector results;
  this->multiply(b,results);
  return results;
}


StdVector SystemMatrix::operator/ (const SystemVector& b)
{
  StdVector results;
  this->solve(b,results);
  return results;
}


void SystemMatrix::dump (const char* fileName, std::streamsize precision,
                         LinAlg::StorageFormat format)
{
  if (format == LinAlg::BINARY)
  {
    std::ofstream fs(fileName,std::ofstream::binary);
    this->dump(fs,format);
  }
  else
  {
    std::ofstream fs(fileName);
    if (precision > 0)
      fs.precision(precision);
    this->dump(fs,format);
  }
}
