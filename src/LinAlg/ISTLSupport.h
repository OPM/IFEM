// $Id$
//==============================================================================
//!
//! \file ISTLSupport.h
//!
//! \date May 12 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief IFEM ISTL support
//!
//==============================================================================

#ifndef _ISTL_SUPPORT_H_
#define _ISTL_SUPPORT_H_

#include <vector>

#ifdef HAS_ISTL
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/overlappingschwarz.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/owneroverlapcopy.hh>
#ifdef HAVE_SUPERLU
#include <dune/istl/superlu.hh>
#endif
#ifdef HAVE_UMFPACK
#include <dune/istl/umfpack.hh>
#endif

namespace ISTL
{
  typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> Mat; //!< A sparse system matrix
  typedef Dune::BlockVector<Dune::FieldVector<double,1>> Vec;  //!< A vector
  typedef Dune::MatrixAdapter<Mat,Vec,Vec> Operator;      //!< A serial matrix operator
  typedef Dune::InverseOperator<Vec, Vec> InverseOperator;     //!< Linear system inversion abstraction
  typedef Dune::Preconditioner<Vec,Vec> Preconditioner;        //!< Preconditioner abstraction

  /*! \brief Wrapper template to avoid memory leaks */
  template<template<class M> class Pre>
  class IOp2Pre : public Dune::InverseOperator2Preconditioner<Pre<ISTL::Mat>, Dune::SolverCategory::sequential> {
    typedef Dune::InverseOperator2Preconditioner<Pre<ISTL::Mat>, Dune::SolverCategory::sequential> SolverType;
  public:
    IOp2Pre(Pre<ISTL::Mat>* iop) : SolverType(*iop)
    {
      m_op.reset(iop);
    }
  protected:
    std::unique_ptr<Pre<ISTL::Mat>> m_op;
  };

#if defined(HAVE_UMFPACK)
  typedef Dune::UMFPack<ISTL::Mat> LUType;
  typedef IOp2Pre<Dune::UMFPack> LU;
#elif defined(HAVE_SUPERLU)
  typedef Dune::SuperLU<ISTL::Mat> LUType;
  typedef IOp2Pre<Dune::SuperLU> LU;
#endif

  // Preconditioner types
  typedef Dune::SeqILU0<Mat, Vec, Vec> ILU0; //!< Sequential ILU(0)
  typedef Dune::SeqILUn<Mat, Vec, Vec> ILUn; //!< Sequential ILU(n)
  typedef  Dune::SeqSOR<Mat, Vec, Vec>  SOR; //!< Sequential SOR
  typedef Dune::SeqSSOR<Mat, Vec,Vec> SSOR; //!< Sequential SSOR
  typedef  Dune::SeqJac<Mat, Vec, Vec>    GJ; //!< Sequential Gauss-Jacobi
  typedef   Dune::SeqGS<Mat, Vec, Vec>     GS; //!< Sequential Gauss-Seidel
  typedef Dune::SeqOverlappingSchwarz<Mat, Vec,
                                 Dune::AdditiveSchwarzMode,
                                 LUType> ASMLU; //!< Sequential additive Schwarz w/ LU subdomain solvers
  typedef Dune::SeqOverlappingSchwarz<Mat, Vec> ASM; //!< Sequential additive Schwarz w/ ILU subdomain solvers
}
#endif

#endif
