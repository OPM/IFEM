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
#include <dune/common/version.hh>
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
  using Mat = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>; //!< A sparse system matrix
  using Vec = Dune::BlockVector<Dune::FieldVector<double,1>>;  //!< A vector
  using Operator = Dune::MatrixAdapter<Mat,Vec,Vec>;           //!< A serial matrix operator
  using InverseOperator = Dune::InverseOperator<Vec, Vec>;     //!< Linear system inversion abstraction
  using Preconditioner = Dune::Preconditioner<Vec,Vec>;        //!< Preconditioner abstraction

  /*! \brief Wrapper template to avoid memory leaks */
  template<template<class M> class Pre>
  class IOp2Pre : public Dune::InverseOperator2Preconditioner<Pre<ISTL::Mat>,
                                                              Dune::SolverCategory::sequential>
  {
    using SolverType = Dune::InverseOperator2Preconditioner<Pre<ISTL::Mat>,
                                                            Dune::SolverCategory::sequential>;

  public:
    explicit IOp2Pre(Pre<ISTL::Mat>* iop) : SolverType(*iop)
    {
      m_op.reset(iop);
    }

  protected:
    std::unique_ptr<Pre<ISTL::Mat>> m_op;
  };

#if defined(HAVE_UMFPACK)
  using LUType = Dune::UMFPack<ISTL::Mat>;
  using LU = IOp2Pre<Dune::UMFPack>;
#elif defined(HAVE_SUPERLU)
  using LUType = Dune::SuperLU<ISTL::Mat>;
  using LU = IOp2Pre<Dune::SuperLU>;
#endif

  // Preconditioner types
  using ILU = Dune::SeqILU<Mat, Vec, Vec>;  //!< Sequential ILU
  using SOR = Dune::SeqSOR<Mat, Vec, Vec>;  //!< Sequential SOR
  using SSOR = Dune::SeqSSOR<Mat, Vec,Vec>; //!< Sequential SSOR
  using GJ = Dune::SeqJac<Mat, Vec, Vec>;   //!< Sequential Gauss-Jacobi
  using GS = Dune::SeqGS<Mat, Vec, Vec>;    //!< Sequential Gauss-Seidel
  using ASMLU = Dune::SeqOverlappingSchwarz<Mat, Vec,
                                            Dune::AdditiveSchwarzMode,
                                            LUType>; //!< Sequential additive Schwarz w/ LU subdomain solvers
  using ASM = Dune::SeqOverlappingSchwarz<Mat, Vec>; //!< Sequential additive Schwarz w/ ILU subdomain solvers
}
#endif

#endif
