// $Id$
//==============================================================================
//!
//! \file GlbL2projector.h
//!
//! \date Oct 16 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief General integrand for L2-projection of secondary solutions.
//!
//==============================================================================

#ifndef _GLB_L2_PROJECTOR_H
#define _GLB_L2_PROJECTOR_H

#include "IntegrandBase.h"
#include "SparseMatrix.h"


/*!
  \brief General integrand for L2-projection of secondary solutions.
*/

class GlbL2 : public Integrand
{
public:
  //! \brief The constructor initializes the projection matrices.
  //! \param[in] p The main problem integrand
  //! \param[in] n Dimension of the L2-projection matrices (number of nodes)
  GlbL2(IntegrandBase& p, size_t n);
  //! \brief Empty destructor.
  virtual ~GlbL2() {}

  //! \brief Defines which FE quantities are needed by the integrand.
  virtual int getIntegrandType() const { return problem.getIntegrandType(); }

  //! \brief Returns a local integral contribution object for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] iEl Global element number (1-based)
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t iEl,
                                          bool neumann) const;

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] fe Nodal and integration point data for current element
  //! \param[in] X0 Cartesian coordinates of the element center
  //! \param[in] nPt Number of integration points in this element
  //! \param elmInt Local integral for element
  virtual bool initElement(const IntVec& MNPC, const FiniteElement& fe,
                           const Vec3& X0, size_t nPt,
                           LocalIntegral& elmInt);

  //! \brief Dummy implementation.
  virtual bool initElement(const IntVec&, LocalIntegral&) { return false; }
  //! \brief Dummy implementation.
  virtual bool initElement(const IntVec&, const IntVec&, size_t,
                           LocalIntegral&) { return false; }
  //! \brief Dummy implementation.
  virtual bool initElementBou(const IntVec&, LocalIntegral&) { return false; }
  //! \brief Dummy implementation.
  virtual bool initElementBou(const IntVec&, const IntVec&, size_t,
                              LocalIntegral&) { return false; }

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt,
                       const FiniteElement& fe, const Vec3& X) const;

  //! \brief Pre-computes the sparsity pattern of the projection matrix \b A.
  //! \param[in] MMNPC Matrix of matrices of nodal point correspondances
  //! \param[in] nel Number of elements
  void preAssemble(const std::vector<IntVec>& MMNPC, size_t nel);

  //! \brief Solves the projection equation system and evaluates nodal values.
  //! \param[out] sField Nodal/control-point values of the projected results.
  bool solve(Matrix& sField);

private:
  IntegrandBase& problem; //!< The main problem integrand
  mutable SparseMatrix A; //!< Left-hand-side matrix of the L2-projection
  mutable StdVector    B; //!< Right-hand-side vectors of the L2-projection
};

#endif
