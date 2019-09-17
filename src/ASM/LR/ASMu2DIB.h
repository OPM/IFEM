// $Id$
//==============================================================================
//!
//! \file ASMu2DIB.h
//!
//! \date Dec 18 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of unstructured 2D spline FE models with immersed boundary.
//!
//==============================================================================

#ifndef _ASM_U2D_IB_H
#define _ASM_U2D_IB_H

#include "ASMu2D.h"

namespace Immersed { class Geometry; }


/*!
  \brief Assembly of unstructured 2D spline FE models with immersed boundaries.
*/

class ASMu2DIB : public ASMu2D
{
public:
  //! \brief Default constructor.
  ASMu2DIB(unsigned char n_s = 2, unsigned char n_f = 1, int max_depth = 5);
  //! \brief Copy constructor.
  ASMu2DIB(const ASMu2DIB& patch, unsigned char n_f = 0);
  //! \brief The destructor deletes the dynamically allocated geometry object.
  virtual ~ASMu2DIB();

  //! \brief Adds a circular hole in the physical geometry.
  //! \param[in] R Hole radius
  //! \param[in] Xc X-coordinate of the hole centre
  //! \param[in] Yc Y-coordinate of the hole centre
  virtual void addHole(double R, double Xc, double Yc);
  //! \brief Adds an oval hole in the physical geometry.
  //! \param[in] R Hole radius
  //! \param[in] X1 X-coordinate of the first circle centre
  //! \param[in] Y1 Y-coordinate of the first circle centre
  //! \param[in] X2 X-coordinate of the second circle centre
  //! \param[in] Y2 Y-coordinate of the second circle centre
  virtual void addHole(double R, double X1, double Y1, double X2, double Y2);

  //! \brief Defines the immersed geometry from a scalar function.
  //! \param[in] f The function to use
  //! \param[in] power Exponent to apply on the specified function
  //! \param[in] threshold Inside/outside threshold
  virtual bool setGeometry(RealFunc* f, double power, double threshold);

  //! \brief Returns an additional geometry to visualize (immersed boundaries).
  virtual ElementBlock* immersedGeometry() const;

  //! \brief Computes the total number of integration points in this patch.
  virtual void getNoIntPoints(size_t& nPt, size_t& nIPt);

  //! \brief Generates the finite element topology data for the patch.
  //! \details This method is overridden in this class, to include the
  //! calculation of the quadrature point parameters (and weights)
  //! according to the Immersed boundary scheme.
  virtual bool generateFEMTopology();

  using ASMu2D::integrate;
  //! \brief Evaluates an integral over the interior patch domain.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrate(Integrand& integrand,
                         GlobalIntegral& glbInt, const TimeDomain& time);

  //! \brief Filters out result point values that are outside physical domain.
  //! \param field Array of result values
  //! \param[in] grid The visualization grid
  virtual void filterResults(Matrix& field, const ElementBlock* grid) const;

private:
  Immersed::Geometry* myGeometry; //!< The physical geometry description

  Real3DMat quadPoints; //!< The Gauss quadrature points for this patch
  int       maxDepth;   //!< Maximum depth up to which to refine each element
};

#endif
