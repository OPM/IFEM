// $Id$
//==============================================================================
//!
//! \file ASMs2DIB.h
//!
//! \date Dec 18 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of structured 2D spline FE models with immersed boundaries.
//!
//==============================================================================

#ifndef _ASM_S2D_IB_H
#define _ASM_S2D_IB_H

#include "ASMs2D.h"


/*!
  \brief Assembly of structured 2D spline FE models with immersed boundaries.
*/

class ASMs2DIB : public ASMs2D
{
public:
  //! \brief Default constructor.
  ASMs2DIB(unsigned char n_s = 2, unsigned char n_f = 1, int max_depth = 5);
  //! \brief Copy constructor.
  ASMs2DIB(const ASMs2DIB& patch, unsigned char n_f = 0);
  //! \brief Empty destructor.
  virtual ~ASMs2DIB() {}

  //! \brief Generates the finite element topology data for the patch.
  //! \details This method is overridden in this class, to include the
  //! calculation of the quadrature point parameters (and weights)
  //! according to the Immersed boundary scheme.
  virtual bool generateFEMTopology();

  //! \brief Evaluates an integral over the interior patch domain.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrate(Integrand& integrand,
                         GlobalIntegral& glbInt, const TimeDomain& time);

private:
  Real3DMat quadPoints; //!< The Gauss quadrature points for this patch
  int       maxDepth;   //!< Maximum depth up to which to refine each element
};

#endif
