// $Id$
//==============================================================================
//!
//! \file AndersonAcceleration.h
//!
//! \date Sep 18 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Anderson acceleration of a fixed point iteration.
//!
//==============================================================================

#ifndef _ANDERSON_ACCELERATION_H
#define _ANDERSON_ACCELERATION_H

#include "MatVec.h"

#include <cstddef>

class SAM;


/*!
  \brief Anderson acceleration of a fixed point iteration.

  \details A fixed point iteration \a x <- \a G(x) converges at a rate set by
  the spectrum of the iteration, and a scheme that is a relaxed Richardson
  iteration underneath, as a projection scheme iterated towards a stationary
  state is, converges slowly. Anderson acceleration replaces the next point
  by the combination of the last few arrival points whose residuals cancel
  best in the least squares sense.

  The vectors combined are nodal solution vectors, so every inner product
  taken here goes through the SAM of the model. A domain decomposed run
  would otherwise give each process its own coefficients, and the combined
  point would not be continuous across the interfaces.
*/

class AndersonAcceleration
{
public:
  //! \brief The constructor sets the parameters of the acceleration.
  //! \param[in] mDepth Number of previous steps to combine, zero to disable
  //! \param[in] reg Tikhonov term, relative to the trace of the system
  //! \param[in] growth Residual growth that drops the history
  explicit AndersonAcceleration(size_t mDepth = 10,
                                double reg = 1.0e-10,
                                double growth = 1.0e4)
    : depth(mDepth), eps(reg), maxGrowth(growth) {}

  //! \brief Returns the number of previous steps combined.
  size_t getDepth() const { return depth; }

  //! \brief Sets the number of previous steps to combine.
  void setDepth(size_t mDepth) { depth = mDepth; }

  //! \brief Forgets the history.
  //! \details To be called whenever the points stop being comparable, which
  //! a change of the mesh makes them.
  void clear();

  //! \brief Accelerates one step of a fixed point iteration.
  //! \param[in] sam Assembly management data, for the inner products
  //! \param[in] x The point the step was taken from
  //! \param gx The point the step arrived at, replaced by the accelerated one
  //! \return \e false if the least squares problem could not be solved
  bool accelerate(const SAM& sam, const Vector& x, Vector& gx);

  //! \brief Returns the number of steps currently in the history.
  size_t getSize() const { return dF.size(); }

private:
  size_t depth;     //!< Number of previous steps to combine
  double eps;       //!< Tikhonov term, relative to the trace of the system
  double maxGrowth; //!< Residual growth that drops the history

  Vectors dF; //!< Differences of the residuals of the steps in the history
  Vectors dG; //!< Differences of the arrival points of those steps
  Vector prevF; //!< Residual of the previous step
  Vector prevG; //!< Arrival point of the previous step
  double bestF = -1.0; //!< Smallest residual norm seen since the last drop
};

#endif
