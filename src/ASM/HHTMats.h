// $Id$
//==============================================================================
//!
//! \file HHTMats.h
//!
//! \date Nov 13 2014
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of the element matrices for a dynamic FEM problem.
//!
//==============================================================================

#ifndef _HHT_MATS_H
#define _HHT_MATS_H

#include "NewmarkMats.h"


/*!
  \brief Class representing the element matrices for a dynamic FEM problem
  based on generalized alpha time integration.
*/

class HHTMats : public NewmarkMats
{
public:
  //! \brief The constructor initializes the time integration parameters.
  HHTMats(double alpha, double a, double b, bool old = false);
  //! \brief Empty destructor.
  virtual ~HHTMats() {}

  //! \brief Returns the element-level Newton matrix.
  virtual const Matrix& getNewtonMatrix() const;
  //! \brief Returns the element-level right-hand-side vector.
  virtual const Vector& getRHSVector() const;

private:
  bool oldHHT; //!< If \e true, used toghether with NewmarkNLSIM
};

#endif
