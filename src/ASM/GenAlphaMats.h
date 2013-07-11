// $Id$
//==============================================================================
//!
//! \file GenAlphaMats.h
//!
//! \date Jul 04 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of the element matrices for a dynamic FEM problem.
//!
//==============================================================================

#ifndef _GEN_ALPHA_MATS_H
#define _GEN_ALPHA_MATS_H

#include "NewmarkMats.h"


/*!
  \brief Class representing the element matrices for a dynamic FEM problem
  based on generalized alpha time integration.
*/

class GenAlphaMats : public NewmarkMats
{
public:
  //! \brief The constructor initializes the time integration parameters.
  GenAlphaMats(double alpha, double a, double b);
  //! \brief Empty destructor.
  virtual ~GenAlphaMats() {}

  //! \brief Returns the element-level Newton matrix.
  virtual const Matrix& getNewtonMatrix() const;
  //! \brief Returns the element-level right-hand-side vector.
  virtual const Vector& getRHSVector() const;
};

#endif
