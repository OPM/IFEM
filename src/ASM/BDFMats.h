// $Id$
//==============================================================================
//!
//! \file BDFMats.h
//!
//! \date Nov 2 2013
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Representation of the element matrices for a dynamic FEM problem.
//!
//==============================================================================

#ifndef _BDF_MATS_H
#define _BDF_MATS_H

#include "NewmarkMats.h"
#include "BDF.h"


/*!
  \brief Class representing the element matrices for a dynamic FEM problem
  using backward difference formulae (BDF).
*/

class BDFMats : public NewmarkMats
{
public:
  //! \brief The constructor initializes the time integration parameters.
  //! param[in] bdfscheme BDF time discretization scheme
  explicit BDFMats(const TimeIntegration::BDFD2& bdfscheme) : bdf(bdfscheme) {}
  //! \brief Empty destructor.
  virtual ~BDFMats() {}

  //! \brief Returns the element-level Newton matrix.
  virtual const Matrix& getNewtonMatrix() const;
  //! \brief Returns the element-level right-hand-side vector.
  virtual const Vector& getRHSVector() const;

protected:
  TimeIntegration::BDFD2 bdf; //!< BDF time integration scheme
};

#endif
