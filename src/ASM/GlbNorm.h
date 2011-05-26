// $Id$
//==============================================================================
//!
//! \file GlbNorm.h
//!
//! \date Dec 09 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of integrated global norm quantities.
//!
//==============================================================================

#ifndef _GLB_NORM_H
#define _GLB_NORM_H

#include "GlobalIntegral.h"
#include <vector>


/*!
  \brief Class representing integrated global norms.
  \details The class is essentially a vector of doubles, but is derived from
  GlobalIntegral such that it may be passed as argument to ASMbase::integrate.
*/

class GlbNorm : public GlobalIntegral
{
public:
  //! \brief Operations to applied after summing element contributions.
  enum FinalOp { NONE, ABS, SQRT };

  //! \brief The constructor initializes a reference to the global norm vector.
  //! \param[in] vec Vector of global norm quantities
  //! \param[in] op Operation to be performed after accumulating element norms
  GlbNorm(std::vector<double>& vec, FinalOp op = NONE);
  //! \brief The destructor applies the operation \a myOp on \a myVals
  virtual ~GlbNorm();

  //! \brief Adds element norm quantities into the global norm object.
  //! \param[in] elmObj Pointer to the element norms to add into \a *this
  //! \param[in] elmId Global number of the element associated with \a *elmObj
  virtual bool assemble(const LocalIntegral* elmObj, int elmId);

private:
  //! \brief Applies the operation \a myOp on the given \a value.
  void applyFinalOp(double& value) const;

  std::vector<double>& myVals; //!< Reference to a vector of global norm values
  FinalOp              myOp;   //!< Operation to be perfomed on summed values
};

#endif
