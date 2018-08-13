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
#include "MatVec.h"
#include "ASMenums.h"


/*!
  \brief Class representing integrated global norms.
  \details The class is essentially a vector of doubles, but is derived from
  GlobalIntegral such that it may be passed as argument to ASMbase::integrate.
*/

class GlbNorm : public GlobalIntegral
{
public:
  //! \brief The constructor initializes a reference to the global norm vector.
  //! \param[in] v Vector of global norm quantities
  //! \param[in] op Operation to be performed after accumulating element norms
  GlbNorm(Vectors& v, ASM::FinalNormOp op = ASM::NONE);
  //! \brief The destructor applies the operation \a myOp on \a myVals.
  virtual ~GlbNorm();

  //! \brief Adds element norm quantities into the global norm object.
  //! \param[in] elmObj Pointer to the element norms to add into \a *this
  //! \param[in] elmId Global number of the element associated with \a *elmObj
  virtual bool assemble(const LocalIntegral* elmObj, int elmId = 0);

  //! \brief Returns \e true if all elements can be assembled in parallel.
  virtual bool threadSafe() const { return delAss; }
  //! \brief Flags that the element values are assembled in the final stage.
  void delayAssembly() { delAss = true; }

private:
  //! \brief Applies the operation \a myOp on the given \a value.
  void applyFinalOp(double& value) const;

  Vectors&         myVals; //!< Reference to a vector of global norm values
  ASM::FinalNormOp myOp;   //!< Operation to be performed on summed values
  bool             delAss; //!< If \e true, element assembly is delayed
};

#endif
