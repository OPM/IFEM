// $Id$
//==============================================================================
//!
//! \file GlbForce.h
//!
//! \date Jun 26 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Representation of integrated global force quantities.
//!
//==============================================================================

#ifndef _GLB_FORCE_H
#define _GLB_FORCE_H

#include "GlobalIntegral.h"
#include "IntegrandBase.h"
#include <vector>
#include "MatVec.h"
#include "ElmNorm.h"


/*!
  \brief Class representing integrated global forces.
  \details The class is essentially a vector of doubles, but is derived from
  GlobalIntegral such that it may be passed as argument to ASMbase::integrate.
*/

class GlbForce : public GlobalIntegral
{
public:
  //! \brief The constructor initializes a reference to the global norm vector.
  //! \param[in] vec Vector of global norm quantities
  //! \param[in] op Operation to be performed after accumulating element norms
  GlbForce(std::vector<double>& vec) : myVals(vec)
  {
  }

  //! \brief Empty destructor
  virtual ~GlbForce() {}

  //! \brief Adds element norm quantities into the global norm object.
  //! \param[in] elmObj Pointer to the element norms to add into \a *this
  //! \param[in] elmId Global number of the element associated with \a *elmObj
  //! \details Dummy stub since forces are added up after integration
  virtual bool assemble(const LocalIntegral* elmObj, int elmId = 0)
  { 
    return true;
  }

  void assemble(const LintegralVec& elms)
  {
    for (size_t i=0;i<elms.size();++i)
      for (size_t j=0;j<myVals.size();++j)
        myVals[j] += (*(static_cast<ElmNorm*>(elms[i])))[j];
  }
private:
  std::vector<double>& myVals; //!< Reference to a vector of global force values
};

#endif
