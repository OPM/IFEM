// $Id$
//==============================================================================
//!
//! \file ReactionsOnly.h
//!
//! \date Nov 7 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Global integral sub-class for reaction force integration.
//!
//==============================================================================

#ifndef _REACTIONS_ONLY_H
#define _REACTIONS_ONLY_H

#include "GlobalIntegral.h"
#include "SystemMatrix.h"


/*!
  \brief Class for assembly of reaction forces only.

  \details This class can be used for linear problems which normally consists
  of a single assembly loop. For non-linear problems, the reaction forces for
  a given solution state is calculated during the subsequent iteration where
  the solution vector of the previous iteration is used. In linear problems we
  therefore need a separate assembly loop where only the reaction forces are
  calculated. This class is provided to facilitate such calculations.
*/

class ReactionsOnly : public GlobalIntegral
{
public:
  //! \brief The constructor initializes the data members.
  //! \param[in] rf The reaction force vector to be assembled
  //! \param[in] sam Data for FE assembly management
  //! \param[in] adm Parallell processing administrator
  //! \param[in] sf Internal force vector to be assembled
  ReactionsOnly(Vector& rf, const SAM* sam, const ProcessAdm& adm,
                Vector* sf = nullptr);
  //! \brief Empty destructor.
  virtual ~ReactionsOnly() {}

  //! \brief Initializes the integrated quantity to zero.
  virtual void initialize(bool);
  //! \brief Finalizes the integrated quantity after element assembly.
  virtual bool finalize(bool);

  //! \brief Adds a LocalIntegral object into a corresponding global object.
  //! \param[in] elmObj The local integral object to add into \a *this.
  //! \param[in] elmId Global number of the element associated with elmObj
  virtual bool assemble(const LocalIntegral* elmObj, int elmId);

  //! \brief Returns \e false if no contributions from the specified patch.
  //! \param[in] pidx 1-based patch index
  //! \param[in] pvec Physical property mapping
  virtual bool haveContributions(size_t pidx,
                                 const std::vector<Property>& pvec) const;

private:
  const SAM*        mySam; //!< Data for FE assembly management
  const ProcessAdm& myAdm; //!< Parallel processing administrator

  Vector&   R; //!< Nodal reaction forces
  StdVector b; //!< Dummy right-hand-side vector
  Vector*   S; //!< Internal force vector
};

#endif
