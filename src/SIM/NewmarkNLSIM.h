// $Id$
//==============================================================================
//!
//! \file NewmarkNLSIM.h
//!
//! \date Jul 4 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Nonlinear solution driver for isogeometric dynamic FEM simulators.
//!
//==============================================================================

#ifndef _NEWMARK_NL_SIM_H
#define _NEWMARK_NL_SIM_H

#include "NewmarkSIM.h"

class SystemVector;


/*!
  \brief Newmark-based solution driver for dynamic isogeometric FEM simulators.
*/

class NewmarkNLSIM : public NewmarkSIM
{
public:
  //! \brief The constructor initializes default solution parameters.
  //! \param sim The FE model
  NewmarkNLSIM(SIMbase& sim);
  //! \brief Empty destructor.
  virtual ~NewmarkNLSIM() {}

  //! \brief Parses a data section from an XML document.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Initializes primary solution vectors and integration parameters.
  //! \param[in] nSol Number of consequtive solutions stored
  virtual void init(size_t nSol = 3);

  //! \brief Advances the time step one step forward.
  //! \param param Time stepping parameters
  //! \param[in] updateTime If \e false, the time parameters are not incremented
  virtual bool advanceStep(TimeStep& param, bool updateTime = true);

protected:
  //! \brief Calculates predicted velocities and accelerations.
  virtual bool predictStep(TimeStep& param);
  //! \brief Updates configuration variables (solution vector) in an iteration.
  virtual bool correctStep(TimeStep& param, bool = false);
  //! \brief Finalizes the right-hand-side vector on the system level.
  virtual void finalizeRHSvector();

private:
  SystemVector* Finert; //!< Actual inertia forces in last converged time step
};

#endif
