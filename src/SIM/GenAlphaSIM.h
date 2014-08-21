// $Id$
//==============================================================================
//!
//! \file GenAlphaSIM.h
//!
//! \date Aug 21 2014
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Generalized-alpha driver for isogeometric dynamic FEM simulators.
//!
//==============================================================================

#ifndef _GEN_ALPHA_SIM_H
#define _GEN_ALPHA_SIM_H

#include "NewmarkSIM.h"


/*!
  \brief Generalized-alpha driver for dynamic isogeometric FEM simulators.
*/

class GenAlphaSIM : public NewmarkSIM
{
public:
  //! \brief The constructor initializes default solution parameters.
  GenAlphaSIM(SIMbase& s);
  //! \brief Empty destructor.
  virtual ~GenAlphaSIM() {}

  //! \brief Parses a data section from an XML document.
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Initializes time integration parameters for the integrand.
  virtual void initPrm();
  //! \brief Initializes primary solution vectors.
  //! \param[in] nSol Number of consequtive solutions stored
  virtual bool initSol(size_t nSol = 3);

  //! \brief Advances the time step one step forward.
  //! \param param Time stepping parameters
  //! \param[in] updateTime If \e false, the time parameters are not incremented
  virtual bool advanceStep(TimeStep& param, bool updateTime = true);

  //! \brief Modifies the current solution vector (used by sub-iterations only).
  virtual void setSolution(const Vector& newSol, int idx);

protected:
  //! \brief Calculates predicted velocities and accelerations.
  virtual bool predictStep(TimeStep& param);
  //! \brief Updates configuration variables (solution vector) in an iteration.
  virtual bool correctStep(TimeStep& param, bool converged);

private:
  double alpha_m;  //!< Generalized-alpha parameter
  double alpha_f;  //!< Generalized-alpha parameter
  Vectors prevSol; //!< Solution of last converged time step
  Vectors tempSol; //!< Intermediate solution
};

#endif
