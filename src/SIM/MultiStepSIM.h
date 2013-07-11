// $Id$
//==============================================================================
//!
//! \file MultiStepSIM.h
//!
//! \date Jul 11 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for multi-step solution drivers.
//!
//==============================================================================

#ifndef _MULTI_STEP_SIM_H
#define _MULTI_STEP_SIM_H

#include "SIMinput.h"
#include "MatVec.h"

class SIMbase;


/*!
  \brief Base class for multi-step solution drivers.
*/

class MultiStepSIM : public SIMinput
{
public:
  //! \brief The constructor initializes the FE model reference.
  //! \param sim The FE model
  MultiStepSIM(SIMbase& sim);
  //! \brief Empty destructor.
  virtual ~MultiStepSIM() {}

  //! \brief Returns a list of prioritized XML-tags.
  virtual const char** getPrioritizedTags() const;

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  bool saveModel(char* fileName);

  //! \brief Saves the converged results to VTF file of a given time step.
  //! \param[in] iStep Time/load step identifier
  //! \param[in] time Current time/load parameter
  //! \param[in] psolOnly If \e true, skip secondary solution field output
  //! \param[in] vecName Optional name of primary solution vector field
  bool saveStep(int iStep, double time,
                bool psolOnly = false, const char* vecName = NULL);

  //! \brief Dumps the primary solution for inspection.
  //! \param[in] iStep Time/load step identifier
  //! \param[in] time Current time/load parameter
  //! \param[in] os The output stream to write the solution to
  //! \param[in] withID If \e true, write node ID and coordinates too
  void dumpStep(int iStep, double time, std::ostream& os,
                bool withID = true) const;

  //! \brief Dumps solution variables at user-defined points.
  //! \param[in] time Current time/load parameter
  //! \param[in] os The output stream to write the solution to
  //! \param[in] precision Number of digits after the decimal point
  //! \param[in] formatted If \e false, write all result points on a single line
  void dumpResults(double time, std::ostream& os,
                   std::streamsize precision = 3, bool formatted = true) const;

  //! \brief Returns a const reference to current solution vector.
  const Vector& getSolution(int idx = 0) const { return solution[idx]; }

protected:
  SIMbase& model;    //!< The isogeometric FE model
  Vector   residual; //!< Residual force vector
  Vector   linsol;   //!< Linear solution vector
  Vectors  solution; //!< Primary solution vectors

  int geoBlk; //!< Running VTF geometry block counter
  int nBlock; //!< Running VTF result block counter
};

#endif
