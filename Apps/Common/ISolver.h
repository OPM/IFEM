// $Id$
//==============================================================================
//!
//! \file ISolver.h
//!
//! \date Oct 14 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Abstract simulation solver interface.
//!
//==============================================================================

#ifndef _I_SOLVER_H_
#define _I_SOLVER_H_

class TimeStep;


/*!
  \brief Abstract solver interface.
  \details To be realized through the SIMSolver template.
*/

class ISolver
{
public:
  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  virtual bool saveModel(char* fileName, int& geoBlk, int& nBlock) = 0;

  //! \brief Saves the converged results of a given time step to VTF file.
  //! \param[in] tp Time stepping parameters
  //! \param[in] nBlock Running VTF block counter
  virtual bool saveStep(const TimeStep& tp, int& nBlock) = 0;

  //! \brief Initializes the simulator time stepping loop.
  virtual bool init(const TimeStep& tp) = 0;

  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep& tp) = 0;

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp) = 0;

  //! \brief Registers fields for output to a data exporter.
  virtual void registerFields(DataExporter& exporter) = 0;

  //! \brief Sets initial conditions for the solution field.
  virtual void setInitialConditions() = 0;

  //! \brief Initializes and sets up field dependencies.
  virtual void setupDependencies() = 0;
};

#endif
