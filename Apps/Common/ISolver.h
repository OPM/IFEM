//==============================================================================
//!
//! \file ISolver.h
//!
//! \date Oct 14 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Abstract simulation solver interface
//!
//==============================================================================

#include "TimeStep.h"


/*!
  \brief Abstract solver interface. To be realized through the SIMSolver template
*/
class ISolver {
  public:
    //! \brief Opens a new VTF-file and writes the model geometry to it.
    //! \param[in] fileName File name used to construct the VTF-file name from
    virtual bool saveModel(char* fileName, int& nBlock) = 0;

    //! \brief Saves the converged results to VTF file of a given time step.
    //! \param[in] iStep Time step identifier
    //! \param[in] time Current time step info
    //! \param[in] nBlock Running VTF block counter
    virtual bool saveStep(const TimeStep& tp, int& nBlock) = 0;

    //! \brief Advances the time step one step forward.
    //! \param[in] tp Time step structure to advance
    //! \return True if new solution step is to be performed
    virtual bool advanceStep(TimeStep& tp) = 0;

    //! \brief Computes the solution for the current time step.
    //! \param[in] tp Time step structure to advance
    //! \return True on success
    virtual bool solveStep(TimeStep& tp) = 0;

    //! \brief Returns a const reference to the time stepping parameters.
    const TimeStep& getTimePrm() const = 0;
};
