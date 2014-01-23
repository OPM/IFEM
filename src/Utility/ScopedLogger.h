//==============================================================================
//!
//! \file ScopedLogger.h
//!
//! \date Sep 20 2013
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Scoped logger.
//!
//==============================================================================
#ifndef SCOPED_LOGGER_H_
#define SCOPED_LOGGER_H_

#include <iostream>

/*! \brief Scoped logging class used for logging function entry/departures
*/

class ScopedLogger {
  public:
    //! \brief Constructor
    //! \param name_ The name of the function
    //! \param _stream The stream to log to
    ScopedLogger(const char* name_, std::ostream& _stream=std::cout);

    //! \brief Destructor
    ~ScopedLogger();
  protected:
    const char* name;     //!< The name of the function
    std::ostream& stream; //!< The stream to log to
    int rank;             //!< Process rank (MPI)
};
#endif
