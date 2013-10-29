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

class ScopedLogger {
  public:
    ScopedLogger(const char* name_, std::ostream& _stream=std::cout);

    ~ScopedLogger();
  protected:
    const char* name;
    std::ostream& stream;
    int rank;
};
#endif
