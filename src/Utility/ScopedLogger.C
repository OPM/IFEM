//==============================================================================
//!
//! \file ScopedLogger.C
//!
//! \date Sep 20 2013
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Scoped logger implementations.
//!
//==============================================================================
#include "ScopedLogger.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

ScopedLogger::ScopedLogger(const char* name_, std::ostream& stream_) :
  name(name_), stream(stream_)
{
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  stream << "[" << rank << "]: Entering \"" << name << "\"" << std::endl;
#else
  rank = -1;
  stream << "Entering \"" << name << "\"" << std::endl;
#endif
}

ScopedLogger::~ScopedLogger()
{
  if (rank == -1)
    stream << "Exiting \"" << name << "\"" << std::endl;
  else
    stream << "[" << rank << "]: Exiting \"" << name << "\"" << std::endl;
}
