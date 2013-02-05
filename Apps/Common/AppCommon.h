//==============================================================================
//!
//! \file AppCommon.h
//!
//! \date Nov 06 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Common helper templates for applications
//!
//==============================================================================

#ifndef APP_COMMON_H_
#define APP_COMMON_H_


#include "HDF5Writer.h"
#include "XMLWriter.h"

namespace SIM {

//! \brief Handle application restarts
//! \param[template] Simulator The top SIM of your application
//! \param[template] Solver The SIMSolver of your application
//! \param[in] Simulator The top SIM instance of your application
//! \param[in] Simulator The SIMSolver instance of your application
//! \param[in] interval The stride in the input file
//! \param[in] steps The number of timesteps to load
//! \param[in] restartfile The file to read from
  template<class Simulator, class Solver>
void handleRestart(Simulator& simulator, Solver& solver,
                   int interval, int steps, const char* restartfile)
{ 
  DataExporter reader(true, interval, steps);
  simulator.registerFields(reader);
  XMLWriter* xml = new XMLWriter(restartfile);
  HDF5Writer* hdf = new HDF5Writer(restartfile,true);
  reader.registerWriter(xml);
  reader.registerWriter(hdf);
  int max = reader.getTimeLevel();
  if (steps == 2 && reader.getTimeLevel() > -1) {
    reader.loadTimeLevel(reader.getTimeLevel()-1,xml,hdf);
    solver.advanceStep();
  }
  if (max > -1 && reader.loadTimeLevel(-1,xml,hdf)) {
    xml->writeTimeInfo(0, interval, steps, solver.getTimePrm());
    solver.fastForward(reader.realTimeLevel(max));
  }
}


}
#endif
