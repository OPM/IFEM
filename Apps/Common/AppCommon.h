// $Id$
//==============================================================================
//!
//! \file AppCommon.h
//!
//! \date Nov 06 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Common helper templates for applications.
//!
//==============================================================================

#ifndef _APP_COMMON_H_
#define _APP_COMMON_H_

#include "HDF5Writer.h"
#include "XMLWriter.h"


namespace SIM
{
  //! \brief Handles application restarts.
  //! \param[in] simulator The top SIMbase instance of your application
  //! \param[in] solver The SIMSolver instance of your application
  //! \param[in] interval The stride in the input file
  //! \param[in] steps The number of timesteps to load
  //! \param[in] restartfile The file to read from
  template<class Simulator, class Solver>
  void handleRestart(Simulator& simulator, Solver& solver,
                     const std::string& restartfile, int interval, int steps)
  {
    DataExporter reader(true,interval,steps);
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

  //! \brief Handles application restarts.
  //TODO: Remove this later. Obsolete wrapper.
  template<class Simulator, class Solver>
  void handleRestart(Simulator& simulator, Solver& solver,
		     int interval, int steps, const char* restartfile)
  { return handleRestart(simulator,solver,restartfile,interval,steps); }

  //! \brief Handles application data output.
  //! \param[in] simulator The top SIMbase instance of your application
  //! \param[in] solver The SIMSolver instance of your application
  //! \param[in] hdf5file The file to save to
  //! \param[in] append Whether or not to append to file
  //! \param[in] interval The stride in the input file
  //! \param[in] steps The number of timesteps to load
  template<class Simulator, class Solver>
  DataExporter* handleDataOutput(Simulator& simulator, Solver& solver,
                                 const std::string& hdf5file,
                                 bool append, int interval, int steps)
  {
    DataExporter* writer = new DataExporter(true,interval,steps);
    simulator.registerFields(*writer);
    XMLWriter* xml = new XMLWriter(hdf5file);
    HDF5Writer* hdf = new HDF5Writer(hdf5file,append);
    writer->registerWriter(xml);
    writer->registerWriter(hdf);
    if (!append)
      writer->dumpTimeLevel(&solver.getTimePrm()); // initial state

    return writer;
  }

  //! \brief Handles application data output.
  //TODO: Remove this later. Obsolete wrapper.
  template<class Simulator, class Solver>
  DataExporter* handleDataOutput(Simulator& simulator, Solver& solver,
				 bool append, const std::string& hdf5file,
				 int interval, int steps)
  { return handleDataOutput(simulator,solver,hdf5file,append,interval,steps); }
}

#endif
