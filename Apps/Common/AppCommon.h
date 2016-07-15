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
#include "IFEM.h"
#include "XMLInputBase.h"

namespace SIM
{
  /*! \brief Base class for input file pre-parsing in applications.
   */

  class AppXMLInputBase : public XMLInputBase
  {
  public:
    int dim = 3; //!< Dimensionality of simulation

  protected:
    //! \brief Parse a tag for input file
    virtual bool parse(const TiXmlElement* elem) override
    {
      if (!strcasecmp(elem->Value(),"geometry")) {
        if (!utl::getAttribute(elem, "dimension", dim))
          utl::getAttribute(elem, "dim", dim);
        const TiXmlElement* child = elem->FirstChildElement();
        for (; child; child= child->NextSiblingElement())
          if (!strcasecmp(child->Value(),"patchfile")) {
            std::string type;
            utl::getAttribute(child,"type",type);
            if (type == "lrspline")
              IFEM::getOptions().discretization = ASM::LRSpline;
          }
      }

      return true;
    }
  };


  //! \brief Handles application restarts.
  //! \param[in] simulator The top SIMbase instance of your application
  //! \param[in] solver The SIMSolver instance of your application
  //! \param[in] restartfile The file to read from
  //! \param[in] interval The stride in the input file
  //! \param[in] steps The number of time steps to load
  template<class Simulator, class Solver>
  void handleRestart(Simulator& simulator, Solver& solver,
                     const std::string& restartfile,
                     int interval = 1, int steps = 1)
  {
    DataExporter reader(true,interval,steps);
    XMLWriter* xml = new XMLWriter(restartfile,solver.getProcessAdm());
    HDF5Writer* hdf = new HDF5Writer(restartfile,solver.getProcessAdm(),true);
    reader.registerWriter(xml);
    reader.registerWriter(hdf);
    simulator.registerFields(reader);
    int max = reader.getTimeLevel();
    // correct loaded level if we stopped in the middle of a "stride" level
    if ((max+1) % (steps+1))
      max -= (max % (steps+1))+steps;
    double time;
    hdf->openFile(max-steps);
    hdf->readDouble(max-steps,"timeinfo","SIMbase-1",time);
    solver.fastForward(time/solver.getTimePrm().time.dt);
    for (int i=steps;i>=0;--i) {
      reader.loadTimeLevel(max-i,xml,hdf);
      solver.postSolve(solver.getTimePrm(),true);
      if (i > 0) solver.advanceStep();
    }
    xml->writeTimeInfo(0, interval, steps, solver.getTimePrm());
  }

  //! \brief Handles application data output.
  //! \param[in] simulator The top SIMbase instance of your application
  //! \param[in] solver The SIMSolver instance of your application
  //! \param[in] hdf5file The file to save to
  //! \param[in] append Whether or not to append to file
  //! \param[in] interval The stride in the output file
  //! \param[in] steps The number of time steps to dump in onw row
  template<class Simulator, class Solver>
  DataExporter* handleDataOutput(Simulator& simulator, Solver& solver,
                                 const std::string& hdf5file,
                                 bool append = false,
                                 int interval = 1, int steps = 1)
  {
    DataExporter* writer = new DataExporter(true,interval,steps);
    XMLWriter* xml = new XMLWriter(hdf5file,solver.getProcessAdm());
    HDF5Writer* hdf = new HDF5Writer(hdf5file,solver.getProcessAdm(),append);
    writer->registerWriter(xml);
    writer->registerWriter(hdf);
    simulator.registerFields(*writer);
    IFEM::registerCallback(*writer);
    return writer;
  }
}

#endif
