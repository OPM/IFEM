// $Id$
//==============================================================================
//!
//! \file SIMoptions.h
//!
//! \date Feb 13 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Class for encapsulation of general simulation options.
//!
//==============================================================================

#ifndef _SIM_OPTIONS_H_
#define _SIM_OPTIONS_H_

#include "ASMenums.h"
#include "LinAlgenums.h"
#include <string>
#include <map>

namespace utl {
  class LogStream;
}
class TiXmlElement;


/*!
  \brief Class for encapsulation of general simulation options.

  \details This class is supposed to contain all global simulation options
  that not already are stored in the model objects themselves. By keeping them
  independent of the model they can be initialized before instanciating
  the application-dependent model object(s). The class is equipped with some
  methods to initialize the options through the parsing of XML-tags.
*/

class SIMoptions
{
public:
  //! \brief The constructor initializes the default input options.
  SIMoptions();

  //! \brief Defines the linear equation solver to be used.
  void setLinearSolver(const std::string& eqsolver);

  //! \brief Parses a subelement of the \a console XML-tag.
  bool parseConsoleTag(const TiXmlElement* elem);
  //! \brief Parses a subelement of the \a discretization XML-tag.
  bool parseDiscretizationTag(const TiXmlElement* elem);
  //! \brief Parses a subelement of the \a eigensolver XML-tag.
  bool parseEigSolTag(const TiXmlElement* elem);
  //! \brief Parses a subelement of the \a resultoutput XML-tag.
  bool parseOutputTag(const TiXmlElement* elem);
  //! \brief Parses the \a restart XML-tag.
  bool parseRestartTag(const TiXmlElement* elem);

  //! \brief Parses obsolete command-line arguments (backward compatibility).
  bool parseOldOptions(int argc, char** argv, int& i);
  //! \brief Returns \e true if the i'th argument is obsolete.
  static bool ignoreOldOptions(int argc, char** argv, int& i);

  //! \brief Returns whether HDF5 output is requested or not.
  bool dumpHDF5(const char* defaultName);

  //! \brief Prints out the simulation options to the given stream.
  utl::LogStream& print(utl::LogStream& os, bool addBlankLine = false) const;

  //! \brief Parses a projection method XML-tag.
  bool parseProjectionMethod(const char* ptype, int version = 1);

public:
  int nGauss[2]; //!< Gaussian quadrature rules

  ASM::Discretization discretization; //!< Spatial discretization option
  LinAlg::MatrixType  solver;         //!< The linear equation solver to use

  int num_threads_SLU; //!< Number of threads for SuperLU_MT

  // Eigenvalue solver options
  int    eig;   //!< Eigensolver method (1,...,5)
  int    nev;   //!< Number of eigenvalues/vectors
  int    ncv;   //!< Number of Arnoldi vectors
  double shift; //!< Eigenvalue shift

  // Output options
  int format;    //!< VTF-file format (-1=NONE, 0=ASCII, 1=BINARY)
  int nViz[3];   //!< Number of visualization points over each knot-span
  int saveInc;   //!< Number of load/time increments between each result output
  double dtSave; //!< Time interval between each result output
  bool pSolOnly; //!< If \e true, don't save secondary solution variables
  bool saveNorms;//!< If \e true, save element norms

  std::string hdf5; //!< Prefix for HDF5-file
  std::string vtf; //!< Prefix for VTF file

  // Restart options
  int         restartInc;  //!< Number of increments between each restart output
  int         restartStep; //!< Index to the actual state to restart from
  std::string restartFile; //!< File to read restart state data from

  bool enableController; //!< Whether or not to enable external program control

  int         printPid;   //!< PID to print info to screen for
  std::string log_prefix; //!< Prefix for process log files

  //! \brief Enum defining the available projection methods.
  enum ProjectionMethod { NONE, GLOBAL, VDSA, QUASI, LEASTSQ,
                          DGL2, CGL2, CGL2_INT, SCR };
  //! \brief Projection method name mapping.
  typedef std::map<ProjectionMethod,std::string> ProjectionMap;

  ProjectionMap project; //!< The projection methods to use
};

#endif
