// $Id$
//==============================================================================
//!
//! \file TimeStep.h
//!
//! \date Feb 11 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for encapsulation of general time stepping parameters.
//!
//==============================================================================

#ifndef TIME_STEP_H_
#define TIME_STEP_H_

#include "TimeDomain.h"
#include <vector>
#include <cstddef>
#include <iostream>
#include <map>
#include <string>

class TiXmlElement;


/*!
  \brief Class for encapsulation of general time stepping parameters.
*/

class TimeStep
{
public:
  //! \brief The constructor initializes the counters to zero.
  TimeStep();
  //! \brief Copy constructor.
  TimeStep(const TimeStep& ts) : iter(time.it) { *this = ts; }
  //! \brief Assigment operator.
  TimeStep& operator=(const TimeStep& ts);

  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  bool parse(char* keyWord, std::istream& is);
  //! \brief Parses a data section from an XML document.
  //! \param[in] elem The XML element to parse
  bool parse(const TiXmlElement* elem);

  //! \brief Returns \e true if the simulation consists of several time steps.
  bool multiSteps() const;
  //! \brief Returns \e true if the given time \a t has been reached.
  bool hasReached(double t) const;

  //! \brief Resets the time step to the specified step.
  //! \return \e false, if the \a istep is passed the end of the simulation
  bool reset(int istep = 0);
  //! \brief Advances the time increments one step further.
  //! \return \e true, if we have reached the end of the simulation
  bool increment();
  //! \brief Restarts current increment with a smaller step size on divergence.
  //! \return \e false Cannot do further cut-back, time step size too small
  bool cutback();

  //! \brief Serialize internal state for restarting purposes.
  //! \param data Container for serialized data
  bool serialize(std::map<std::string,std::string>& data);
  //! \brief Set internal state from a serialized state.
  //! \param[in] data Container for serialized data
  bool deSerialize(const std::map<std::string,std::string>& data);

  int        step; //!< Time step counter
  int&       iter; //!< Iteration counter
  TimeDomain time; //!< Time domain data

  double starTime; //!< Start time of simulation
  double stopTime; //!< Stop time of simulation
  double maxCFL;   //!< CFL restriction on time step size (0.0: no restriction)

private:
  int    niter;      //!< Number of iterations in previous time step
  int    nInitStep;  //!< Number of fixed timesteps in the beginning
  int    maxStep;    //!< Maximum number of time steps
  double dtMin;      //!< Minimum time increment size
  double dtMax;      //!< Maximun time increment size
  double f1;         //!< Scale factor for increased time step size
  double f2;         //!< Scale factor for reduced time step size

  typedef std::pair<std::vector<double>,double> Step; //!< Time step definition
  typedef std::vector<Step> TimeSteps;                //!< Time step container

  TimeSteps           mySteps; //!< Time step definitions
  TimeSteps::iterator stepIt;  //!< Running iterator over the time steps

  size_t lstep; //!< Local step counter, i.e., within current \a *stepIt
};

#endif
