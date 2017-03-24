// $Id$
//==============================================================================
//!
//! \file TimeStep.C
//!
//! \date Feb 11 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for encapsulation of general time stepping parameters.
//!
//==============================================================================

#include "TimeStep.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml.h"
#include <sstream>
#ifdef HAS_CEREAL
#include <cereal/cereal.hpp>
#include <cereal/archives/binary.hpp>
#endif


TimeStep::TimeStep () : step(0), iter(time.it), lstep(0)
{
  starTime = time.t = 0.0;
  stopTime = time.dt = 1.0;
  dtMin = dtMax = 1.0;
  nInitStep = 0;
  maxCFL = 0.0;
  f1 = 1.5;
  f2 = 0.25;
  maxStep = niter = 0;
  stepIt = mySteps.end();
}


TimeStep& TimeStep::operator= (const TimeStep& ts)
{
  step = ts.step;
  lstep = ts.lstep;
  starTime = ts.starTime;
  stopTime = ts.stopTime;
  dtMin = ts.dtMin;
  dtMax = ts.dtMax;
  maxCFL = ts.maxCFL;
  nInitStep = ts.nInitStep;
  f1 = ts.f1;
  f2 = ts.f2;
  maxStep = ts.maxStep;
  niter = ts.niter;

  time = ts.time;
  mySteps = ts.mySteps;
  stepIt = mySteps.begin();
  return *this;
}


bool TimeStep::parse (char* keyWord, std::istream& is)
{
  if (strncasecmp(keyWord,"TIME_STEPPING",13))
    return true;

  int nstep = atoi(keyWord+13);
  if (nstep < 1) nstep = 1;

  double dt;
  mySteps.resize(nstep);
  for (int i = 0; i < nstep; i++)
  {
    std::istringstream cline(utl::readLine(is));
    if (i == 0) cline >> starTime;
    cline >> mySteps[i].second >> dt;
    if (i == 0 && cline.good()) cline >> maxCFL;
    if (cline.fail() || cline.bad())
      return false;

    if (dt > 1.0 && ceil(dt) == dt)
    {
      // The number of time steps are specified
      dt = (mySteps[i].second - (i == 0 ? starTime : mySteps[i-1].second))/dt;
      mySteps[i].first.push_back(dt);
    }
    else while (!cline.fail() && !cline.bad())
    {
      // The time step size(s) is/are specified
      mySteps[i].first.push_back(dt);
      cline >> dt;
    }
  }

  this->reset();

  time.t = starTime;
  return true;
}


bool TimeStep::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"timestepping"))
    return true;

  utl::getAttribute(elem,"start",starTime);
  utl::getAttribute(elem,"end",stopTime);
  utl::getAttribute(elem,"dtMin",dtMin);
  utl::getAttribute(elem,"dtMax",dtMax);
  utl::getAttribute(elem,"maxCFL",maxCFL);
  utl::getAttribute(elem,"nInitStep",nInitStep);
  utl::getAttribute(elem,"maxStep",maxStep);
  utl::getAttribute(elem,"f1",f1);
  utl::getAttribute(elem,"f2",f2);

  if (f1 < 1.0) f1 = 1.0;
  if (f2 > 1.0) f2 = 1.0;

  const TiXmlElement* child = elem->FirstChildElement("step");
  for (; child; child = child->NextSiblingElement())
  {
    double start = 0.0, end = 0.0, dt = 0.0;
    utl::getAttribute(child,"start",start);
    utl::getAttribute(child,"end",end);
    if (mySteps.empty())
      starTime = start;
    else
      start = mySteps.back().second;

    if (child->FirstChild())
    {
      std::vector<double> timeStep;
      std::istringstream cline(child->FirstChild()->Value());
      cline >> dt;
      if (dt > 1.0 && ceil(dt) == dt)
        // The number of steps are specified
        timeStep.push_back((end-start)/dt);

      else while (!cline.fail() && !cline.bad())
      {
        // The time step size(s) is/are specified
        timeStep.push_back(dt);
        cline >> dt;
      }

      bool graded = false;
      utl::getAttribute(child,"graded",graded);
      if (graded && timeStep.size() == 2)
      {
        // Geometric grading
        double dt  = timeStep.front();
        double dt2 = timeStep.back();
        double eta = 1.0 - (dt - dt2)/(end - start);
        timeStep.resize(1);
        for (double T = start+dt; T < end; T += dt)
        {
          dt *= eta;
          timeStep.push_back(T+dt <= end ? dt : end-T);
        }
        IFEM::cout <<"\tGeometric graded time increments in ["
                   << start <<","<< end <<"]:\n\t"<< timeStep.front();
        for (size_t i = 1; i < timeStep.size(); i++)
          IFEM::cout << (i%10 ? " " : "\n\t") << timeStep[i];
        IFEM::cout << std::endl;
      }

      if (!timeStep.empty())
        mySteps.push_back(std::make_pair(timeStep,end));
    }
  }

  if (!mySteps.empty())
    this->reset();
  else if (!utl::getAttribute(elem,"dt",time.dt))
    time.dt = stopTime - starTime;

  time.t = starTime;
  return true;
}


bool TimeStep::multiSteps () const
{
  const double epsT = 1.0e-6;
  return starTime + (1.0+epsT)*time.dt < stopTime;
}


bool TimeStep::hasReached (double t) const
{
  const double epsT = 1.0e-6;
  return time.t + epsT*std::max(epsT,time.dt) > t;
}


bool TimeStep::reset (int istep)
{
  lstep = step = 0;
  stepIt = mySteps.begin();
  stopTime = mySteps.back().second;
  time.dt = mySteps.front().first.front();
  time.dtn = time.t = time.CFL = 0.0;
  for (int i = 0; i < istep; i++)
    if (!this->increment())
      return false;

  return true;
}


bool TimeStep::increment ()
{
  time.dtn = time.dt;

  if (maxCFL > 0.0 && time.CFL > 1.0e-12) {
    // Increase CFL by a given factor
    if (step > nInitStep) {
      double dt = maxCFL/time.CFL;
      if (dt > time.dt*f1)
        time.dt *= f1;
      else {
        time.dt = dt;
        maxCFL *= f1;
      }
    }
  }
  else if (stepIt != mySteps.end())
    if (++lstep <= stepIt->first.size())
      time.dt = stepIt->first[lstep-1];

  if (dtMin < dtMax && maxCFL <= 0.0 && step > 1)
  {
    // Adjust the time step size based on the number of iterations in last step
    if (iter <= 4 && niter <= 4)
      time.dt *= f1;
    else if (iter > 10)
      time.dt *= f2;

    if (time.dt < dtMin)
      time.dt = dtMin;
    else if (time.dt > dtMax)
      time.dt = dtMax;
  }

  if (this->hasReached(stopTime))
  {
    if (step > 0)
      IFEM::cout <<"\n  Time integration completed."<< std::endl;
    return false;
  }
  else if (step++ == maxStep && maxStep > 0)
  {
    IFEM::cout <<"\n  ** Terminating, maximum number of time steps reached."
               << std::endl;
    return false;
  }

  niter = iter;
  time.t += time.dt;

  if (stepIt != mySteps.end())
  {
    if (stepIt->first.size() <= lstep)
      stepIt->first.push_back(time.dt);

    if (this->hasReached(stepIt->second))
    {
      if (time.t != stepIt->second)
      {
        // Adjust the size of the last time step
        time.dt += stepIt->second - time.t;
        time.t = stepIt->second;
        stepIt->first.back() = time.dt;
      }
      lstep = 0;
      ++stepIt;
    }
  }

  if (time.t > stopTime)
  {
    // Adjust the size of the last time step
    time.dt += stopTime - time.t;
    time.t = stopTime;
  }

  return true;
}


bool TimeStep::cutback ()
{
  if (time.dt <= dtMin)
  {
    // Already reached the minimum step size, cannot do further cut-back
    IFEM::cout <<" *** Iterations diverged, terminating..."<< std::endl;
    return false;
  }

  // Reduce the time step size by the factor f2 < 1.0 or 0.5
  double dt = time.dt;
  time.dt *= (f2 < 1.0 ? f2 : 0.5);
  if (time.dt < dtMin)
    time.dt = dtMin;

  iter = 0;
  niter = 10;
  time.first = 'c';
  time.t += time.dt - dt;
  if (stepIt != mySteps.end())
    stepIt->first.back() = time.dt;

  IFEM::cout <<"  ** Iterations diverged, trying cut-back with dt="
             << time.dt << std::endl;
  return true;
}


#ifdef HAS_CEREAL
/*!
  \brief Serialize to/from archive.
  \param ar Input or output archive
*/

template<class T> void doSerializeOps(T& ar, TimeStep& tp)
{
  ar(tp.step);
  ar(tp.starTime);
  ar(tp.maxCFL);
  ar(tp.time.t);
  ar(tp.time.dt);
  ar(tp.time.dtn);
  ar(tp.time.CFL);
  ar(tp.time.first);
}
#endif


bool TimeStep::serialize(std::map<std::string,std::string>& data)
{
#ifdef HAS_CEREAL
  std::ostringstream str;
  cereal::BinaryOutputArchive ar(str);
  doSerializeOps(ar,*this);
  data.insert(std::make_pair("TimeStep", str.str()));
  return true;
#endif
  return false;
}


bool TimeStep::deSerialize(const std::map<std::string,std::string>& data)
{
#ifdef HAS_CEREAL
  std::stringstream str;
  auto it = data.find("TimeStep");
  if (it != data.end()) {
    str << it->second;
    cereal::BinaryInputArchive ar(str);
    doSerializeOps(ar,*this);
    return true;
  }
#endif
  return false;
}
