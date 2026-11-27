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
#include "tinyxml2.h"
#include <algorithm>
#include <sstream>
#ifdef HAS_CEREAL
#include <cereal/cereal.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#endif


TimeStep::TimeStep () : step(0), iter(time.it), lstep(0), endStep(0)
{
  starTime = time.t = 0.0;
  stopTime = time.dt = 1.0;
  dtMin = dtMax = 1.0;
  nInitStep = 0;
  maxCFL = 0.0;
  f1 = 1.5;
  f2 = 0.25;
  maxStep = niter = 0;
  maxCut = nCut = 0;
  targetIts = stepIts = 0;
  dtSched = dtFull = 1.0;
  dtCut = dtHint = 0.0;
  limited = false;
  stepIt = mySteps.end();
}


TimeStep& TimeStep::operator= (const TimeStep& ts)
{
  step = ts.step;
  lstep = ts.lstep;
  endStep = ts.endStep;
  starTime = ts.starTime;
  stopTime = ts.stopTime;
  dtMin = ts.dtMin;
  dtMax = ts.dtMax;
  maxCFL = ts.maxCFL;
  nInitStep = ts.nInitStep;
  f1 = ts.f1;
  f2 = ts.f2;
  maxStep = ts.maxStep;
  maxCut = ts.maxCut;
  nCut = ts.nCut;
  targetIts = ts.targetIts;
  stepIts = ts.stepIts;
  niter = ts.niter;
  iter = ts.iter;
  dtSched = ts.dtSched;
  dtCut = ts.dtCut;
  dtHint = ts.dtHint;
  dtFull = ts.dtFull;
  limited = ts.limited;
  taken = ts.taken;

  time = ts.time;
  mySteps = ts.mySteps;
  stepIt = mySteps.begin() + (ts.stepIt - ts.mySteps.begin());
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


bool TimeStep::parse (const tinyxml2::XMLElement* elem)
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
  utl::getAttribute(elem,"maxCutbacks",maxCut);
  utl::getAttribute(elem,"target_iterations",targetIts);

  if (f1 < 1.0) f1 = 1.0;
  if (f2 > 1.0) f2 = 1.0;

  const tinyxml2::XMLElement* child = elem->FirstChildElement("step");
  for (; child; child = child->NextSiblingElement("step"))
  {
    double start = 0.0, end = 0.0;
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
      double dt = 0.0;
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
          IFEM::cout << ((i%10) ? " " : "\n\t") << timeStep[i];
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
  if (t > 0.0)
    return time.t + epsT*std::max(epsT,time.dt) > t;
  else if (t < 0.0)
    return time.t < 0.0 && time.t - epsT*epsT < t;
  else
    return true;
}


bool TimeStep::reset (int istep)
{
  nCut = 0;
  if (mySteps.empty()) {
    // Without time step definitions, the time and step size limits are kept
    step = istep;
    return true;
  }

  // The replay of the steps starts without limits, and a suggested size of
  // the step after the last one taken is kept for it
  if (static_cast<size_t>(step) == taken.dt.size())
    taken.dtHint = dtHint;
  dtCut = dtHint = 0.0;
  limited = false;
  lstep = endStep = step = 0;
  stepIt = mySteps.begin();
  if (istep == 0)
    stopTime = mySteps.back().second;
  time.dt = mySteps.front().first.front();
  time.dtn = time.t = time.CFL = 0.0;
  for (int i = 0; i < istep; i++)
    if (!this->increment())
      return false;

  return true;
}


void TimeStep::scheduleIncrement ()
{
  if (limited && mySteps.empty() && maxCFL <= 0.0)
    time.dt = dtSched; // Continue from the step size before it was limited

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
  {
    if (++lstep <= stepIt->first.size())
      time.dt = stepIt->first[lstep-1];
    else if (step > 0)
      time.dt = dtSched; // Beyond the definition, as scheduled before
  }

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

  if (dtHint > 0.0 && mySteps.empty() && maxCFL <= 0.0)
  {
    // Without time step definitions, a suggested step size replaces it,
    // within the step size bounds
    time.dt = dtHint;
    if (dtMin < dtMax)
      time.dt = std::max(std::min(dtHint,dtMax),dtMin);
    dtHint = 0.0;
  }

  // Limit the scheduled step size after a cut-back, or as suggested,
  // a suggestion not below the minimum step size
  dtSched = time.dt;
  double dtLim = dtHint;
  if (dtLim > 0.0 && dtMin < dtMax && dtLim < dtMin)
    dtLim = dtMin;
  if (dtCut > 0.0)
  {
    // Let the step size grow back to the scheduled one after a cut-back
    dtCut *= f1;
    if (dtCut < dtSched)
      dtLim = dtLim > 0.0 ? std::min(dtLim,dtCut) : dtCut;
    else
      dtCut = 0.0;
  }
  dtHint = 0.0;
  if ((limited = dtLim > 0.0 && dtLim < dtSched))
    time.dt = dtLim;
}


bool TimeStep::increment ()
{
  if (targetIts > 0 && stepIts > 0)
    this->suggestNextStep(dtFull*targetIts/stepIts);
  stepIts = 0;

  time.dtn = time.dt;

  const bool replay = static_cast<size_t>(step) < taken.dt.size();
  const size_t lprev = lstep;
  if (replay)
  {
    time.dt = taken.dt[step];
    dtSched = taken.sched[step];
    dtHint = 0.0;
    if (stepIt != mySteps.end() && taken.defined[step])
      ++lstep;
  }
  else
    this->scheduleIncrement();
  const char defined = lstep != lprev;

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
  nCut = 0;
  endStep = 0;
  dtFull = time.dt;
  time.t += time.dt;

  if (stepIt != mySteps.end() && this->hasReached(stepIt->second))
  {
    if (time.t != stepIt->second)
    {
      // Adjust the size of the last time step
      time.dt += stepIt->second - time.t;
      time.t = stepIt->second;
    }
    endStep = lstep+1;
    lstep = 0;
    ++stepIt;
  }

  if (time.t > stopTime)
  {
    // Adjust the size of the last time step
    time.dt += stopTime - time.t;
    time.t = stopTime;
  }

  if (!replay)
  {
    // Record the step, which a replay of the steps repeats
    taken.dt.push_back(time.dt);
    taken.sched.push_back(dtSched);
    taken.defined.push_back(defined);
    taken.dtCut = dtCut;
    taken.limited = limited;
    taken.niter = niter;
    taken.dtHint = 0.0;
  }
  else if (static_cast<size_t>(step) == taken.dt.size())
  {
    // The replay continues as the last step taken did
    dtCut = taken.dtCut;
    limited = taken.limited;
    niter = taken.niter;
    dtHint = taken.dtHint;
  }

  return true;
}


bool TimeStep::cutback ()
{
  // Without a maximum number of cut-backs, the minimum step size is the limit
  const bool useMin = maxCut < 1 || dtMin < dtMax;
  if ((useMin && time.dt <= dtMin) || (maxCut > 0 && nCut >= maxCut))
  {
    // Already reached the minimum step size, cannot do further cut-back
    IFEM::cout <<" *** Iterations diverged, terminating..."<< std::endl;
    return false;
  }

  // Reduce the time step size to the suggested one, if smaller,
  // or else by the factor f2 < 1.0 or 0.5
  double dt = time.dt;
  if (dtHint > 0.0 && dtHint < dt)
    time.dt = dtHint;
  else
    time.dt *= (f2 < 1.0 ? f2 : 0.5);
  dtHint = 0.0;
  if (useMin && time.dt < dtMin)
    time.dt = dtMin;

  ++nCut;
  iter = 0;
  stepIts = 0;
  niter = 10;
  time.first = 'c';
  time.t += time.dt - dt;
  dtFull = time.dt;

  // The reduced step no longer ends the time step definition it ended
  if (endStep > 0)
  {
    --stepIt;
    lstep = endStep-1;
    endStep = 0;
  }

  if (maxCut > 0)
  {
    // Limit the following steps, which grow back to the scheduled step size
    dtCut = time.dt;
    limited = true;
  }
  else
    dtSched = time.dt; // Continue with the reduced step size

  // The step taken is the reduced one, and a replay of the steps ends with it
  if (step > 0)
  {
    taken.dt.truncate(step);
    taken.sched.truncate(step);
    taken.defined.truncate(step);
    taken.dt.setBack(time.dt);
    taken.sched.setBack(dtSched);
    taken.dtCut = dtCut;
    taken.limited = limited;
    taken.niter = niter;
    taken.dtHint = 0.0;
  }

  IFEM::cout <<"  ** Iterations diverged, trying cut-back with dt="
             << time.dt << std::endl;
  return true;
}


void TimeStep::suggestNextStep (double dt)
{
  if (dt > 0.0 && (dtHint <= 0.0 || dt < dtHint))
    dtHint = dt;
}


bool TimeStep::finished () const
{
  if (step >= maxStep && maxStep > 0)
    return true;

  // Also check negative stopTime in case of path-following algorithms
  return this->hasReached(stopTime) || this->hasReached(-stopTime);
}

#ifdef HAS_CEREAL
//! \brief Serializes TimeStep data \a tp to/from the archive \a ar.
template<class T> void doSerializeOps (T& ar, TimeStep& tp)
{
  ar(tp.step);
  ar(tp.maxCFL);
  ar(tp.time.t);
  ar(tp.time.dt);
  ar(tp.time.dtn);
  ar(tp.time.CFL);
  ar(tp.time.first);
}


//! \brief Saves the runs of values \a v to the archive \a ar.
template<class Archive, class T>
void saveRuns (Archive& ar, const T& v)
{
  ar(v.vals,v.ends);
}


//! \brief Loads the runs of values \a v from the archive \a ar.
template<class Archive, class T>
bool loadRuns (Archive& ar, T& v)
{
  ar(v.vals,v.ends);
  if (v.ends.size() != v.vals.size())
    return false;

  for (size_t i = 0; i < v.ends.size(); i++)
    if (v.ends[i] <= (i > 0 ? v.ends[i-1] : 0))
      return false;

  return true;
}
#endif


bool TimeStep::serialize (std::map<std::string,std::string>& data) const
{
#ifdef HAS_CEREAL
  std::ostringstream str;
  {
    cereal::BinaryOutputArchive ar(str);
    doSerializeOps(ar,*const_cast<TimeStep*>(this));
  }
  data.insert(std::make_pair("TimeStep",str.str()));

  // The position in the time step definitions, the step size limits, and the
  // time steps taken, for a replay of the steps. Their sizes are run-length
  // encoded, as they mostly repeat a scheduled size.
  std::ostringstream str2;
  {
    cereal::BinaryOutputArchive ar(str2);
    const int seg = stepIt - mySteps.begin();
    ar(seg,lstep,iter,niter,dtSched,dtCut,dtHint,limited);
    saveRuns(ar,taken.dt);
    saveRuns(ar,taken.sched);
    saveRuns(ar,taken.defined);
    ar(taken.dtCut,taken.limited,taken.niter,taken.dtHint);
  }
  data.emplace("TimeStep::schedule",str2.str());

  std::ostringstream str3;
  {
    cereal::BinaryOutputArchive ar(str3);
    ar(stepIts,dtFull);
  }
  data.emplace("TimeStep::iterations",str3.str());
  return true;
#else
  return false;
#endif
}


bool TimeStep::deSerialize (const std::map<std::string,std::string>& data)
{
#ifdef HAS_CEREAL
  std::map<std::string,std::string>::const_iterator it = data.find("TimeStep");
  if (it != data.end()) {
    std::stringstream str(it->second);
    cereal::BinaryInputArchive ar(str);
    doSerializeOps(ar,*this);

    it = data.find("TimeStep::schedule");
    if (it == data.end())
    {
      taken = Taken();
      return this->reset(step);
    }

    std::stringstream str2(it->second);
    cereal::BinaryInputArchive ar2(str2);
    int seg = 0;
    ar2(seg,lstep,iter,niter,dtSched,dtCut,dtHint,limited);
    if (seg < 0 || seg > static_cast<int>(mySteps.size()))
      return false;

    if (!loadRuns(ar2,taken.dt) || !loadRuns(ar2,taken.sched) ||
        !loadRuns(ar2,taken.defined))
      return false;
    ar2(taken.dtCut,taken.limited,taken.niter,taken.dtHint);
    if (taken.sched.size() != taken.dt.size() ||
        taken.defined.size() != taken.dt.size())
      return false;

    stepIts = 0;
    dtFull = time.dt;
    if ((it = data.find("TimeStep::iterations")) != data.end())
    {
      std::stringstream str3(it->second);
      cereal::BinaryInputArchive ar3(str3);
      ar3(stepIts,dtFull);
    }

    stepIt = mySteps.begin() + seg;
    endStep = 0;
    nCut = 0;
    return true;
  }
#endif
  return false;
}
