// $Id$
//==============================================================================
//!
//! \file Profiler.C
//!
//! \date May 18 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Simple profiling of computational tasks.
//!
//==============================================================================

#include "Profiler.h"
#ifdef PARALLEL_PETSC
#include "petscsys.h"
#endif
#include <sys/time.h>

#include "LinAlgInit.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif


Profiler* utl::profiler = 0;


Profiler::Profiler (const std::string& name) : myName(name), nRunners(0)
{
  this->start("Total");

  allCPU = allWall = 0.0;

  // Update pointer to current profiler (it should only be one at any time)
  if (utl::profiler) delete utl::profiler;
  utl::profiler = this;
  LinAlgInit::increfs();
}


Profiler::~Profiler ()
{
  this->stop("Total");
  this->report(std::cout);
  LinAlgInit::decrefs();
}


//! \brief Returns the current wall time in seconds and resolution in microsec.

static double WallTime ()
{
#ifdef USE_OPENMP
  return omp_get_wtime();
#endif
  timeval tmpTime;
  gettimeofday(&tmpTime,NULL);
  return tmpTime.tv_sec + tmpTime.tv_usec/1.0e6;
}


void Profiler::start (const std::string& funcName)
{
  Profile& p = myTimers[funcName];
  if (p.running) return;

  p.running = true;
  nRunners++;
  p.nCalls++;
  p.startWall = WallTime();
  p.startCPU = clock();
}


void Profiler::stop (const std::string& funcName)
{
  clock_t stopCPU = clock();
  double stopWall = WallTime();

  std::map<std::string,Profile>::iterator it = myTimers.find(funcName);
  if (it == myTimers.end())
    std::cerr <<" *** No matching timer for "<< funcName << std::endl;
  else if (it->second.running)
  {
    // Accumulate consumed CPU and wall time by this task
    Profile& p = it->second;
    double deltaCPU  = double(stopCPU - p.startCPU)/double(CLOCKS_PER_SEC);
    double deltaWall = stopWall - p.startWall;
    p.running = false;
    p.totalCPU  += deltaCPU;
    p.totalWall += deltaWall;
    if (--nRunners == 1)
    {
      // This is a "main" task, accumulate the total time for all main tasks
      allCPU  += deltaCPU;
      allWall += deltaWall;
    }
  }
}


static bool use_ms = false; //!< Print mean times in microseconds?


//! \brief Global stream operator printing a Profile instance.

std::ostream& operator<<(std::ostream& os, const Profiler::Profile& p)
{
  // First the CPU time
  os.width(10);
  os << p.totalCPU;
  if (p.nCalls > 1)
  {
    os.width(9);
    os << (use_ms ? 1000.0 : 1.0)*p.totalCPU/p.nCalls <<" |";
  }
  else
    os <<"          |";

  // Then the wall time and the number of invokations (if more than one)
  os.width(10);
  os << p.totalWall;
  if (p.nCalls > 1)
  {
    os.width(9);
    os << (use_ms ? 1000.0 : 1.0)*p.totalWall/p.nCalls <<" |";
    os.width(6);
    os << p.nCalls;
  }
  else
    os <<"          |";

  return os <<'\n';
}


void Profiler::report (std::ostream& os) const
{
  if (myTimers.empty()) return;

  use_ms = true; // Print mean times in microseconds by default
  std::map<std::string,Profile>::const_iterator it;
  for (it = myTimers.begin(); it != myTimers.end(); it++)
  {
    // Make sure the task has stopped profiling (in case of exceptions)
    if (it->second.running) const_cast<Profiler*>(this)->stop(it->first);
    if (it->second.nCalls > 1)
      if (it->second.totalWall/it->second.nCalls >= 100.0)
	use_ms = false; // Print mean times in seconds
  }

  // Find the time for "other" tasks, i.e., the difference between
  // the measured total time and the sum of all the measured tasks
  Profile other;
  std::map<std::string,Profile>::const_iterator tit = myTimers.find("Total");
  if (tit != myTimers.end())
  {
    if (!tit->second.haveTime()) return; // Nothing to report, run in zero time
    other.totalCPU  = tit->second.totalCPU  - allCPU;
    other.totalWall = tit->second.totalWall - allWall;
  }

  int myId = -1;
#ifdef PARALLEL_PETSC
  MPI_Comm_rank(PETSC_COMM_WORLD,&myId);
#endif

  // Print a table with timing results, all tasks with zero time are ommitted
  const char* Ms = (use_ms ? "Mean(ms)" : "Mean(s) ");
  os <<"\n==============================================================="
     <<"\n===   Profiling results for "<< myName;
  if (myId >= 0) os <<" on processor "<< myId;
  os <<"\n================================================================="
     <<"\n                      |       CPU time     |      Wall time     |"
     <<"\nTask                  |  Total(s)  "<<Ms<<"|  Total(s)  "<<Ms<<"| calls"
     <<"\n----------------------+--------------------+--------------------+------"
     <<"\n";
  os.precision(2);
  os.flags(std::ios::fixed|std::ios::right);
  for (it = myTimers.begin(); it != myTimers.end(); it++)
    if (it != tit && it->second.haveTime())
      if (it->first.size() >= 22)
	os << it->first.substr(0,22) <<'|'<< it->second;
      else
	os << it->first << std::string(22-it->first.size(),' ')
	   <<'|'<< it->second;

  // Finally, print the "other" and "total" times
  if (other.haveTime())
    os <<"Other                 |"<< other;
  if (tit != myTimers.end())
    os <<"----------------------+--------------------+--------------------+------"
       <<"\nTotal time            |"<< tit->second;
  os <<"======================================================================="
     << std::endl;
}
