// $Id$
//==============================================================================
//!
//! \file Profiler.h
//!
//! \date May 18 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Simple profiling of computational tasks.
//!
//==============================================================================

#ifndef _PROFILER_H
#define _PROFILER_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <ctime>


/*!
  \brief Simple class for profiling of computational tasks.

  \details The profiler measures the CPU time and the wall clock time for
  computational tasks identified be the user. Each task may be invoked an
  arbitrary number of times, and the average time consumption is then also
  recorded along with the number of invokations.

  The profiling results are printed in a nicely formatted table when the
  profiler object goes out of scope, typically at the end of the program.
*/

class Profiler
{
public:
  //! \brief The constructor initializes the profiler object.
  //! \param[in] name Program name to be printed in the profiling report header.
  //!
  //! \details The constructor also updates the global static pointer
  //! utl::profiler to point to \a *this, deleting any already pointed-to
  //! object first. This means, only one Profiler object can exist at any time.
  Profiler(const std::string& name);
  //! \brief The destructor prints the profiling report to the console.
  ~Profiler();

  //! \brief Starts profiling of task \a funcName and increments \a nRunners.
  void start(const std::string& funcName);
  //! \brief Stops profiling of task \a funcName and decrements \a nRunners.
  void stop(const std::string& funcName);

  //! \brief Prints a profiling report for all tasks that have been measured.
  void report(std::ostream& os) const;
  //! \brief Clears the profiler.
  void clear() { myTimers.clear(); allCPU = allWall = 0.0; nRunners = 0; }

private:
  //! \brief Stores profiling data for one computational task.
  struct Profile
  {
    clock_t startCPU;  //!< The last starting CPU time of this task
    double  startWall; //!< The last starting wall clock time of this task
    double  totalCPU;  //!< Total CPU time consumed by this task so far
    double  totalWall; //!< Total wall clock time consumed by this task so far
    size_t  nCalls;    //!< Number of invokations of this task
    bool    running;   //!< Flag indicating if this task is currently running

    //! \brief The constructor initializes the total times to zero.
    Profile() : nCalls(0), running(false) { totalCPU = totalWall = 0.0; }
    //! \brief Checks if this profile item have any timing to report.
    bool haveTime() const { return totalCPU >= 0.005 || totalWall >= 0.005; }
  };

  //! \brief Global stream operator printing a Profile instance.
  friend std::ostream& operator<<(std::ostream& os, const Profile& p);

  std::string myName; //!< Name of this profiler

  typedef std::map<std::string,Profile> ProfileMap; //!< Map of profilers

  ProfileMap              myTimers;  //!< The task profiles with names
  std::vector<ProfileMap> myMTimers; //!< Task profiles with names and thread iD

  double allCPU;  //!< Accumulated CPU time from all "main" tasks
  double allWall; //!< Accumulated wall clock time of all "main" tasks

  //! \details When \a nRunners is 1, we are only measuring the total time.
  //! When \a nRunners is 2, we are also measuring a "main" task, which total
  //! time is subtracted from the measured total time to find the reminder of
  //! all tasks that are not being measured. When \a nRunners is 3 or larger,
  //! the measured task is already included in another "main" task,
  //! and therefore its time is not included when calculating the "other" times.
  size_t nRunners; //!< Number of tasks currently running
};


namespace utl
{
  extern Profiler* profiler; //!< Pointer to the one and only profiler object.

  //! \brief Convenience class to profile the local scope.
  class prof
  {
    const char* name; //!< Name tag on the local scope to profile
  public:
    //! \brief The constructor starts the profiling of the named task.
    prof(const char* tag) : name(tag) { if (profiler) profiler->start(name); }
    //! \brief The destructor stops the profiling.
    ~prof() { if (profiler) profiler->stop(name); }
  };
}


//! \brief Macro to add profiling of the local scope.
#define PROFILE(label) utl::prof _prof(label)

#if PROFILE_LEVEL >= 1
#define PROFILE1(label) PROFILE(label)
#else
//! \brief Macro to add level 1 profiling of the local scope.
#define PROFILE1(label)
#endif

#if PROFILE_LEVEL >= 2
#define PROFILE2(label) PROFILE(label)
#else
//! \brief Macro to add level 2 profiling of the local scope.
#define PROFILE2(label)
#endif

#if PROFILE_LEVEL >= 3
#define PROFILE3(label) PROFILE(label)
#else
//! \brief Macro to add level 3 profiling of the local scope.
#define PROFILE3(label)
#endif

#if PROFILE_LEVEL >= 4
#define PROFILE4(label) PROFILE(label)
#else
//! \brief Macro to add level 4 profiling of the local scope.
#define PROFILE4(label)
#endif

#endif
