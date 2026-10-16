// $Id$
//==============================================================================
//!
//! \file SIMStat.h
//!
//! \date Oct 16 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Iterates a time stepping simulator towards its stationary state.
//!
//==============================================================================

#ifndef _SIM_STAT_H
#define _SIM_STAT_H

#include "AndersonAcceleration.h"

#include "IFEM.h"
#include "MatVec.h"
#include "Profiler.h"
#include "SAM.h"
#include "TimeStep.h"
#include "Utilities.h"

#include <tinyxml2.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>


/*!
  \brief Describes which DOFs to measure a solution vector using.
*/

struct ConvGroup
{
  std::string name;      //!< What the simulator calls this part of the solution
  std::vector<int> dofs; //!< Nodal degrees of freedom it covers
};


/*!
  \brief Iterates a time stepping simulator towards its stationary state.

  \details A time stepping scheme reaches the stationary state by being
  iterated in pseudo time until the solution stops changing. The time
  derivative is identically zero at such a fixed point.
  The iteration is an accelerated fixed point iteration.

  The simulator underneath has to provide:

  - \a solveTimeStep(tp), one step of the scheme without the post-processing
    that follows a converged step.
  - \a postSolve(tp) for post-processing.
  - \a convergenceGroups(), describing the parts of the solution the change
    should be measured on,
  - \a inputContext, the XML tag the simulator reads its own options from.

  and, from SIMsolution, at least two solution vectors, this step's and the
  one before it.
*/

template<class SIM>
class SIMStat : public SIM
{
public:
  using SIM::SIM; // Import constructors from base class

  //! \brief Iterates the scheme towards the stationary state.
  //! \param tp Time stepping parameters
  bool solveStep(TimeStep& tp)
  {
    PROFILE1("SIMStat::solveStep");

    const std::vector<ConvGroup> groups = this->convergenceGroups();
    if (groups.empty()) {
      std::cerr <<" *** SIMStat: The simulator named no part of its solution"
                <<" to measure the change on."<< std::endl;
      return false;
    }

    // Named in the input, or the one the simulator puts first
    size_t iGroup = 0;
    if (!cNorm.empty() && cNorm != "all") {
      const auto it = std::find_if(groups.begin(), groups.end(),
                                   [this](const ConvGroup& g)
                                   { return g.name == cNorm; });
      if (it == groups.end()) {
        std::cerr <<" *** SIMStat: There is no \""<< cNorm <<"\" to measure"
                  <<" the change on. The simulator offers";
        for (const ConvGroup& g : groups)
          std::cerr <<" \""<< g.name <<"\"";
        std::cerr <<" and \"all\"."<< std::endl;
        return false;
      }
      iGroup = it - groups.begin();
    }

    if (this->getSolutions().size() < 2) {
      std::cerr <<" *** SIMStat: The simulator keeps "
                << this->getSolutions().size() <<" solution vector, and a"
                <<" time stepping scheme needs the step before it as well."
                << std::endl;
      return false;
    }

    const int outerIter = tp.iter;
    tp.time.dt = tp.time.dtn = pseudoDt;

    // The per-iteration solution summary would drown the convergence history,
    // so it is printed once, after the last iteration
    const int msgL = SIM::msgLevel;
    SIM::msgLevel = std::min(msgL,0);

    bool converged = false;
    for (tp.iter = 0; tp.iter <= maxIter && !converged; tp.iter++)
    {
      // The point this step is taken from, kept rather than read back from
      // index one afterwards. The first iteration does not advance, so what
      // index one holds then is whatever the caller left there and not the
      // point the iteration starts at
      const Vector was(this->getSolutions().front());

      if ((tp.iter > 0 && !this->advanceStep(tp)) || !this->solveTimeStep(tp)) {
        SIM::msgLevel = msgL;
        return false;
      }

      // Change over the pseudo time step, field by field. The solution of
      // this simulator is read directly, since the adaptive loop redirects
      // index zero to a vector of its own which is only synchronized after
      // the solve.
      const Vector& now = this->getSolutions().front();

      // The norms go through the SAM, which owns each shared degree of
      // freedom once and reduces over the processes. A norm summed here
      // would be a different number on every process, and the processes
      // would then leave the loop at different iterations
      const SAM& sam = *this->getSAM();
      std::vector<double> dX(groups.size()), nX(groups.size());
      Vector dv(now.size()), nv(now.size());
      for (size_t i = 0; i < groups.size(); i++) {
        dv.fill(0.0);
        nv.fill(0.0);
        for (int d : groups[i].dofs) {
          dv(d) = now(d) - was(d);
          nv(d) = was(d);
        }
        dX[i] = sam.norm2(dv,'A');
        nX[i] = sam.norm2(nv,'A');
      }

      double diff, rel;
      if (cNorm == "all") {
        double d2 = 0.0, n2 = 0.0;
        for (size_t i = 0; i < groups.size(); i++) {
          d2 += dX[i]*dX[i];
          n2 += nX[i]*nX[i];
        }
        diff = sqrt(d2); rel = sqrt(n2);
      }
      else {
        diff = dX[iGroup]; rel = nX[iGroup];
      }
      if (rel == 0.0)
        rel = 1.0;
      diff /= rel;

      if (msgL >= 0) {
        IFEM::cout <<"  iter="<< tp.iter <<" conv="<< diff;
        for (size_t i = 0; i < groups.size(); i++)
          IFEM::cout <<" "<< groups[i].name <<"_abs="<< dX[i]
                     <<" "<< groups[i].name <<"_rel="<< nX[i];
        IFEM::cout << std::endl;
      }

      converged = diff < rTol && tp.iter > 0;

      // The change measured above is the residual of the fixed point
      // iteration, so it is taken before the point is moved by the
      // acceleration
      if (!converged && anderson.getDepth() > 0)
        anderson.accelerate(sam, was, this->theSolutions().front());
    }

    SIM::msgLevel = msgL;
    tp.iter = outerIter;

    if (!converged)
      IFEM::cout <<"  ** The stationary state was not reached in "<< maxIter
                 <<" iterations."<< std::endl;

    if (SIM::msgLevel >= 1)
      this->printSolutionSummary(this->getSolutions().front(), 1,
                                 groups[iGroup].name.c_str(), 0);

    return this->postSolve(tp);
  }

  using SIM::init;
  //! \brief Initializes the simulator for the mesh it now has.
  //! \param[in] tp Time stepping parameters
  //! \details The points of the acceleration history belong to the mesh they
  //! were computed on, so they are dropped here.
  bool init(const TimeStep& tp)
  {
    anderson.clear();
    return this->SIM::init(tp);
  }

protected:
  //! \brief Parses a data section from an XML element.
  //! \param[in] elem The XML element to parse
  bool parse(const tinyxml2::XMLElement* elem) override
  {
    if (strcasecmp(elem->Value(),this->inputContext.c_str()))
      return this->SIM::parse(elem);

    utl::getAttribute(elem,"stationary_norm",cNorm);
    utl::getAttribute(elem,"stationary_tol",rTol);
    utl::getAttribute(elem,"dt",pseudoDt);
    utl::getAttribute(elem,"maxits",maxIter);

    int aDepth = 0;
    if (utl::getAttribute(elem,"anderson",aDepth)) {
      if (aDepth < 0) {
        std::cerr <<" *** SIMStat: A negative Anderson depth, "<< aDepth
                  <<", was asked for."<< std::endl;
        return false;
      }
      anderson.setDepth(aDepth);
    }

    // Which group this names is checked where the groups exist, the model
    // not having been preprocessed yet
    IFEM::cout <<"\tIterating to the stationary state on the change in "
               << (cNorm.empty() ? std::string("the primary field") :
                   cNorm == "all" ? std::string("every field") : "the " + cNorm)
               <<", tolerance "<< rTol
               <<",\n\t  at most "<< maxIter <<" iterations of step "
               << pseudoDt << std::endl;
    if (anderson.getDepth() > 0)
      IFEM::cout <<"\tAccelerating the iteration over its last "
                 << anderson.getDepth() <<" steps."<< std::endl;
    else
      IFEM::cout <<"\tNot accelerating the iteration."<< std::endl;

    return this->SIM::parse(elem);
  }

  double rTol = 1.0e-6;  //!< Relative change accepted as stationary
  double pseudoDt = 1.0; //!< Step size of the pseudo time stepping
  int maxIter = 1000;    //!< Maximum number of iterations

  //! \brief Part of the solution to measure the change on.
  //! \details Empty for the one the simulator puts first, which is the one
  //! it considers primary, or "all" for every field at once.
  std::string cNorm;

  //! \brief Acceleration of the iteration.
  //! \details On by default. Without it the iteration is a relaxed
  //! Richardson iteration for the Schur complement and needs thousands of
  //! steps, which makes the driver of little practical use.
  AndersonAcceleration anderson;
};

#endif
