// $Id$
//==============================================================================
//!
//! \file SIMSolverTS.h
//!
//! \date Feb 1 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Time-slab adaptive SIM solver class template.
//!
//==============================================================================

#ifndef _SIM_SOLVER_TS_H_
#define _SIM_SOLVER_TS_H_

#include "SIMSolver.h"


/*!
  \brief Template class for time-slab adaptive simulator drivers.
  \details This template can be instanciated over any type implementing the
  ISolver interface. It provides a time stepping loop with data output,
  with mesh refinements at fixed time intervals. The refinement is based on the
  solution state a given number of increments ahead of the refinement time.
  The solution is then recomputed on the refined mesh after each refinement,
  starting from the time at which the refinement took place.
*/

template<class T1> class SIMSolverTS : public SIMSolver<T1>
{
public:
  //! \brief The constructor forwards to the parent class constructor.
  SIMSolverTS(T1& s1) : SIMSolver<T1>(s1)
  {
    nForward = nPredict = maxPred = 1;
    maxRef = 2;
    minFrac = 0.1;
    beta = -1.0;
  }

  //! \brief Empty destructor.
  virtual ~SIMSolverTS() {}

  //! \brief Solves the problem up to the final time.
  virtual int solveProblem(char* infile, const char* heading, bool saveInit)
  {
    // Perform some initial refinement to resolve geometric features, etc.
    if (!this->S1.initialRefine(beta,minFrac,maxRef))
      return 4;

    // Save the initial FE model (and state) to VTF and HDF5 for visualization
    int geoBlk = 0, nBlock = 0;
    if (!this->saveState(geoBlk,nBlock,true,infile,saveInit))
      return 2;

    this->printHeading(heading);

    // Adaptive loop
    for (int iStep = 1; true; iStep++)
    {
      if (iStep > 1)
        this->S1.saveState(); // Save current solution state internally
      int tranStep = this->tp.step; // Time step of next solution transer
      int lastRef = 0, refElms = 0;

      // Prediction cycle loop
      for (int iPred = 0; iPred < maxPred; iPred++)
      {
        if (iPred > 0)
        {
          // Reset time step counter to the last saved state
          this->tp.reset(tranStep);
          // Save current solution state internally for the refined mesh
          this->S1.saveState();
        }

        // Solve for (up to) nPredict steps without saving results on disk
        for (size_t i = 0; i < nPredict; i++)
          if (!this->advanceStep()) // Final time reached
            // Save updated FE model if mesh has been refined
            // Save current results to VTF and HDF5 and then exit
            return this->saveState(geoBlk,nBlock,refElms > 0) ? 0 : 2;
          else if (!this->S1.solveStep(this->tp))
            return 3;

        // Adapt the mesh based on current solution, and
        // transfer the old solution variables onto the new mesh
        IFEM::cout <<"\n  >>> Paused for mesh refinement at time="
                   << this->tp.time.t << std::endl;
        lastRef = this->S1.adaptMesh(beta,minFrac,maxRef);
        if (lastRef < 0)
          return 4; // Something went wrong, bailing
        else if (lastRef > 0)
          refElms += lastRef; // Total number of refined elements
        else // No new mesh refinement, exit the prediction cycles
          break;
      }

      if (refElms == 0)
        IFEM::cout <<"  No refinement, resume from current state"<< std::endl;

      if (lastRef == 0)
      {
        // The mesh is sufficiently refined at this state.
        // Save the current results to VTF and HDF5, and continue.
        // Note: We then miss the results from the preceding nPredict-1 steps.
        if (!this->saveState(geoBlk,nBlock,refElms > 0))
          return 2;
      }
      else
      {
        // The mesh was refined, and we continue (at most) nForward steps
        // without checking for further refinement at this step
        this->tp.reset(tranStep);
        IFEM::cout <<"\n  >>> Resuming simulation from time="<< this->tp.time.t
                   <<" on the new mesh."<< std::endl;

        // Solve for each time step up to final time,
        // but only up to nForward steps on this mesh
        for (size_t j = 0; j < nForward; j++)
          if (!this->advanceStep())
            return 0; // Final time reached, we're done
          else if (!this->S1.solveStep(this->tp))
            return 3;
          else if (!this->saveState(geoBlk,nBlock,j < 1))
            return 2;
      }
    }
  }

protected:
  using SIMSolver<T1>::parse;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),"adaptive"))
      return this->SIMSolver<T1>::parse(elem);

    const char* value = nullptr;
    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      if ((value = utl::getValue(child,"forward_steps")))
        nForward = atoi(value);
      else if ((value = utl::getValue(child,"predict_steps")))
        nPredict = atoi(value);
      else if ((value = utl::getValue(child,"beta")))
        beta = atof(value);
      else if ((value = utl::getValue(child,"nrefinements")))
        maxRef = atoi(value);
      else if ((value = utl::getValue(child,"max_prediction")))
        maxPred = atoi(value);
      else if ((value = utl::getValue(child,"min_frac")))
        minFrac = atof(value);

    IFEM::cout <<"\tNumber of trial steps before refinement: "<< nPredict
               <<"\n\tNumber of steps between each refinement: "<< nForward
               <<"\n\tMax. number of trial cycles at refinement: "<< maxPred
               <<"\n\tMax. number of refinements of an element: "<< maxRef
               <<"\n\tRelative refinement threshold: "<< minFrac;
    if (beta > 0.0)
      IFEM::cout <<"\n\tPercentage of elements to refine: "<< beta;
    IFEM::cout << std::endl;
    return true;
  }

private:
  size_t nForward; //!< Number of steps to advance on new mesh
  size_t nPredict; //!< Number of steps to use for prediction
  int    maxPred;  //!< Maximum number of prediction cycles
  int    maxRef;   //!< Number of refinements to allow for a single element
  double minFrac;  //!< Element-level refinement threshold
  double beta;     //!< Percentage of elements to refine
};

#endif
