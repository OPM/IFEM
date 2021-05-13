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
  \details This template can be instantiated over any type implementing the
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
    preRef = 0;
    refTol = 0.0;
    minFrac = 0.1;
    aStart = beta = -1.0;
    dumpGrid = 0;
    lastRLev = -1;
  }

  //! \brief Empty destructor.
  virtual ~SIMSolverTS() {}

  //! \brief Reads solver data from the specified input file.
  virtual bool read(const char* file)
  {
    if (!this->SIMadmin::read(file))
      return false;
    else if (preRef < 1)
      return true;

    // Perform some pre-refinement to resolve geometric features, etc.
    return this->S1.preRefine(maxRef,preRef,refTol);
  }

  //! \brief Solves the problem up to the final time.
  virtual int solveProblem(char* infile, const char* heading)
  {
    if (this->tp.step == 0)
    {
      // Perform some initial refinement to resolve geometric features, etc.
      if (maxRef > preRef && !this->S1.initialRefine(beta,minFrac,maxRef))
        return 4;

      // Dump the initial (refined) mesh
      if (maxRef > preRef || preRef > 0)
        this->dumpMesh(infile,false);
    }

    // Save the initial FE model (and state) to VTF and HDF5 for visualization
    int geoBlk = 0, nBlock = 0;
    if (!this->saveState(geoBlk,nBlock,true,infile,preRef<1))
      return 2;

    this->printHeading(heading);

    // Pre-adaptive loop, solve for a certain time period on the initial mesh
    while (this->tp.time.t < aStart && this->advanceStep())
      if (!this->S1.solveStep(this->tp))
        return 3;
      else if (!this->saveState(geoBlk,nBlock))
        return 2;

    // Activate printing of time level for result dumps
    SIMSolver<T1>::dumpLog = true;

    // Adaptive loop
    while (!this->tp.finished())
    {
      this->S1.saveState(); // Save current solution state internally
      int tranStep = this->tp.step; // Time step of next solution transfer
      int lastRef = 1, refElms = 0;

      // Prediction cycle loop, disable staggering cycles
      if (maxPred < 0) this->S1.enableStaggering(false);
      for (int iPred = 0; iPred < abs(maxPred) && lastRef > 0; iPred++)
      {
        lastRef = 0;
        if (iPred > 0)
        {
          // Reset time step counter to the last saved state
          this->tp.reset(tranStep);
          // Save current solution state internally for the refined mesh
          this->S1.saveState();
        }

        // Solve for (up to) nPredict steps without saving results on disk
        for (size_t i = 0; i < nPredict && this->advanceStep(); i++)
          if (!this->S1.solveStep(this->tp))
            return 3;

        if (maxPred > 0 && this->S1.stopped())
          break; // Terminated due to other user-criterion
        else if (this->tp.finished())
          break; // Stop time reached

        // Adapt the mesh based on current solution, and
        // transfer the old solution variables onto the new mesh
        IFEM::cout <<"\n  >>> Paused for mesh refinement at time="
                   << this->tp.time.t << std::endl;
        if ((lastRef = this->S1.adaptMesh(beta,minFrac,maxRef)))
          this->dumpMesh(infile,false);

        refElms += lastRef; // Total number of refined elements
      }

      if (lastRef < 0)
        return 4; // Something went wrong, bailing
      else if (refElms == 0 && !this->tp.finished() && !this->S1.stopped())
        IFEM::cout <<"  No refinement, resume from current state"<< std::endl;

      if (lastRef == 0 && maxPred > 0)
      {
        // The mesh is sufficiently refined at this state.
        // Save the current results to VTF and HDF5, and continue.
        // Note: We then miss the results from the preceding nPredict-1 steps.
        if (!this->saveState(geoBlk,nBlock,refElms > 0))
          return 2;
      }
      else
      {
        // Either the mesh was refined or we were doing prediction steps
        // without staggering cycles. Now continue (at most) nForward steps
        // without checking for further refinement at those steps.
        this->tp.reset(tranStep);
        IFEM::cout <<"\n  >>> Resuming simulation from time="<< this->tp.time.t;
        if (lastRef > 0)
        {
          IFEM::cout <<" on the new mesh."<< std::endl;
          if (SIMSolver<T1>::restartAdm)
            lastRLev = SIMSolver<T1>::restartAdm->getTimeLevel()+1;
        }
        else
        {
          IFEM::cout <<" on current mesh."<< std::endl;
          this->S1.restoreState(); // Solution transfer was not performed
        }

        // Solve for each time step up to final time,
        // but only up to nForward steps on this mesh
        this->S1.enableStaggering(true);
        for (size_t j = 0; j < nForward && this->advanceStep(); j++)
          if (!this->S1.solveStep(this->tp))
            return 3;
          else if (!this->saveState(geoBlk,nBlock,j < 1))
            return 2;
      }
    }

    return this->dumpMesh(infile) ? 0 : 2;
  }

  //! \brief Serialize internal state for restarting purposes.
  //! \param data Container for serialized data
  virtual bool serialize(HDF5Restart::SerializeData& data)
  {
    if (lastRLev < 0 && SIMSolver<T1>::restartAdm && this->S1.getRefined())
      lastRLev = SIMSolver<T1>::restartAdm->getTimeLevel()+1;

    if (lastRLev >= 0)
      data[this->S1.getName()+"::basis"] = std::to_string(lastRLev);
    return this->SIMSolver<T1>::serialize(data);
  }

protected:
  using SIMSolver<T1>::parse;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),"adaptive"))
      return this->SIMSolver<T1>::parse(elem);

    utl::getAttribute(elem,"start",aStart);
    utl::getAttribute(elem,"dump",dumpGrid);
    utl::getAttribute(elem,"dumpLog",SIMSolver<T1>::dumpLog);

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
      else if ((value = utl::getValue(child,"prerefine")))
      {
        if (utl::getAttribute(child,"levels",preRef))
          this->S1.parsePreref(child);
        else
          preRef = atoi(value);
        utl::getAttribute(child,"tol",refTol);
      }
      else if (!strcasecmp(child->Value(),"prerefine"))
      {
        utl::getAttribute(child,"levels",preRef);
        utl::getAttribute(child,"tol",refTol);
      }

    if (preRef > 0)
      IFEM::cout <<"\tNumber of pre-refinement cycles: "<< preRef
                 <<"\n\tPre-refinement tolerance: "<< refTol <<"\n";
    IFEM::cout <<"\tNumber of trial steps before refinement: "<< nPredict
               <<"\n\tNumber of steps between each refinement: "<< nForward
               <<"\n\tMax. number of trial cycles at refinement: "<< maxPred
               <<"\n\tMax. number of refinements of an element: "<< maxRef
               <<"\n\tRelative refinement threshold: "<< minFrac;
    if (beta > 0.0)
      IFEM::cout <<"\n\tPercentage of elements to refine: "<< beta;
    if (aStart > 0.0)
      IFEM::cout <<"\n\tMesh adaptation starting at time: "<< aStart;

    IFEM::cout << std::endl;
    return true;
  }

  //! \brief Dumps the refined LR-mesh to file.
  bool dumpMesh(char* infile, bool done = true) const
  {
    if (dumpGrid < 1)
      return true;
    else if (dumpGrid == 1) // Dump the final grid only
      return done ? this->S1.dumpMesh(strcat(strtok(infile,"."),".lr")) : true;
    else if (done)
      return true;

    // Dump grid after each refinement step
    char fileName[128];
    static int gridNo = 1;
    sprintf(fileName,"%s_m%d.lr",strtok(infile,"."),++gridNo);
    return this->S1.dumpMesh(fileName);
  }

private:
  size_t nForward; //!< Number of steps to advance on new mesh
  size_t nPredict; //!< Number of steps to use for prediction
  int    maxPred;  //!< Maximum number of prediction cycles
  int    maxRef;   //!< Number of refinements to allow for a single element
  int    preRef;   //!< Number of pre-refinement cycles
  double refTol;   //!< Pre-refinement threshold
  double minFrac;  //!< Element-level refinement threshold
  double beta;     //!< Percentage of elements to refine
  double aStart;   //!< Time from where to check for mesh adaptation
  char   dumpGrid; //!< Option for mesh output: 0=none, 1=last, 2=all
  int    lastRLev; //!< Restart time level for last mesh refinement
};

#endif
