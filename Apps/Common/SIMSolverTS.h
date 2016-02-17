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
    nForward = nPredict = 1;
    maxRef = 2;
    minFrac = 0.1;
    beta = -1.0;
  }

  //! \brief Empty destructor.
  virtual ~SIMSolverTS() {}

  //! \brief Solves the problem up to the final time.
  virtual int solveProblem(char* infile, DataExporter* exporter = nullptr,
                           const char* heading = nullptr, bool = true)
  {
    // Perform some initial refinement to resolve geometric features, etc.
    if (!this->S1.initialRefine(beta,minFrac,maxRef))
      return 4;

    // Save the initial FE model to VTF file for visualization
    int geoBlk = 0, nBlock = 0;
    if (!this->S1.saveModel(infile,geoBlk,nBlock))
      return 1;

    // Save the initial configuration also
    if (!this->S1.saveStep(this->tp,nBlock))
      return 2;

    this->printHeading(heading);

    // Adaptive loop
    for (int iStep = 1; true; iStep++)
    {
      // Save the current time and solution state internally
      int sStep = this->tp.step;
      this->S1.saveState();

      // Solve for (up to) nPredict steps without storing any results on disk
      bool finished = false;
      for (size_t i = 0; i < nPredict; i++)
        if ((finished = !this->advanceStep()))
          break; // Final time reached, exit time stepping loop
        else if (!this->S1.solveStep(this->tp))
          return 3;

      int newElms = 0;
      if (!finished)
      {
	// Adapt the mesh based on current solution, and
	// transfer the old solution variables onto the new mesh
	IFEM::cout <<"\n  >>> Paused for mesh refinement at time="
                   << this->tp.time.t << std::endl;
	if ((newElms = this->S1.adaptMesh(beta,minFrac,maxRef)) < 0)
	  return 4;
	else if (newElms == 0)
	  IFEM::cout <<"  No refinement, resume from current state"<< std::endl;
      }

      if (newElms == 0)
      {
        // The mesh is unchanged, so just continue after saving the last state.
        // Note: We then miss the results from the preceding nPredict-1 steps.
        if (!this->S1.saveStep(this->tp,nBlock))
          return 2;
        else if (exporter)
          exporter->dumpTimeLevel(&this->tp, iStep < 2);
        if (finished)
          return 0;
        else
          continue;
      }

      // Save the updated FE model to VTF
      if (!this->S1.saveModel(nullptr,geoBlk,nBlock))
        return 1;

      this->tp.reset(sStep);
      IFEM::cout <<"\n  >>> Resuming simulation from time="<< this->tp.time.t
                 <<" on the new mesh."<< std::endl;

      // Solve for each time step up to final time,
      // but only up to nForward steps on this mesh
      for (size_t j = 0; j < nForward; j++)
        if (!this->advanceStep())
          return 0; // Final time reached, we're done
        else if (!this->S1.solveStep(this->tp))
          return 3;
        else if (!this->S1.saveStep(this->tp,nBlock))
          return 2;
        else if (exporter)
          exporter->dumpTimeLevel(&this->tp, j < 1);
    }

    return 0;
  }

protected:
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
      else if ((value = utl::getValue(child,"min_frac")))
        minFrac = atof(value);

    IFEM::cout <<"\tNumber of trial steps before refinement: "<< nPredict
               <<"\n\tNumber of steps between each refinement: "<< nForward
               <<"\n\tMax. number of refinements of an element: "<< maxRef
               <<"\n\tRelative refinement threshold on |c|: "<< minFrac;
    if (beta > 0.0)
      IFEM::cout <<"\n\tPercentage of elements to refine: "<< beta;
    IFEM::cout << std::endl;
    return true;
  }

private:
  size_t nForward; //!< Number of steps to advance on new mesh
  size_t nPredict; //!< Number of steps to use for prediction
  int    maxRef;   //!< Number of refinements to allow for a single element
  double minFrac;  //!< Element-level refinement threshold on |c|
  double beta;     //!< Percentage of elements to refine
};

#endif
