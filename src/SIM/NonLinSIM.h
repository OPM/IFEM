// $Id: NonLinSIM.h,v 1.15 2011-01-28 15:28:54 kmo Exp $
//==============================================================================
//!
//! \file NonLinSIM.h
//!
//! \date Jun 1 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Nonlinear solution driver for isogeometric FEM simulators.
//!
//==============================================================================

#ifndef _NONLIN_SIM_H
#define _NONLIN_SIM_H

#include "SIMinput.h"
#include "SIMenums.h"
#include "TimeDomain.h"
#include "MatVec.h"

class SIMbase;


/*!
  \brief Nonlinear solution driver for isogeometric FEM simulators.
  \details This class contains data and methods for computing the nonlinear
  solution to a FE problem based on splines/NURBS basis functions, through
  Newton-Raphson iterations.
*/

class NonLinSIM : public SIMinput
{
public:
  //! \brief The constructor initialized default solution parameters.
  //! \param sim Pointer to the spline FE model
  NonLinSIM(SIMbase* sim = 0);
  //! \brief The destructor frees the dynamically allocated FE model object.
  virtual ~NonLinSIM();

  /*!
    \brief A struct for nonlinear solution parameters.
  */
  struct SolvePrm
  {
    int  step;        //!< Load/time step counter
    int& iter;        //!< Iteration counter
    int  maxit;       //!< Maximum number of iterations
    int  nupdat;      //!< Number of iterations with updated tangent
    double startTime; //!< Start (pseudo)time of simulation
    double stopTime;  //!< Stop (pseudo)time of simulation
    RealArray tInc;   //!< Time (or pseudo time) increment size(s)
    double convTol;   //!< Relative convergence tolerance
    double divgLim;   //!< Relative divergence limit
    double refNorm;   //!< Reference energy norm used in convergence checks
    double alpha;     //!< Iteration acceleration parameter (line search)
    double eta;       //!< Line searce tolerance
    TimeDomain time;  //!< Time domain data

    //! \brief The constructor initializes the counters to zero.
    SolvePrm() : iter(time.it) { step = maxit = 0; refNorm = alpha = 1.0; }
    //! \brief Returns \e true if the simulation consists of several load steps.
    bool multiSteps() const;
    //! \brief Increments the time/load step.
    //! \return \e true, if we have reached the end of the simulation.
    bool increment();
  };

  //! \brief Initializes solution parameters object with values read from file.
  virtual void init(SolvePrm& param);
  //! \brief Advances the time/load step one step forward.
  virtual bool advanceStep(SolvePrm& param);

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[in] format Format of VTF-file (0=ASCII, 1=BINARY)
  //! \param[in] nViz Number of visualization points over a knot-span
  bool saveModel(char* fileName, int format, int* nViz);

  //! \brief Solves the nonlinear equations by Newton-Raphson iterations.
  //! \param param Solution algorithm parameters
  //! \param[in] mode Solution mode to use for this step
  bool solveStep(SolvePrm& param, SIM::SolutionMode mode = SIM::STATIC);

  //! \brief Computes and prints some solution norm quantities.
  //! \param[in] time Parameters for nonlinear/time-dependent simulations.
  //! \param[in] compName Solution name to be used in the norm output
  bool solutionNorms(const TimeDomain& time,
		     const char* compName = "displacement");

  //! \brief Saves the converged results to VTF file of a given load/time step.
  //! \param[in] iStep Load/time step identifier
  //! \param[in] time Current time/load parameter
  //! \param[in] nViz Number of visualization points over each knot-span
  //! \param[in] psolOnly If \e true, skip secondary solution field output
  bool saveStep(int iStep, double time, int* nViz, bool psolOnly = false);

  //! \brief Dumps the primary solution to given stream for inspection.
  //! \param[in] iStep Load/time step identifier
  //! \param[in] time Current time/load parameter
  //! \param[in] os The output streeam to write to.
  //! \param[in] withID If \e true, write node ID and coordinates too
  void dumpStep(int iStep, double time, std::ostream& os,
		bool withID = true) const;

protected:
  //! \brief Convergence status enum.
  enum ConvStatus { NONE, CONVERGED, DIVERGED };

  //! \brief Perform line search to accelerate convergence.
  virtual bool lineSearch(SolvePrm& param);
  //! \brief Checks whether the nonlinear iterations have converged or diverged.
  virtual ConvStatus checkConvergence(SolvePrm& param);
  //! \brief Updates configuration variables (solution vector) in an iteration.
  virtual bool updateConfiguration(SolvePrm& param);

  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);

private:
  SIMbase* model; //!< The isogeometric FE model

protected:
  Vector              linsol;   //!< Linear solution vector
  Vector              residual; //!< Residual force vector
  std::vector<Vector> solution; //!< Total solution vectors

private:
  // Nonlinear solution algorithm parameters
  double startTime; //!< Start time of the simulation
  double stopTime;  //!< Stop time of the simulation
  RealArray tInc;   //!< Time increment size(s)
  double convTol;   //!< Relative convergence tolerance
  double divgLim;   //!< Relative divergence limit
  double eta;       //!< Line search tolerance
  int    maxit;     //!< Maximum number of iterations in a time/load step
  int    nupdat;    //!< Number of iterations with updated tangent

  // Post-processing attributes
  int    nBlock; //!< Running VTF result block counter
  Vector gNorm;  //!< Global norms
  Matrix eNorm;  //!< Element norms
};

#endif
