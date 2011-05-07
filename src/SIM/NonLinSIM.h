// $Id$
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
#include "SIMparameters.h"
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

  //! \brief A class for nonlinear solution parameters.
  class SolvePrm : public SIMparameters
  {
  public:
    //! \brief Default constructor.
    SolvePrm() : refNorm(1.0), alpha(1.0) {}

    int    maxit;   //!< Maximum number of iterations
    int    nupdat;  //!< Number of iterations with updated tangent
    double convTol; //!< Relative convergence tolerance
    double divgLim; //!< Relative divergence limit
    double refNorm; //!< Reference energy norm used in convergence checks
    double alpha;   //!< Iteration acceleration parameter (line search)
    double eta;     //!< Line search tolerance
  };

  //! \brief Initializes the solution parameters with values read from file.
  //! \param param Solution algorithm parameters
  //! \param[in] initVal Initial values of the primary solution
  virtual void init(SolvePrm& param, const Vector& initVal = Vector());
  //! \brief Advances the time/load step one step forward.
  virtual bool advanceStep(SolvePrm& param);

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[in] format Format of VTF-file (0=ASCII, 1=BINARY)
  //! \param[in] nViz Number of visualization points over a knot-span
  bool saveModel(char* fileName, int format, int* nViz);

  //! \brief Sets the initial guess in the Newton-Raphson iterations.
  //! \param value The initial guess to use
  void setInitialGuess(const Vector& value) { solution.front() = value; }

  //! \brief Solves the nonlinear equations by Newton-Raphson iterations.
  //! \param param Solution algorithm parameters
  //! \param[in] mode Solution mode to use for this step
  //! \param[in] compName Solution name to be used in the norm output
  //! \param[in] energyNorm If \e true, integrate energy norm of the solution
  bool solveStep(SolvePrm& param, SIM::SolutionMode mode = SIM::STATIC,
		 const char* compName = "displacement",
		 bool energyNorm = false);

  //! \brief Enum describing the norm used for convergence checks.
  enum CNORM { L2, ENERGY };

  //! \brief Sets the current norm in used in iteration error estimates.
  //! param[in] norm The wanted error norm
  void setErrorNorm(CNORM norm) { iteNorm = norm; }

  //! \brief Computes and prints some solution norm quantities.
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param[in] compName Solution name to be used in the norm output
  //! \param[in] energyNorm If \e true, integrate energy norm of the solution
  bool solutionNorms(const TimeDomain& time,
		     const char* compName, bool energyNorm);

  //! \brief Saves the converged results to VTF file of a given load/time step.
  //! \param[in] iStep Load/time step identifier
  //! \param[in] time Current time/load parameter
  //! \param[in] nViz Number of visualization points over each knot-span
  //! \param[in] psolOnly If \e true, skip secondary solution field output
  //! \param[in] vecName Optional name of primary solution vector field
  bool saveStep(int iStep, double time, int* nViz, bool psolOnly = false,
		const char* vecName = 0);

  //! \brief Dumps the primary solution for inspection.
  //! \param[in] iStep Load/time step identifier
  //! \param[in] time Current time/load parameter
  //! \param[in] os The output stream to write the solution to
  //! \param[in] withID If \e true, write node ID and coordinates too
  void dumpStep(int iStep, double time, std::ostream& os,
		bool withID = true) const;

  //! \brief Dumps solution variables at user-defined points.
  //! \param[in] time Current time/load parameter
  //! \param[in] os The output stream to write the solution to
  //! \param[in] precision Number of digits after the decimal point
  bool dumpResults(double time, std::ostream& os, int precision = 3) const;

  //! \brief Returns a const reference to current solution vector.
  const Vector& getSolution() const { return solution.front(); }

  //! \brief Projects the secondary solution onto the spline basis.
  //! \details The secondary solution, defined through the Integrand class,
  //! is projected onto the spline basis to obtain the control point values
  //! of the secondary solution. These values then overwrite the current
  //! primary solution vector.
  bool project();

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
  Vector  linsol;   //!< Linear solution vector
  Vector  residual; //!< Residual force vector
  Vectors solution; //!< Total solution vectors

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
  CNORM  iteNorm;   //!< The norm type used to measure the residual

  // Post-processing attributes
  int    nBlock; //!< Running VTF result block counter
  Vector gNorm;  //!< Global norms
  Matrix eNorm;  //!< Element norms
};

#endif
