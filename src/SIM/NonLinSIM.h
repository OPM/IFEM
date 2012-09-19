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
#include "MatVec.h"

class SIMbase;
class TimeStep;
class TimeDomain;


/*!
  \brief Nonlinear solution driver for isogeometric FEM simulators.
  \details This class contains data and methods for computing the nonlinear
  solution to a FE problem based on splines/NURBS basis functions, through
  Newton-Raphson iterations.
*/

class NonLinSIM : public SIMinput
{
public:
  //! \brief Convergence status enum.
  enum ConvStatus { FAILURE, OK, SLOW, CONVERGED, DIVERGED };
  //! \brief Enum describing the norm used for convergence checks.
  enum CNORM { L2, ENERGY };
  //! \brief Enum describing reference norm options.
  enum NormOp { MAX, ALL };

  //! \brief The constructor initializes default solution parameters.
  //! \param sim Pointer to the spline FE model
  //! \param[in] n Which type of iteration norm to use in convergence checks
  NonLinSIM(SIMbase* sim = 0, CNORM n = ENERGY);
  //! \brief The destructor frees the dynamically allocated FE model object.
  virtual ~NonLinSIM();

  //! \brief Initializes the primary solution vectors.
  //! \param[in] initVal Initial values of the primary solution
  virtual void init(const RealArray& initVal = RealArray());

  //! \brief Sets the initial guess in the Newton-Raphson iterations.
  //! \param[in] value Initial values of the primary solution
  void setInitialGuess(const RealArray& value);

  //! \brief Advances the time/load step one step forward.
  //! \param param Time stepping parameters
  //! \param[in] updateTime If \e false, the time parameters are not incremented
  virtual bool advanceStep(TimeStep& param, bool updateTime = true);

  //! \brief Solves the nonlinear equations by Newton-Raphson iterations.
  //! \param[in] zero_tolerance Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  ConvStatus solve(double zero_tolerance = 1.0e-8, std::streamsize outPrec = 0);

  //! \brief Solves the nonlinear equations by Newton-Raphson iterations.
  //! \param param Time stepping parameters
  //! \param[in] mode Solution mode to use for this step
  //! \param[in] energyNorm If \e true, integrate energy norm of the solution
  //! \param[in] zero_tolerance Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  ConvStatus solveStep(TimeStep& param, SIM::SolutionMode mode = SIM::STATIC,
                       bool energyNorm = false, double zero_tolerance = 1.0e-8,
                       std::streamsize outPrec = 0);

  //! \brief Computes and prints some solution norm quantities.
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param[in] energyNorm If \e true, integrate energy norm of the solution
  //! \param[in] zero_tolerance Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  virtual bool solutionNorms(const TimeDomain& time,
			     bool energyNorm = false,
			     double zero_tolerance = 1.0e-8,
			     std::streamsize outPrec = 0);

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  bool saveModel(char* fileName);

  //! \brief Saves the converged results to VTF file of a given load/time step.
  //! \param[in] iStep Load/time step identifier
  //! \param[in] time Current time/load parameter
  //! \param[in] psolOnly If \e true, skip secondary solution field output
  //! \param[in] vecName Optional name of primary solution vector field
  virtual bool saveStep(int iStep, double time,
			bool psolOnly = false, const char* vecName = 0);

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
  void dumpResults(double time, std::ostream& os,
                   std::streamsize precision = 3) const;

  //! \brief Returns a const reference to current solution vector.
  const Vector& getSolution(int i = 0) const { return solution[i]; }

protected:
  //! \brief Checks whether the nonlinear iterations have converged or diverged.
  virtual ConvStatus checkConvergence(TimeStep& param);
  //! \brief Updates configuration variables (solution vector) in an iteration.
  virtual bool updateConfiguration(TimeStep& param);
  //! \brief Performs line search to accelerate convergence.
  virtual bool lineSearch(TimeStep& param);

public:
  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);

  //! \brief Parses a data section from an XML document.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Returns a list of prioritized XML-tags.
  virtual const char** getPrioritizedTags() const;

protected:
  SIMbase* model;    //!< The isogeometric FE model
  Vectors  solution; //!< Primary solution vectors of the last increments
  Vector   linsol;   //!< Linear solution vector
  Vector   residual; //!< Residual force vector

  // Nonlinear solution algorithm parameters
  CNORM  iteNorm; //!< The norm type used to measure the residual
  NormOp refNopt; //!< Reference norm option
  double refNorm; //!< Reference norm value used in convergence checks
  double convTol; //!< Relative convergence tolerance
  double divgLim; //!< Relative divergence limit
  double eta;     //!< Line search tolerance
  double alpha;   //!< Iteration acceleration parameter (for line search)
  int    maxit;   //!< Maximum number of iterations in a time/load step
  int    nupdat;  //!< Number of iterations with updated tangent
  int    prnSlow; //!< How many DOFs to print out on slow convergence

  // Post-processing attributes
  int    geoBlk; //!< Running VTF geometry block counter
  int    nBlock; //!< Running VTF result block counter
  Matrix eNorm;  //!< Element norm values
};

#endif
