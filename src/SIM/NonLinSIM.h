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

#ifndef _NON_LIN_SIM_H
#define _NON_LIN_SIM_H

#include "MultiStepSIM.h"


/*!
  \brief Nonlinear quasi-static solution driver for isogeometric FEM simulators.
  \details This class contains data and methods for computing the nonlinear
  solution to a quasi-static FE problem based on splines/NURBS basis functions,
  through Newton-Raphson iterations.
*/

class NonLinSIM : public MultiStepSIM
{
public:
  //! \brief Enum describing the norm used for convergence checks.
  enum CNORM { NONE, L2, L2SOL, ENERGY };

  //! \brief The constructor initializes default solution parameters.
  //! \param sim Pointer to the spline FE model
  //! \param[in] n Which type of iteration norm to use in convergence checks
  NonLinSIM(SIMbase& sim, CNORM n = ENERGY);
  //! \brief The destructor prints out the slow-converging nodes, if any.
  virtual ~NonLinSIM();

  //! \brief Defines which type of iteration norm to use in convergence checks.
  void setConvNorm(CNORM n) { if ((iteNorm = n) == NONE) fromIni = true; }

  //! \brief Initializes the primary solution vectors.
  //! \param[in] nSol Number of consequtive solutions stored in core
  //! \param[in] initVal Initial values of the primary solution
  void init(size_t nSol, const RealArray& initVal);

  //! \brief Initializes some integration parameters for the integrand.
  virtual void initPrm();

  //! \brief Advances the load step one step forward.
  //! \param param Time stepping parameters
  //! \param[in] updateTime If \e false, the time parameters are not incremented
  virtual bool advanceStep(TimeStep& param, bool updateTime = true);

  //! \brief Solves the nonlinear equations by Newton-Raphson iterations.
  //! \param[in] zero_toler Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  SIM::ConvStatus solve(double zero_toler = 1.0e-8,
                        std::streamsize outPrec = 0);

  //! \brief Solves the nonlinear equations by Newton-Raphson iterations.
  //! \param param Time stepping parameters
  //! \param[in] mode Solution mode to use for this step
  //! \param[in] zero_tolerance Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  virtual SIM::ConvStatus solveStep(TimeStep& param,
                                    SIM::SolutionMode mode = SIM::STATIC,
                                    double zero_tolerance = 1.0e-8,
                                    std::streamsize outPrec = 0);

  //! \brief Solves the linearized system of current iteration.
  //! \param[in] param Time stepping parameters
  SIM::ConvStatus solveIteration(TimeStep& param);

  //! \brief Returns the maximum number of iterations.
  int getMaxit() const { return maxit; }

  //! \brief Returns whether this solution driver is linear or not.
  virtual bool isLinear() const { return iteNorm == NONE; }

protected:
  //! \brief Prints out the worst DOFs when slow convergence is detected.
  //! \param os The output stream to print to
  //! \param[in] eps Only print DOF values larger than this tolerance
  void printWorst(utl::LogStream& os, double eps);

  //! \brief Checks whether the nonlinear iterations have converged or diverged.
  virtual SIM::ConvStatus checkConvergence(TimeStep& param);
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

protected:
  bool   fromIni; //!< If \e true, always solve from initial configuration
  CNORM  iteNorm; //!< The norm type used to measure the residual
  double rTol;    //!< Relative convergence tolerance
  double aTol;    //!< Absolute convergence tolerance
  double divgLim; //!< Relative divergence limit
  double eta;     //!< Line search tolerance
  double alpha;   //!< Iteration acceleration parameter (for line search)
  double alphaO;  //!< Final line search acceleration scaling (for output only)
  int    maxit;   //!< Maximum number of iterations in a load step
  int    maxIncr; //!< Maximum number of iterations with increasing norm
  int    nupdat;  //!< Number of iterations with updated tangent
  int    prnSlow; //!< How many DOFs to print out on slow convergence

  std::map<int,int> slowNodes; //!< Nodes for which slow convergence is detected

public:
  static const char* inputContext; //!< Input file context for solver parameters
};

#endif
