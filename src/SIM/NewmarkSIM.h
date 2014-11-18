// $Id$
//==============================================================================
//!
//! \file NewmarkSIM.h
//!
//! \date Jul 4 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Newmark solution driver for isogeometric dynamic FEM simulators.
//!
//==============================================================================

#ifndef _NEWMARK_SIM_H
#define _NEWMARK_SIM_H

#include "MultiStepSIM.h"


/*!
  \brief Newmark-based solution driver for dynamic isogeometric FEM simulators.
*/

class NewmarkSIM : public MultiStepSIM
{
public:
  //! \brief The constructor initializes default solution parameters.
  NewmarkSIM(SIMbase& sim);
  //! \brief Empty destructor.
  virtual ~NewmarkSIM() {}

  using MultiStepSIM::parse;
  //! \brief Parses a data section from an XML document.
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Prints out problem-specific data to the given stream.
  virtual void printProblem(std::ostream& os) const;

  //! \brief Initializes primary solution vectors and integration parameters.
  virtual bool initSol(size_t nSol = 3);

  //! \brief Advances the time step one step forward.
  //! \param param Time stepping parameters
  //! \param[in] updateTime If \e false, the time parameters are not incremented
  virtual bool advanceStep(TimeStep& param, bool updateTime = true);

  //! \brief Solves the dynamic equations by a predictor/multi-corrector method.
  //! \param param Time stepping parameters
  //! \param[in] zero_tolerance Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  virtual SIM::ConvStatus solveStep(TimeStep& param,
                                    SIM::SolutionMode = SIM::STATIC,
                                    double zero_tolerance = 1.0e-8,
                                    std::streamsize outPrec = 0);

protected:
  //! \brief Computes and prints some solution norm quantities.
  //! \param[in] zero_tolerance Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  virtual bool solutionNorms(double zero_tolerance = 1.0e-8,
                             std::streamsize outPrec = 0);

  //! \brief Checks whether the corrector iterations have converged or diverged.
  SIM::ConvStatus checkConvergence(TimeStep& param);
  //! \brief Calculates predicted velocities and accelerations.
  virtual bool predictStep(TimeStep& param);
  //! \brief Updates configuration variables (solution vector) in an iteration.
  virtual bool correctStep(TimeStep& param, bool = false);
  //! \brief Finalizes the right-hand-side vector on the system level.
  virtual void finalizeRHSvector(bool) {}

public:
  //! \brief Returns a const reference to current velocity vector.
  const Vector& getVelocity() const { return solution[solution.size()-2]; }
  //! \brief Returns a const reference to current velocity vector.
  const Vector& getAcceleration() const { return solution[solution.size()-1]; }

  //! \brief Dumps solution variables at user-defined points.
  //! \param[in] time Current time/load parameter
  //! \param[in] os The output stream to write the solution to
  //! \param[in] precision Number of digits after the decimal point
  //! \param[in] formatted If \e false, write all result points on a single line
  virtual void dumpResults(double time, std::ostream& os,
                           std::streamsize precision = 3,
                           bool formatted = true) const;

protected:
  // Time integration parameters
  double alpha1; //!< Mass-proportional damping parameter
  double alpha2; //!< Stiffness-proportional damping parameter
  double beta;   //!< Newmark time integration parameter
  double gamma;  //!< Newmark time integration parameter

  // Solution algorithm parameters
  char   predictor; //!< Predictor type flag
  char   rotUpd;    //!< Option for how to update of nodal rotations
  int    maxit;     //!< Maximum number of iterations in a time step
  double convTol;   //!< Convergence tolerance
  double divgLim;   //!< Relative divergence limit
  unsigned short int cNorm; //!< Option for which convergence norm to use
};

#endif
