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
  explicit NewmarkSIM(SIMbase& sim);

  using MultiStepSIM::parse;
  //! \brief Parses a data section from an XML document.
  virtual bool parse(const tinyxml2::XMLElement* elem);

  //! \brief Prints out problem-specific data to the log stream.
  //! \param[in] stopInputTimer If \e true, stop file input timer before print
  virtual void printProblem(bool stopInputTimer = false) const;

  //! \brief Initializes time integration parameters for the integrand.
  virtual void initPrm();

  //! \brief Calculates initial accelerations.
  bool initAcc(double zero_tolerance = 1.0e-8, std::streamsize outPrec = 0);

  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep& param, bool updateTime = true);

  //! \brief Solves the dynamic equations by a predictor/multi-corrector method.
  //! \param param Time stepping parameters
  //! \param[in] zero_tolerance Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  virtual SIM::ConvStatus solveStep(TimeStep& param,
                                    SIM::SolutionMode = SIM::STATIC,
                                    double zero_tolerance = 1.0e-8,
                                    std::streamsize outPrec = 0);

  //! \brief Solves the linearized system of current iteration.
  //! \param[in] param Time stepping parameters
  SIM::ConvStatus solveIteration(TimeStep& param);

  //! \brief Returns the maximum number of iterations.
  int getMaxit() const { return maxit; }

  //! \brief Returns whether this is a dynamic solution driver or not.
  virtual bool isDynamic() const { return true; }

protected:
  //! \brief Computes and prints some solution norm quantities.
  //! \param[in] zero_tolerance Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  virtual bool solutionNorms(const TimeDomain&, double zero_tolerance = 1.0e-8,
                             std::streamsize outPrec = 0);

  //! \brief Checks whether the corrector iterations have converged or diverged.
  virtual SIM::ConvStatus checkConvergence(TimeStep& param);
  //! \brief Calculates predicted velocities and accelerations.
  virtual bool predictStep(TimeStep& param);
  //! \brief Updates configuration variables (solution vector) in an iteration.
  virtual bool correctStep(TimeStep& param, bool = false);
  //! \brief Finalizes the right-hand-side vector on the system level.
  virtual void finalizeRHSvector(bool);

public:
  //! \brief Returns a const reference to current velocity vector.
  const Vector& getVelocity() const { return solution[1]; }
  //! \brief Returns a const reference to current velocity vector.
  const Vector& getAcceleration() const { return solution[2]; }

  using MultiStepSIM::saveStep;
  //! \brief Saves the converged solution to VTF file for a given time step.
  //! \param[in] iStep Time step identifier
  //! \param[in] param Time stepping parameters
  bool saveStep(int iStep, TimeStep& param);

  //! \brief Dumps solution variables at user-defined points.
  //! \param[in] time Current time
  //! \param[in] os The output stream to write the solution to
  //! \param[in] precision Number of digits after the decimal point
  //! \param[in] formatted If \e false, write all result points on a single line
  virtual void dumpResults(double time, utl::LogStream& os,
                           std::streamsize precision = 3,
                           bool formatted = true) const;

protected:
  // Time integration parameters
  double alpha1; //!< Mass-proportional damping parameter
  double alpha2; //!< Stiffness-proportional damping parameter
  double beta;   //!< Newmark time integration parameter
  double gamma;  //!< Newmark time integration parameter

  // Solution algorithm parameters
  bool   solveDisp; //!< If \e true, use incremental displacements as unknowns
  char   predictor; //!< Predictor type flag
  int    maxit;     //!< Maximum number of iterations in a time step
  int    maxIncr;   //!< Maximum number of iterations with increasing norm
  int    nupdat;    //!< Number of iterations with updated tangent
  int    saveIts;   //!< Time step for which iteration result should be saved
  bool   saveVelAc; //!< If \e true, save velocity and acceleration to VTF
  double rTol;      //!< Relative convergence tolerance
  double aTol;      //!< Absolute convergence tolerance
  double divgLim;   //!< Relative divergence limit

  unsigned short int cNorm; //!< Option for which convergence norm to use

  std::string factor_file; //!< File for storage of factorized Newton matrix
  bool write_factor;       //!< If \e true, write new factorized Newton matrix
  bool read_factor;        //!< If \e true, read factorized Newton matrix

public:
  static const char* inputContext; //!< Input file context for solver parameters
};

#endif
