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

#include "SIMinput.h"
#include "SIMenums.h"
#include "MatVec.h"

class SIMbase;
class TimeStep;


/*!
  \brief Newmark-based solution driver for dynamic isogeometric FEM simulators.
*/

class NewmarkSIM : public SIMinput
{
public:
  //! \brief The constructor initializes default solution parameters.
  //! \param sim The spline FE model
  NewmarkSIM(SIMbase& sim);
  //! \brief Empty destructor.
  virtual ~NewmarkSIM() {}

  //! \brief Parses a data section from an input stream (depreciated).
  virtual bool parse(char*, std::istream&) { return false; }
  //! \brief Parses a data section from an XML document.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Initializes primary solution vectors and integration parameters.
  void init();

  //! \brief Advances the time step one step forward.
  //! \param param Time stepping parameters
  //! \param[in] updateTime If \e false, the time parameters are not incremented
  virtual bool advanceStep(TimeStep& param, bool updateTime = true);

  //! \brief Solves the dynamic equations by a predictor/multi-corrector method.
  //! \param param Time stepping parameters
  //! \param[in] zero_tolerance Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  SIM::ConvStatus solveStep(TimeStep& param, double zero_tolerance = 1.0e-8,
                            std::streamsize outPrec = 0);

  //! \brief Computes and prints some solution norm quantities.
  //! \param[in] zero_tolerance Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  virtual bool solutionNorms(double zero_tolerance = 1.0e-8,
                             std::streamsize outPrec = 0);

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  bool saveModel(char* fileName);

  //! \brief Saves the converged results to VTF file of a given time step.
  //! \param[in] iStep Time step identifier
  //! \param[in] time Current time parameter
  //! \param[in] psolOnly If \e true, skip secondary solution field output
  //! \param[in] vecName Optional name of primary solution vector field
  bool saveStep(int iStep, double time,
                bool psolOnly = false, const char* vecName = NULL);

  //! \brief Returns a const reference to current solution vector.
  const Vector& getSolution(int idx = 0) const { return solution[idx]; }

protected:
  //! \brief Checks whether the corrector iterations have converged or diverged.
  SIM::ConvStatus checkConvergence(TimeStep& param);
  //! \brief Calculates predicted velocities and accelerations.
  virtual bool predictStep(TimeStep& param);
  //! \brief Updates configuration variables (solution vector) in an iteration.
  virtual bool correctStep(TimeStep& param);

protected:
  SIMbase& model;    //!< The isogeometric FE model
  Vectors  solution; //!< Primary solution vectors of the last increments
  Vector   linsol;   //!< Linear solution vector
  Vector   residual; //!< Residual force vector

  // Time integration parameters
  double alpha1; //!< Mass-proportional damping parameter
  double alpha2; //!< Stiffness-proportional damping parameter
  double beta;   //!< Newmark time integration parameter
  double gamma;  //!< Newmark time integration parameter

  // Solution algorithm parameters
  char   predictor; //!< Predictor type flag
  int    maxit;     //!< Maximum number of iterations in a time step
  double convTol;   //!< Convergence tolerance
  double divgLim;   //!< Relative divergence limit
  double refNorm;   //!< Reference norm value used in convergence checks

  // Post-processing attributes
  int geoBlk; //!< Running VTF geometry block counter
  int nBlock; //!< Running VTF result block counter
};

#endif
