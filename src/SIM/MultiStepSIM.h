// $Id$
//==============================================================================
//!
//! \file MultiStepSIM.h
//!
//! \date Jul 11 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for multi-step solution drivers.
//!
//==============================================================================

#ifndef _MULTI_STEP_SIM_H
#define _MULTI_STEP_SIM_H

#include "SIMinput.h"
#include "SIMenums.h"
#include "MatVec.h"

class SIMbase;
class SIMoutput;
class TimeStep;
struct TimeDomain;


/*!
  \brief Base class for multi-step solution drivers.
*/

class MultiStepSIM : public SIMinput
{
protected:
  //! \brief Enum describing reference norm options.
  enum NormOp { MAX, ALL };

  //! \brief The constructor initializes the FE model reference.
  //! \param sim The FE model
  MultiStepSIM(SIMbase& sim);

public:
  //! \brief Empty destructor.
  virtual ~MultiStepSIM() {}

  //! \brief Prints out problem-specific data to the log stream.
  virtual void printProblem() const;

  //! \brief Returns a list of prioritized XML-tags.
  virtual const char** getPrioritizedTags() const;

  //! \brief Initializes time integration parameters for the integrand.
  virtual void initPrm() {}
  //! \brief Initializes the primary solution vectors.
  //! \param[in] nSol Number of consequtive solutions stored in core
  virtual bool initSol(size_t nSol = 1);

  //! \brief Allocates the FE system matrices.
  //! \param[in] withRF Whether nodal reaction forces should be computed or not
  bool initEqSystem(bool withRF = true);

  //! \brief Advances the time/load step one step forward.
  //! \param param Time stepping parameters
  //! \param[in] updateTime If \e false, the time parameters are not incremented
  virtual bool advanceStep(TimeStep& param, bool updateTime = true);

  //! \brief Solves the FE equations at current time/load step.
  //! \param param Time stepping parameters
  //! \param[in] zero_tolerance Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  virtual SIM::ConvStatus solveStep(TimeStep& param,
                                    SIM::SolutionMode = SIM::STATIC,
                                    double zero_tolerance = 1.0e-8,
                                    std::streamsize outPrec = 0) = 0;

protected:
  //! \brief Computes and prints some solution norm quantities.
  //! \param[in] zero_tolerance Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  virtual bool solutionNorms(const TimeDomain&, double zero_tolerance = 1.0e-8,
                             std::streamsize outPrec = 0);

  //! \brief Returns the last step that was save to VTF
  int getLastSavedStep() const { return lastSt; }

public:
  //! \brief Initializes the geometry block counter.
  void setStartGeo(int gID);

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  bool saveModel(char* fileName);

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param gBlock Running geometry block counter
  //! \param rBlock Running result block counter
  //! \param[in] fileName File name used to construct the VTF-file name from
  bool saveModel(int& gBlock, int& rBlock, char* fileName = nullptr);

  //! \brief Saves the converged results to VTF file of a given time/load step.
  //! \param[in] iStep Time/load step identifier
  //! \param[in] time Current time/load parameter
  //! \param[in] vecName Optional name of primary solution vector field
  bool saveStep(int iStep, double time, const char* vecName = nullptr);

  //! \brief Saves the converged results to VTF file of a given time/load step.
  //! \param[in] iStep Time/load step identifier
  //! \param rBlock Running result block counter
  //! \param[in] time Current time/load parameter
  //! \param[in] vecName Optional name of primary solution vector field
  bool saveStep(int iStep, int& rBlock, double time,
                const char* vecName = nullptr);

  //! \brief Saves the converged solution to VTF file of a given time/load step.
  //! \param[in] iStep Time/load step identifier
  //! \param rBlock Running result block counter
  //! \param[in] vecName Name of primary solution vector field
  bool saveStep(int iStep, int& rBlock, const char* vecName);

  //! \brief Dumps the primary solution for inspection.
  //! \param[in] iStep Time/load step identifier
  //! \param[in] time Current time/load parameter
  //! \param[in] os The output stream to write the solution to
  //! \param[in] withID If \e true, write node ID and coordinates too
  void dumpStep(int iStep, double time, utl::LogStream& os,
                bool withID = true) const;

  //! \brief Dumps solution variables at user-defined points.
  //! \param[in] time Current time/load parameter
  //! \param[in] os The output stream to write the solution to
  //! \param[in] precision Number of digits after the decimal point
  //! \param[in] formatted If \e false, write all result points on a single line
  virtual void dumpResults(double time, utl::LogStream& os,
                           std::streamsize precision = 3,
                           bool formatted = true) const;

  //! \brief Returns whether a points result file has been defined or not.
  bool hasPointResultFile() const;
  //! \brief Saves point-wise solution to file for a given time/load step.
  //! \param[in] time Time/load step parameter
  //! \param[in] step Time/load step counter
  bool savePoints(double time, int step) const;

  //! \brief Returns a const reference to the solution vectors.
  const Vectors& getSolutions() const { return solution; }
  //! \brief Returns a const reference to current solution vector.
  const Vector& getSolution(int idx = 0) const { return solution[idx]; }
  //! \brief Modifies the current solution vector (used by sub-iterations only).
  virtual void setSolution(const Vector& s, int idx = 0) { solution[idx] = s; }

  //! \brief Returns a const reference to the FEmodel.
  const SIMoutput& getModel() const { return model; }

  //! \brief Enum describing sub-iteration status.
  enum SubIt { ITER = 0, FIRST = 1, LAST = 2, NONE = 3 };

  //! \brief Updates the sub-iteration flag.
  void setSubIteration(SubIt flag) { subiter = flag; }
  //! \brief Returns the sub-iteration flag.
  SubIt getSubIteration() const { return subiter; }

protected:
  SIMoutput& model;    //!< The isogeometric FE model
  Vector     residual; //!< Residual force vector
  Vector     linsol;   //!< Linear solution vector
  Vectors    solution; //!< Primary solution vectors

  NormOp refNopt; //!< Reference norm option
  double refNorm; //!< Reference norm value used in convergence checks
  SubIt  subiter; //!< Subiteration flag
  size_t nRHSvec; //!< Number of right-hand-side vectors to assemble
  char   rotUpd;  //!< Option for how to update of nodal rotations

  int geoBlk; //!< Running VTF geometry block counter
  int nBlock; //!< Running VTF result block counter

private:
  int lastSt; //!< The last step that was saved to VTF
};

#endif
