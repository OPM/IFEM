// $Id$
//==============================================================================
//!
//! \file AdaptiveSIM.h
//!
//! \date Sep 22 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Adaptive solution driver for isogeometric FEM simulators.
//!
//==============================================================================

#ifndef _ADAPTIVE_SIM_H
#define _ADAPTIVE_SIM_H

#include "SIMinput.h"
#include "SystemMatrix.h"

class SIMbase;


/*!
  \brief Nonlinear solution driver for isogeometric FEM simulators.
  \details This class contains data and methods for computing the nonlinear
  solution to a FE problem based on splines/NURBS basis functions, through
  Newton-Raphson iterations.
*/

class AdaptiveSIM : public SIMinput
{
public:
  //! \brief The constructor initialized default solution parameters.
  //! \param sim Pointer to the spline FE model
  AdaptiveSIM(SIMbase* sim = 0);
  //! \brief The destructor frees the dynamically allocated FE model object.
  virtual ~AdaptiveSIM();

  //! \brief Solves the nonlinear equations by Newton-Raphson iterations.
  //! \param param Solution algorithm parameters
  //! \param[in] mode Solution mode to use for this step
  //! \param[in] compName Solution name to be used in the norm output
  //! \param[in] energyNorm If \e true, integrate energy norm of the solution
  //! \param[in] zero_tolerance Truncate norm values small than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  bool solveStep(const char* inputfile, SystemMatrix::Type solver, int iStep);

  //! \brief Computes and prints some solution norm quantities.
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param[in] compName Solution name to be used in the norm output
  //! \param[in] energyNorm If \e true, integrate energy norm of the solution
  //! \param[in] zero_tolerance Truncate norm values small than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  bool adaptMesh(int iStep);

  //! \brief Prints out the global norms to given stream
  static void printNorms(const Vector& norms, std::ostream& os);

protected:

  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);

private:
  SIMbase* model; //!< The isogeometric FE model

  int    nStep;
  double stopTol;
  double beta;

  Vector linsol; //!< Linear solution vector
  Vector gNorm;  //!< Global norms
  Matrix eNorm;  //!< Element norms
};

#endif
