// $Id$
//==============================================================================
//!
//! \file NonlinearDriver.h
//!
//! \date Jul 15 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Nonlinear driver for isogeometric finite deformation FEM analysis.
//!
//==============================================================================

#ifndef _NONLINEAR_DRIVER_H
#define _NONLINEAR_DRIVER_H

#include "NonLinSIM.h"
#include "TimeStep.h"

class DataExporter;


/*!
  \brief Nonlinear driver for isogeometric finite deformation FEM analysis.
  \details This class is derived from NonLinSIM of the IFEM kernel.
  It reimplements the \a solutionNorms method to also compute the energy norm
  and other norms of the stress field. In addition, it has the method
  \a solveProblem to manage the pseudo-time step loop.
*/

class NonlinearDriver : public NonLinSIM
{
public:
  //! \brief Default constructor.
  //! \param sim Pointer to the spline FE model
  NonlinearDriver(SIMbase* sim = 0) : NonLinSIM(sim) {}
  //! \brief Empty destructor.
  virtual ~NonlinearDriver() {}

protected:
  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);

  //! \brief Parses a data section from an XML document.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Computes and prints some solution norm quantities.
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param[in] energyNorm If \e true, integrate energy norm of the solution
  //! \param[in] zero_tol Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  virtual bool solutionNorms(const TimeDomain& time, bool energyNorm,
                             double zero_tol, std::streamsize outPrec);

public:
  //! \brief Invokes the main pseudo-time stepping simulation loop.
  //! \param[in] skip2nd If \e true, output only primary solution variables
  //! \param[in] energyNorm If \e true, integrate energy norm of the solution
  //! \param writer HDF5 results exporter
  //! \param oss Output stream for additional ASCII result output
  //! \param[in] dtDump Time increment for dumo of ASCII results
  //! \param[in] zero_tol Truncate norm values smaller than this to zero
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  int solveProblem(bool skip2nd, bool energyNorm,
                   DataExporter* writer, std::ostream* oss, double dtDump,
                   double zero_tol, std::streamsize outPrec);

private:
  TimeStep params; //!< Time stepping parameters
};

#endif
