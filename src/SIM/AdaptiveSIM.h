// $Id$
//==============================================================================
//!
//! \file AdaptiveSIM.h
//!
//! \date Sep 22 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Adaptive solution driver for linear static FEM simulators.
//!
//==============================================================================

#ifndef _ADAPTIVE_SIM_H
#define _ADAPTIVE_SIM_H

#include "SIMadmin.h"
#include "AdaptiveSetup.h"


/*!
  \brief Adaptive solution driver for linear static FEM simulators.
  \details This class contains data and methods for solving linear static FE
  problems adaptively, based on element error norms as refinement indicators.
*/

class AdaptiveSIM : public SIMadmin, public AdaptiveSetup
{
public:
  //! \brief The constructor initializes default adaptation parameters.
  //! \param sim The FE model
  //! \param[in] sa If \e true, this is a stand-alone driver
  explicit AdaptiveSIM(SIMoutput& sim, bool sa = true);
  //! \brief Empty destructor.
  virtual ~AdaptiveSIM() {}

  //! \brief Initializes the \a projs and \a prefix arrays.
  //! \param[in] normGroup Index to the norm group to base mesh adaptation on
  bool initAdaptor(size_t normGroup = 0);

  //! \brief Assembles and solves the linear FE equations on current mesh.
  //! \param[in] inputfile File to read model parameters from after refinement
  //! \param[in] iStep Refinement step counter
  //! \param[in] withRF Whether nodal reaction forces should be computed or not
  //! \param[in] precision Number of digits after decimal point
  bool solveStep(const char* inputfile, int iStep, bool withRF = false,
                 std::streamsize precision = 6);

  //! \brief Refines the current mesh based on the element norms.
  //! \param[in] iStep Refinement step counter
  //! \param[in] outPrec Number of digits after the decimal point in norm print
  bool adaptMesh(int iStep, std::streamsize outPrec = 0);

  //! \brief Writes current mesh and results to the VTF-file.
  //! \param[in] infile File name used to construct the VTF-file name from
  //! \param[in] iStep  Refinement step identifier
  bool writeGlv(const char* infile, int iStep);

  //! \brief Accesses the solution of the linear system.
  const Vector& getSolution(size_t idx = 0) const { return solution[idx]; }
  //! \brief Accesses the projections.
  const Vector& getProjection(size_t idx = 0) const { return projs[idx]; }
  //! \brief Access all the projections.
  const Vectors& getProjections() const { return projs; }
  //! \brief Access the calculated element-wise norms.
  const Matrix& getEnorm() const { return eNorm; }

  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);
  //! \brief Parses a data section from an XML document.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Saves point results to output file for a given time step.
  //! \param[in] time Load/time step parameter
  //! \param[in] step Load/time step counter
  //! \details By default it just forwards to the underlying model
  virtual bool savePoints(double time, int iStep) const;

protected:
  //! \brief Assembles and solves the linear FE equation system.
  virtual bool assembleAndSolveSystem();

private:
  Vectors gNorm; //!< Global norms
  Vectors dNorm; //!< Dual global norms
  Matrix  eNorm; //!< Element norms
  Matrix  fNorm; //!< Dual element norms

  int geoBlk; //!< Running VTF geometry block counter
  int nBlock; //!< Running VTF result block counter

  std::vector<Vector>      projs;  //!< Projected secondary solutions
  std::vector<Vector>      projd;  //!< Projected dual solutions
  std::vector<std::string> prefix; //!< Norm prefices for VTF-output

protected:
  Vectors solution; //!< All solutions (including Galerkin projections)
};

#endif
