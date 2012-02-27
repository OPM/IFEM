// $Id$
//==============================================================================
//!
//! \file AdaptiveSIM.h
//!
//! \date Sep 22 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Adaptive solution driver for linear isogeometric FEM simulators.
//!
//==============================================================================

#ifndef _ADAPTIVE_SIM_H
#define _ADAPTIVE_SIM_H

#include "SIMinput.h"
#include "SystemMatrix.h"

class SIMbase;


/*!
  \brief Adaptive solution driver for linear isogeometric FEM simulators.
  \details This class contains data and methods for solving linear FE problems
  adaptively based on element error norms as refinement indicators.
*/

class AdaptiveSIM : public SIMinput
{
public:
  //! \brief The constructor initializes default adaptation parameters.
  //! \param sim Pointer to the spline FE model
  AdaptiveSIM(SIMbase* sim = 0);
  //! \brief The destructor frees the dynamically allocated FE model object.
  virtual ~AdaptiveSIM();

  //! \brief Assembles and solves the linear FE equations on current mesh.
  //! \param[in] inputfile File to read model parameters from after refinement
  //! \param[in] solver The linear equation solver to use
  //! \param[in] iStep Refinement step counter
  bool solveStep(const char* inputfile, SystemMatrix::Type solver, int iStep);

  //! \brief Refines the current mesh based on the element norms.
  //! \param[in] iStep Refinement step counter
  bool adaptMesh(int iStep);

  //! \brief Writes current mesh and results to the VTF-file.
  //! \param[in] infile File name used to construct the VTF-file name from
  //! \param[in] format Format of VTF-file (0=ASCII, 1=BINARY)
  //! \param[in] nViz   Number of visualization points over a knot-span
  //! \param[in] iStep  Refinement step identifier
  //! \param     nBlock Running result block counter
  bool writeGlv(const char* infile, int format, const int* nViz,
		int iStep, int& nBlock);

  //! \brief Prints out the global norms to given stream.
  static std::ostream& printNorms(const Vector& norms, const Matrix& eNorm,
				  std::ostream& os);

protected:

  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);

  //! \brief Parses a data section from an XML document.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem);

private:
  SIMbase* model; //!< The isogeometric FE model

  bool   storeMesh; //!< Creates a series of .eps files for intermediate steps
  double beta;      //!< Refinement percentage in each step
  double errTol;    //!< Global error stop tolerance
  int    maxStep;   //!< Maximum number of adaptive refinements
  int    maxDOFs;   //!< Maximum number of degrees of freedom
  int    symmetry;  //!< Always refine a multiplum of this
  int    scheme;    //!< Refinement scheme: 0=fullspan, 1=minspan, 2=isotropic_elements, 3=isotropic_functions
  int    knot_mult; //!< Knotline multiplicity

  std::vector<int> options; //!< Mesh refinement options

  Vector linsol; //!< Linear solution vector
  Vector gNorm;  //!< Global norms
  Matrix eNorm;  //!< Element norms
};

#endif
