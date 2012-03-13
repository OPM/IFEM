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
#include "MatVec.h"

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

  //! \brief Initializes index to element norms to base the mesh adaption on.
  //! \param[in] indxProj Index to projection method to base mesh adaption on
  //! \param[in] nNormProj Number of element norms per projection method
  bool initAdaptor(size_t indxProj, size_t nNormProj);

  //! \brief Assembles and solves the linear FE equations on current mesh.
  //! \param[in] inputfile File to read model parameters from after refinement
  //! \param[in] iStep Refinement step counter
  bool solveStep(const char* inputfile, int iStep);

  //! \brief Refines the current mesh based on the element norms.
  //! \param[in] iStep Refinement step counter
  bool adaptMesh(int iStep);

  //! \brief Writes current mesh and results to the VTF-file.
  //! \param[in] infile File name used to construct the VTF-file name from
  //! \param[in] iStep  Refinement step identifier
  //! \param     nBlock Running result block counter
  bool writeGlv(const char* infile, int iStep, int& nBlock);

  //! \brief Prints out the global norms to given stream.
  static std::ostream& printNorms(const Vector& norms, const Matrix& eNorm,
				  std::ostream& os, size_t adaptor = 4);

  //! \brief Accesses the solution of the linear system.
  Vector& getSolution() { return linsol; }

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

  bool   storeMesh; //!< Creates a series of eps-files for intermediate steps
  double beta;      //!< Refinement percentage in each step
  double errTol;    //!< Global error stop tolerance
  int    maxStep;   //!< Maximum number of adaptive refinements
  int    maxDOFs;   //!< Maximum number of degrees of freedom
  int    symmetry;  //!< Always refine a multiplum of this
  int    knot_mult; //!< Knotline multiplicity

  //! Refinement scheme: 0=fullspan, 1=minspan, 2=isotropic_elements,
  //! 3=isotropic_functions
  int    scheme;

  size_t adaptor; //!< Index to the element norm to base mesh adaption on
  Vector linsol;  //!< Linear solution vector
  Vector gNorm;   //!< Global norms
  Matrix eNorm;   //!< Element norms
};

#endif
