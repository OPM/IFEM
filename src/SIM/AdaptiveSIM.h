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
#include "MatVec.h"

class SIMoutput;


/*!
  \brief Adaptive solution driver for linear static FEM simulators.
  \details This class contains data and methods for solving linear static FE
  problems adaptively, based on element error norms as refinement indicators.
*/

class AdaptiveSIM : public SIMadmin
{
public:
  //! \brief The constructor initializes default adaptation parameters.
  //! \param sim The FE model
  //! \param sa If \e true, this is a stand-alone driver
  AdaptiveSIM(SIMoutput& sim, bool sa = true);
  //! \brief The destructor frees the norm prefices.
  virtual ~AdaptiveSIM();

  //! \brief Sets the norm group/index of the norm to base mesh adaptation on.
  void setAdaptationNorm(size_t normGroup, size_t normIdx = 0);
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

  //! \brief Prints out the global norms to the log stream.
  void printNorms(size_t w = 36) const;

  //! \brief Accesses the solution of the linear system.
  const Vector& getSolution(size_t idx = 0) const { return solution[idx]; }
  //! \brief Accesses the projections.
  const Vector& getProjection(size_t idx = 0) const { return projs[idx]; }

  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);
  //! \brief Parses a data section from an XML document.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem);

private:
  SIMoutput& model; //!< The isogeometric FE model
  bool       alone; //!< If \e false, this class is wrapped by SIMSolver

  bool   storeMesh;    //!< Creates a series of eps-files for intermediate steps
  bool   linIndepTest; //!< Test mesh for linear independence after refinement
  double beta;         //!< Refinement percentage in each step
  double errTol;       //!< Global error stop tolerance
  double condLimit;    //!< Upper limit on condition number
  double rCond;        //!< Actual reciprocal condition number
  int    maxStep;      //!< Maximum number of adaptive refinements
  int    maxDOFs;      //!< Maximum number of degrees of freedom
  int    symmetry;     //!< Always refine a multiplum of this
  int    knot_mult;    //!< Knotline multiplicity
  int    maxTjoints;   //!< Maximum number of hanging nodes on one element
  double maxAspRatio;  //!< Maximum element aspect ratio
  bool   closeGaps;    //!< Split elements with a hanging node on each side

  //! Threshold flag for how to interpret the refinement percentage, \a beta
  enum { NONE, MAXIMUM, AVERAGE, MINIMUM, TRUE_BETA, DORFEL } threshold;

  //! Refinement scheme: 0=fullspan, 1=minspan, 2=isotropic_elements,
  //! 3=isotropic_functions
  int scheme;

  size_t  adaptor;  //!< Norm group to base the mesh adaptation on
  size_t  adNorm;   //!< Which norm to base the mesh adaptation on
  Vectors solution; //!< All solutions (galerkin projections)
  Vectors gNorm;    //!< Global norms
  Matrix  eNorm;    //!< Element norms

  int geoBlk; //!< Running VTF geometry block counter
  int nBlock; //!< Running VTF result block counter

  std::vector<Vector>      projs;  //!< Projected secondary solutions
  std::vector<const char*> prefix; //!< Norm prefices for VTF-output
};

#endif
