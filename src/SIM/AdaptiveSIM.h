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
class SIMoutput;
class TimeStep;
class NormBase;


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
  AdaptiveSIM(SIMbase* sim = NULL);
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
  //! \param[in] nNormProj Number of element norms per projection method
  bool writeGlv(const char* infile, int iStep, size_t nNormProj);

  //! \brief Prints out the global norms to the log stream.
  void printNorms(size_t w = 36) const;

  //! \brief Accesses the solution of the linear system.
  const Vector& getSolution(size_t idx = 0) const { return solution[idx]; }

  //! \brief Accesses the projections.
  const Vector& getProjection(size_t idx) const { return projs[idx]; }

  //! \brief Accesses the norm prefices.
  const char** getNormPrefixes() { return &prefix[0]; }

  //! \brief Initializes the projections.
  void setupProjections();

  //! \brief Get number of norms in adaptor group
  int getNoNorms() const;

  //! \brief Dummy time stepping advance (no adaptive + time stepping yet)
  bool advanceStep(TimeStep& tp) const { return false; }

protected:
  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);
  //! \brief Parses a data section from an XML document.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem);

private:
  SIMoutput* model; //!< The isogeometric FE model

  bool   storeMesh;    //!< Creates a series of eps-files for intermediate steps
  bool   linIndepTest; //!< Test mesh for linear independence after refinement
  double beta;         //!< Refinement percentage in each step
  double errTol;       //!< Global error stop tolerance
  int    maxStep;      //!< Maximum number of adaptive refinements
  int    maxDOFs;      //!< Maximum number of degrees of freedom
  int    symmetry;     //!< Always refine a multiplum of this
  int    knot_mult;    //!< Knotline multiplicity
  int    maxTjoints;   //!< Maximum number of hanging nodes on one element
  double maxAspRatio;  //!< Maximum element aspect ratio
  bool   closeGaps;    //!< Split elements with a hanging node on each side
  bool   trueBeta;     //!< Beta measured in solution space dimension increase
  bool   threashold;   //!< Beta generates a threashold error/max(error)=beta

  //! Refinement scheme: 0=fullspan, 1=minspan, 2=isotropic_elements,
  //! 3=isotropic_functions
  int    scheme;

  size_t  adaptor;  //!< Norm group to base the mesh adaption on
  size_t  adNorm;   //!< Which norm to adapt based on
  Vectors solution; //!< All solutions (galerkin projections)
  Vectors gNorm;    //!< Global norms
  Matrix  eNorm;    //!< Element norms

  int geoBlk; //!< Running VTF geometry block counter
  int nBlock; //!< Running VTF result block counter

  std::vector<Vector>      projs;  //!< Projected secondary solutions
  std::vector<const char*> prefix; //!< Norm prefices for VTF-output
};

#endif
