// $Id$
//==============================================================================
//!
//! \file AdaptiveSetup.h
//!
//! \date May 7 2019
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Adaptive solution setup for linear and nonlinear FEM simulators.
//!
//==============================================================================

#ifndef _ADAPTIVE_SETUP_H
#define _ADAPTIVE_SETUP_H

#include "MatVec.h"

class SIMoutput;
class TiXmlElement;
namespace LR { struct RefineData; }


/*!
  \brief Adaptive solution setup for linear and nonlinear FEM simulators.
  \details This class contains parameters for controlling the mesh refinement
  in an adaptive simulation procedure, with methods for parsing the parameters
  from the input file and to calculate the actual refinement indicators to be
  passed to the LR-spline based mesh refinement methods. It also has methods
  for printing out the global error norms that are used in the mesh-refinement
  calculation, and for dumping the refined meshes to files in different formats.
*/

class AdaptiveSetup
{
public:
  //! \brief The constructor initializes default adaptation parameters.
  //! \param sim The FE model
  //! \param[in] sa If \e true, this is a stand-alone driver
  explicit AdaptiveSetup(SIMoutput& sim, bool sa = true);
  //! \brief Empty destructor.
  virtual ~AdaptiveSetup() {}

  //! \brief Sets the norm group/index of the norm to base mesh adaptation on.
  void setAdaptationNorm(size_t g, size_t i = 0) { adaptor = g; adNorm = i; }

  //! \brief Initializes the norm group parameters.
  //! \param[in] normGroup Index to the norm group to base mesh adaptation on
  bool initPrm(size_t normGroup);

  //! \brief Calculates mesh refinement control data based on error estimates.
  //! \param[out] prm Mesh refinement control data
  //! \param[in] iStep Refinement step counter
  //! \param[in] gNorm Global norms
  //! \param[in] refIn Element refinement indicators (element error norms)
  //! \return Number of elements to be refined
  //! \return If zero, no refinement needed
  //! \return Negative value on error
  int calcRefinement(LR::RefineData& prm, int iStep,
                     const Vectors& gNorm, const Vector& refIn) const;

  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  bool parse(char* keyWord, std::istream& is);
  //! \brief Parses a data section from an XML document.
  //! \param[in] elem The XML element to parse
  bool parse(const TiXmlElement* elem);

  //! \brief Prints out global norms to the log stream.
  //! \param[in] gNorm Global norms
  //! \param[in] dNorm Global dual norms
  //! \param[in] eNorm Element norms
  //! \param[in] w Field width for global norm labels
  //! \param[in] printModelNorms True to print norms for model
  void printNorms(const Vectors& gNorm, const Vectors& dNorm,
                  const Matrix& eNorm, size_t w = 36,
                  bool printModelNorms = true) const;

  //! \brief Dumps current mesh to external file(s) for inspection.
  //! \param[in] iStep Current refinement step (1=initial grid)
  bool writeMesh(int iStep) const;

  //! \brief Returns the row-index of the element norm to use for adaptation.
  size_t eIdx() const { return eRow; }

protected:
  SIMoutput& model; //!< The isogeometric FE model

  size_t adaptor; //!< Norm group to base the mesh adaptation on
  size_t adNorm;  //!< Which norm to base the mesh adaptation on
  size_t eRow;    //!< Row-index in \a eNorm of the norm to use for adaptation
  double rCond;   //!< Actual reciprocal condition number of the last mesh

private:
  bool   alone;      //!< If \e false, this class is wrapped by SIMSolver
  bool   linIndep;   //!< Test mesh for linear independence after refinement
  double beta;       //!< Refinement percentage in each step
  double errTol;     //!< Global error stop tolerance
  double condLimit;  //!< Upper limit on condition number
  int    maxStep;    //!< Maximum number of adaptive refinements
  int    maxDOFs;    //!< Maximum number of degrees of freedom
  int    knot_mult;  //!< Knotline multiplicity
  int    maxTjoints; //!< Maximum number of hanging nodes on one element
  double maxAspect;  //!< Maximum element aspect ratio
  bool   closeGaps;  //!< Split elements with a hanging node on each side
  double symmEps;    //!< Epsilon used for symmetrized selection method

  //! \brief Enum defining the refinement threshold flag values.
  enum Threshold { NONE=0, MAXIMUM, AVERAGE, MINIMUM, TRUE_BETA,
                   DORFEL, SYMMETRIZED };

  //! \brief Enum defining available refinement scheme options.
  enum RefScheme { FULLSPAN=0, MINSPAN=1,
                   ISOTROPIC_FUNCTION=2,
                   ISOTROPIC_ELEMENT=3 };

  Threshold   threshold; //!< Flag for how to interpret the parameter \a beta
  RefScheme   scheme;    //!< The actual refinement scheme to use

  int         storeMesh; //!< Flag telling what kind of mesh output we want
  std::string mshPrefix; //!< Prefix for output files with refined meshes
  std::string errPrefix; //!< Prefix for text files with refinement indicators
};

#endif
