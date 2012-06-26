// $Id$
//==============================================================================
//!
//! \file SIMLinElBeamC1.h
//!
//! \date Sep 16 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for NURBS-based FE analysis of C1-continuous beams.
//!
//==============================================================================

#ifndef _SIM_LIN_EL_BEAM_C1_H
#define _SIM_LIN_EL_BEAM_C1_H

#include "SIM1D.h"

class Material;


/*!
  \brief Driver class for isogeometric FEM analysis of C1-continuous beams.
*/

class SIMLinElBeamC1 : public SIM1D
{
public:

  /*!
    \brief Struct defining a nodal point load.
  */
  struct PointLoad
  {
    size_t patch; //!< Patch index [0,nPatch>
    int    inod;  //!< Local node number of the closest node
    double xi;    //!< Parameter of the point
    Vec3   X;     //!< Spatial coordinates of the point
    double pload; //!< Load magnitude
    // \brief Default constructor.
    PointLoad() : patch(0), inod(0) { xi = pload = 0.0; }
  };

  typedef std::vector<PointLoad> PloadVec; //!< Point load container

  //! \brief Default constructor.
  SIMLinElBeamC1();
  //! \brief Empty destructor.
  virtual ~SIMLinElBeamC1() {}

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \param[in] ignored Indices of patches to ignore in the analysis
  //! \param[in] fixDup Merge duplicated FE nodes on patch interfaces?
  virtual bool preprocess(const std::vector<int>& ignored = std::vector<int>(),
			  bool fixDup = false);

protected:
  //! \brief Parses a data section from the input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);
  //! \brief Parses a data section from an XML element.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Initializes material properties for integration of interior terms.
  //! \param[in] propInd Physical property index
  virtual bool initMaterial(size_t propInd);
  //! \brief Initializes the body load properties for current patch.
  //! \param[in] patchInd 1-based patch index
  virtual bool initBodyLoad(size_t patchInd);

  //! \brief Assembles the nodal point loads, if any.
  virtual bool assembleDiscreteTerms(const IntegrandBase* problem);

  //! \brief Computes problem-dependet external energy contributions.
  virtual double externalEnergy(const Vectors& psol) const;

private:
  std::vector<Material*> mVec;    //!< Material data
  RealArray              tVec;    //!< Beam thickness data
  PloadVec               myLoads; //!< Nodal point loads
};

#endif
