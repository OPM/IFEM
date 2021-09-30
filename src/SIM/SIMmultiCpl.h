// $Id$
//==============================================================================
//!
//! \file SIMmultiCpl.h
//!
//! \date Feb 12 2021
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Monolithic coupling of multiple simulators.
//!
//==============================================================================

#ifndef _SIM_MULTI_CPL_H_
#define _SIM_MULTI_CPL_H_

#include "SIMadmin.h"
#include "TopologySet.h"

class SIMoutput;
class SIMinput;


/*!
  \brief Class for monolithic coupled simulators.
*/

class SIMmultiCpl : public SIMadmin
{
public:
  //! \brief The constructor initializes the array of base %SIM objects.
  explicit SIMmultiCpl(const std::vector<SIMoutput*>& sims);
  //! \brief The destructor deletes the base %SIM objects.
  virtual ~SIMmultiCpl();

  using SIMadmin::parse;
  //! \brief Parses a data section from an XML document.
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \param[in] ignored Indices of patches to ignore in the analysis
  //! \param[in] fixDup Merge duplicated FE nodes on patch interfaces?
  virtual bool preprocess(const std::vector<int>& ignored, bool fixDup);

private:
  //! \brief Parses a connection definition from an XML element.
  bool parseConnection(const TiXmlElement* elem);

protected:
  //! \brief Struct defining an inter-simulator coupling.
  struct SIMcoupling
  {
    SIMinput* mstSim = nullptr; //!< Simulator of the master point
    SIMinput* slvSim = nullptr; //!< Simulator of the slave point
    TopEntity::const_iterator master; //!< Master point specification
    TopEntity::const_iterator slave;  //!< Slave point specification
  };

  std::vector<SIMcoupling> myCpl;  //!< Inter-sim coupling definition container
  std::vector<SIMoutput*>  mySims; //!< Dimension-specific simulator container
};

#endif
