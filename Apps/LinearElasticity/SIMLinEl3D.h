// $Id$
//==============================================================================
//!
//! \file SIMLinEl3D.h
//!
//! \date Dec 08 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for 3D NURBS-based linear elastic FEM analysis.
//!
//==============================================================================

#ifndef _SIM_LIN_EL_3D_H
#define _SIM_LIN_EL_3D_H

#include "SIM3D.h"
#include "SIMenums.h"

class Elasticity;
class Material;


/*!
  \brief Driver class for 3D isogeometric FEM analysis of elasticity problems.
  \details The class incapsulates data and methods for solving linear elasticity
  problems using NURBS-based finite elements. It reimplements the parse method
  of the parent class.
*/

class SIMLinEl3D : public SIM3D
{
public:
  //! \brief Default constructor.
  //! \param[in] checkRHS If \e true, ensure the model is in a right-hand system
  SIMLinEl3D(bool checkRHS = false) : SIM3D(checkRHS) { aCode = 0; }
  //! \brief The destructor frees the dynamically allocated material properties.
  virtual ~SIMLinEl3D();

  //! \brief Initializes the property containers of the model.
  virtual void clearProperties();

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented inserting a call to \a getIntegrand.
  //! This makes sure the integrand has been allocated in case of minimum input.
  //! It also resolves inhomogeneous boundary condition fields in case they are
  //! derived from the analytical solution.
  virtual bool preprocess(const std::vector<int>& ignored, bool fixDup);

  //! \brief Print norms to stream
  std::ostream& printNorms(const Vectors& norms, std::ostream& os);
private:
  //! \brief Returns the actual integrand.
  Elasticity* getIntegrand();

protected:
  //! \brief Parses a data section from the input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);

  //! \brief Parses a data section from an XML element
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Initializes material properties for integration of interior terms.
  //! \param[in] propInd Physical property index
  virtual bool initMaterial(size_t propInd);
  //! \brief Initializes the body load properties for current patch.
  //! \param[in] patchInd 1-based patch index
  virtual bool initBodyLoad(size_t patchInd);

  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  virtual bool initNeumann(size_t propInd);

protected:
  std::vector<Material*> mVec; //!< Material data

private:
  int aCode; //!< Analytical BC code (used by destructor)
};

#endif
