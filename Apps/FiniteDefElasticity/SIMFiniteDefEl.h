// $Id$
//==============================================================================
//!
//! \file SIMFiniteDefEl.h
//!
//! \date Dec 18 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution drivers for NURBS-based finite deformation analysis.
//!
//==============================================================================

#ifndef _SIM_FINITE_DEF_EL_H
#define _SIM_FINITE_DEF_EL_H

#include "SIMLinEl2D.h"
#include "SIMLinEl3D.h"
#include "SIMContact.h"
#include "NLoptions.h"


/*!
  \brief Driver class for 2D isogeometric finite deformation analysis.
*/

class SIMFiniteDefEl2D : public SIMLinEl2D, private SIMContact
{
public:
  //! \brief Default constructor.
  //! \param[in] options Additional solution formulation options
  SIMFiniteDefEl2D(const std::vector<int>& options = std::vector<int>());
  //! \brief The destructor deletes the additional integrands.
  virtual ~SIMFiniteDefEl2D();

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

  //! \brief Creates a property set for contact condition on an entity set.
  //! \param[in] slave Name of the slave entity set
  //! \param[out] code Property code associated with the contact set
  virtual bool createContactSet(const std::string& slaveSet, int& code);

public:
  //! \brief Updates the time-dependent in-homogeneous Dirichlet coefficients.
  //! \param[in] time Current time
  //! \param[in] prevSol Pointer to previous primary solution in DOF-order
  virtual bool updateDirichlet(double time, const Vector* prevSol);

  //! \brief Updates Problem-dependent state based on current solution.
  //! \param[in] solution Current primary solution vector
  virtual bool updateConfiguration(const Vector& solution);

private:
  NLoptions nlo; //!< Input options defining the nonlinear formulation
};


/*!
  \brief Driver class for 3D isogeometric finite deformation analysis.
*/

class SIMFiniteDefEl3D : public SIMLinEl3D, private SIMContact
{
public:
  //! \brief Default constructor.
  //! \param[in] checkRHS If \e true, ensure the model is in a right-hand system
  //! \param[in] options Additional solution formulation options
  SIMFiniteDefEl3D(bool checkRHS = false,
		   const std::vector<int>& options = std::vector<int>());
  //! \brief The destructor deletes the additional integrands.
  virtual ~SIMFiniteDefEl3D();

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

  //! \brief Creates a property set for contact condition on an entity set.
  //! \param[in] slave Name of the slave entity set
  //! \param[out] code Property code associated with the contact set
  virtual bool createContactSet(const std::string& slaveSet, int& code);

public:
  //! \brief Updates the time-dependent in-homogeneous Dirichlet coefficients.
  //! \param[in] time Current time
  //! \param[in] prevSol Pointer to previous primary solution in DOF-order
  virtual bool updateDirichlet(double time, const Vector* prevSol);

  //! \brief Updates Problem-dependent state based on current solution.
  //! \param[in] solution Current primary solution vector
  virtual bool updateConfiguration(const Vector& solution);

private:
  NLoptions nlo; //!< Input options defining the nonlinear formulation
};

#endif
