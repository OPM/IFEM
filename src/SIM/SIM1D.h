// $Id$
//==============================================================================
//!
//! \file SIM1D.h
//!
//! \date Feb 04 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for 1D NURBS-based FEM analysis.
//!
//==============================================================================

#ifndef _SIM_1D_H
#define _SIM_1D_H

#include "SIMgeneric.h"


/*!
  \brief Driver class for 1D NURBS-based FEM solver.
*/

class SIM1D : public SIMgeneric
{
public:
  //! \brief Enum announcing the dimensionality (used for template writing).
  enum { dimension = 1 };

  //! \brief Default constructor.
  //! \param[in] n1 Dimension of the primary solution field
  explicit SIM1D(unsigned char n1 = 1, bool = false);
  //! \brief Constructor used for mixed problems.
  //! \param[in] unf Dimension of the primary solution fields
  explicit SIM1D(const CharVec& unf, bool = false);
  //! \brief Constructor that also initializes the integrand pointer.
  //! \param[in] itg Pointer to the integrand of the problem to solve
  //! \param[in] n Dimension of the primary solution field
  explicit SIM1D(IntegrandBase* itg, unsigned char n = 1);
  //! \brief Empty destructor.
  virtual ~SIM1D() {}

  //! \brief Returns the number of parameter dimensions in the model.
  virtual unsigned short int getNoParamDim() const { return 1; }

  //! \brief Evaluates the primary solution at the given point.
  //! \param[in] psol Primary solution vector
  //! \param[in] u Parameter of the point to evaluate at
  //! \param[in] deriv Derivative order of the solution
  //! \param[in] patch 1-based patch index contining the evaluation point
  //! \return Evaluated solution values
  Vector getSolution(const Vector& psol, double u,
                     int deriv = 0, int patch = 1) const;

private:
  //! \brief Parses a subelement of the \a geometry XML-tag.
  bool parseGeometryTag(const tinyxml2::XMLElement* elem);
  //! \brief Parses a subelement of the \a boundaryconditions XML-tag.
  bool parseBCTag(const tinyxml2::XMLElement* elem);

protected:
  //! \brief Parses a data section from an XML document.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const tinyxml2::XMLElement* elem);

  //! \brief Parses a dimension-specific subelement of the \a geometry XML-tag.
  virtual bool parseGeometryDimTag(const tinyxml2::XMLElement* elem)
  {
    return this->parseGeometryTag(elem);
  }

  //! \brief Parses or generates app-specific explicit knots for refinement.
  virtual bool parseXi(const tinyxml2::XMLElement*, RealArray&) const
  { return false; }

  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);

  //! \brief Preprocesses a user-defined Dirichlet boundary property.
  //! \param[in] patch 1-based index of the patch to receive the property
  //! \param[in] lndx Local index of the boundary item to receive the property
  //! \param[in] ldim Dimension of the boundary item to receive the property
  //! \param[in] dirs Which local DOFs to constrain
  //! \param[in] code In-homogeneous Dirichlet condition property code
  //! \param ngnod Total number of global nodes in the model (might be updated)
  virtual bool addConstraint(int patch, int lndx, int ldim,
                             int dirs, int code, int& ngnod, char = 1);

  //! \brief Returns a FEM model generator for a default single-patch model.
  //! \param[in] geo XML element containing geometry definition
  virtual ModelGenerator* getModelGenerator(const tinyxml2::XMLElement* geo) const;

  //! \brief Reads a patch from given input stream.
  //! \param[in] isp The input stream to read from
  //! \param[in] pchInd 0-based index of the patch to read
  //! \param[in] unf Number of unknowns per basis function for each field
  //! \param[in] whiteSpace For message formatting
  virtual ASMbase* readPatch(std::istream& isp, int pchInd, const CharVec& unf,
                             const char* whiteSpace) const;

  //! \brief Connects two patches.
  //! \param[in] interface Patch interface definition
  virtual bool connectPatches(const ASM::Interface& interface, bool = true);

protected:
  unsigned char nf; //!< Number of scalar fields
};

#endif
