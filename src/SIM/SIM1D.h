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
  SIM1D(unsigned char n1 = 1, bool = false);
  //! \brief Constructor used for mixed problems.
  //! \param[in] unf Dimension of the primary solution fields
  SIM1D(const CharVec& unf, bool = false);
  //! \brief Constructor that also initializes the integrand pointer.
  //! \param[in] itg Pointer to the integrand of the problem to solve
  //! \param[in] n Dimension of the primary solution field
  SIM1D(IntegrandBase* itg, unsigned char n = 1);
  //! \brief The destructor deletes the twist angle function.
  virtual ~SIM1D() { delete twist; }

  //! \brief Returns the number of parameter dimensions in the model.
  unsigned short int getNoParamDim() const { return 1; }

  //! \brief Reads a patch from given input stream.
  //! \param[in] isp The input stream to read from
  //! \param[in] pchInd 0-based index of the patch to read
  //! \param[in] unf Number of unknowns per basis function for each field
  virtual ASMbase* readPatch(std::istream& isp, int pchInd,
                             const CharVec& unf) const;

  //! \brief Evaluates the primary solution at the given point.
  //! \param[in] psol Primary solution vector
  //! \param[in] u Parameter of the point to evaluate at
  //! \param[in] deriv Derivative order of the solution
  //! \param[in] patch 1-based patch index contining the evaluation point
  //! \return Evaluated solution values
  Vector getSolution(const Vector& psol, double u,
                     int deriv = 0, int patch = 1) const;

  //! \brief Updates the nodal rotations for problems with rotational DOFs.
  //! \param[in] incSol Incremental solution to update the rotations with
  //! \param[in] alpha Scaling factor for the incremental solution.
  //! If 0.0, reinitialize the rotations from unity
  virtual bool updateRotations(const Vector& incSol, double alpha);

private:
  //! \brief Parses a subelement of the \a geometry XML-tag.
  bool parseGeometryTag(const TiXmlElement* elem);
  //! \brief Parses a subelement of the \a boundaryconditions XML-tag.
  bool parseBCTag(const TiXmlElement* elem);

protected:
  //! \brief Parses the twist angle description along the curve.
  bool parseTwist(const TiXmlElement* elem);

  //! \brief Parses a data section from an XML document.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);

  //! \brief Reads patches from given input stream.
  //! \param[in] isp The input stream to read from
  //! \param[out] patches Array of patches that were read
  //! \param[in] whiteSpace For message formatting
  virtual bool readPatches(std::istream& isp, PatchVec& patches,
                           const char* whiteSpace);

  //! \brief Preprocesses a user-defined Dirichlet boundary property.
  //! \param[in] patch 1-based index of the patch to receive the property
  //! \param[in] lndx Local index of the boundary item to receive the property
  //! \param[in] ldim Dimension of the boundary item to receive the property
  //! \param[in] dirs Which local DOFs to constrain
  //! \param[in] code In-homogeneous Dirichlet condition property code
  //! \param[in] basis Which basis to apply the constraint to (mixed methods)
  virtual bool addConstraint(int patch, int lndx, int ldim,
                             int dirs, int code, int&, char basis = 1);

  //! \brief Creates the computational FEM model from the spline patches.
  //! \details Reimplemented to account for twist angle in beam problems.
  virtual bool createFEMmodel(char = 'y');

  //! \brief Creates G2 representation of a line.
  //! \param geo XML block with geometry definition
  std::string createDefaultG2(const TiXmlElement* geo) const;

  //! \brief Creates a default single-patch geometry.
  virtual ASMbase* createDefaultGeometry(const TiXmlElement* geo) const;

protected:
  unsigned char nf; //!< Number of scalar fields

private:
  RealFunc* twist; //!< Twist angle (for 3D beam problems only)
  Vec3      XZp;   //!< Local Z-direction vector
};

#endif
