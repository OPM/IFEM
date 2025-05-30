// $Id$
//==============================================================================
//!
//! \file SIM2D.h
//!
//! \date Feb 04 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for 2D NURBS-based FEM analysis.
//!
//==============================================================================

#ifndef _SIM_2D_H
#define _SIM_2D_H

#include "SIMgeneric.h"


/*!
  \brief Driver class for 2D NURBS-based FEM solver.
  \details The class implements the parse method of the parent class,
  and can be used for any 2D continuum or shell problem.
*/

class SIM2D : public SIMgeneric
{
public:
  //! \brief Enum announcing the dimensionality (used for template writing).
  enum { dimension = 2 };

  //! \brief Default constructor.
  //! \param[in] n1 Dimension of the primary solution field
  //! \param[in] check If \e true, ensure the model is in a right-hand system
  explicit SIM2D(unsigned char n1 = 2, bool check = false);
  //! \brief Constructor used for mixed problems.
  //! \param[in] unf Dimension of the primary solution fields
  //! \param[in] check If \e true, ensure the model is in a right-hand system
  explicit SIM2D(const CharVec& unf, bool check = false);
  //! \brief Constructor that also initializes the integrand pointer.
  //! \param[in] itg Pointer to the integrand of the problem to solve
  //! \param[in] n Dimension of the primary solution field
  //! \param[in] check If \e true, ensure the model is in a right-hand system
  explicit SIM2D(IntegrandBase* itg, unsigned char n = 2, bool check = false);
  //! \brief Empty destructor.
  virtual ~SIM2D() {}

  //! \brief Returns whether a mixed formulation is used (used by HDF5 output).
  virtual bool mixedProblem() const { return nf.size() > 1; }

  //! \brief Returns the number of parameter dimensions in the model.
  virtual unsigned short int getNoParamDim() const { return 2; }

  //! \brief Creates the FE model by copying the given patches.
  //! \param[in] patches List of patches to borrow the grid from
  //! \param[in] g2ln Global-to-local node number mapping for the borrowed grid
  virtual void clonePatches(const PatchVec& patches,
                            const std::map<int,int>& g2ln);

  using SIMgeneric::getSolution;
  //! \brief Evaluates the primary solution at the given point.
  //! \param[in] psol Primary solution vector
  //! \param[in] u First parameter of the point to evaluate at
  //! \param[in] v Second parameter of the point to evaluate at
  //! \param[in] deriv Derivative order of the solution
  //! \param[in] patch 1-based patch index contining the evaluation point
  //! \return Evaluated solution values
  RealArray getSolution(const Vector& psol, double u, double v,
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

  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);

  //! \brief Parses a dimension-specific subelement of the \a geometry XML-tag.
  virtual bool parseGeometryDimTag(const tinyxml2::XMLElement* elem)
  {
    return this->parseGeometryTag(elem);
  }

  //! \brief Reads global node data for a patch from given input stream.
  //! \param[in] isn The input stream to read from
  //! \param[in] pchInd 0-based index of the patch to read node data for
  //! \param[in] basis The basis to read node data for (mixed FEM)
  //! \param[in] oneBased If \e true the read node numbers are assumed
  //! one-based. If \e false they are assumed to be zero-based.
  virtual bool readNodes(std::istream& isn, int pchInd, int basis = 0,
                         bool oneBased = false);

  //! \brief Reads node numbers from given input stream.
  //! \param[in] isn The file stream to read from
  virtual void readNodes(std::istream& isn);

  //! \brief Preprocesses a user-defined Dirichlet boundary property.
  //! \param[in] patch 1-based index of the patch to receive the property
  //! \param[in] lndx Local index of the boundary item to receive the property
  //! \param[in] ldim Dimension of the boundary item to receive the property
  //! \param[in] dirs Which local DOFs to constrain
  //! \param[in] code In-homogeneous Dirichlet condition property code
  //! \param ngnod Total number of global nodes in the model (might be updated)
  //! \param[in] basis Which basis to apply the constraint to (mixed methods)
  //! \param[in] ovrD If \e true, override conflicting Dirichlet conditions
  virtual bool addConstraint(int patch, int lndx, int ldim, int dirs, int code,
                             int& ngnod, char basis = 1, bool ovrD = false);

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
  //! \param[in] ifc Patch interface definition
  //! \param[in] coordCheck If \e false, do not check for matching coordinates
  virtual bool connectPatches(const ASM::Interface& ifc,
                              bool coordCheck = true);

  //! \brief Writes out the additional functions to VTF-file.
  virtual bool writeAddFuncs(int iStep, int& nBlock, int idBlock, double time);

private:
  std::map<std::string,RealFunc*> myAddScalars; //!< Scalar functions to output

protected:
  CharVec nf;         //!< Number of scalar fields
  bool    checkRHSys; //!< Check if all patches are in a right-hand system
};

#endif
