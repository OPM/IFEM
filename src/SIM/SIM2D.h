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

#include "SIMbase.h"


/*!
  \brief Driver class for 2D NURBS-based FEM solver.
  \details The class implements the parse method of the parent class, and can
  be used for any 2D continuum or shell problem.
*/

class SIM2D : public SIMbase
{
public:
  //! \brief Default constructor.
  //! \param[in] n1 Dimension of the primary solution field
  //! \param[in] n2 Dimension of the second solution field (mixed method)
  SIM2D(unsigned char n1 = 2, unsigned char n2 = 0);
  //! \brief Empty destructor.
  virtual ~SIM2D() {}

  //! \brief Defines the spatial numerical integration scheme to use.
  //! \param[in] ng Number of Gauss points in each parameter direction
  virtual void setQuadratureRule(size_t ng);

  //! \brief Creates the FE model by copying the given patches.
  //! \param[in] patches List of patches to borrow the grid from
  //! \param[in] g2ln Global-to-local node number mapping for the borrowed grid
  void clonePatches(const FEModelVec& patches, const std::map<int,int>& g2ln);

protected:
  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);

  //! \brief Parses a data section from an XML document.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Reads patches from given input stream.
  //! \param[in] isp The input stream to read from
  virtual bool readPatches(std::istream& isp);
  //! \brief Reads a patch from given input stream.
  //! \param[in] isp The input stream to read from
  //! \param[in] pchInd 0-based index of the patch to read
  virtual bool readPatch(std::istream& isp, int pchInd);
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
  void readNodes(std::istream& isn);

  //! \brief Refines a list of elements.
  //! \param[in] elements 1-based indices of the elements to refine
  //! \param[in] options Input options to refinement algorithm
  //! \param[in] fName Optional mesh output file (Encapsulated PostScript)
  virtual bool refine(const std::vector<int>& elements,
		      const std::vector<int>& options,
		      const char* fName = 0);

  //! \brief Preprocesses a user-defined Dirichlet boundary property.
  //! \param[in] patch 1-based index of the patch to receive the property
  //! \param[in] lndx Local index of the boundary item to receive the property
  //! \param[in] ldim Dimension of the boundary item to receive the property
  //! \param[in] dirs Which local DOFs to constrain
  //! \param[in] code In-homegeneous Dirichlet condition property code
  virtual bool addConstraint(int patch, int lndx, int ldim,
			     int dirs, int code = 0);

private:
  bool isRefined; //!< If \e true, the model has been adaptively refined
  //! \brief Parse subtags of the <geometry> block from the input file
  bool parseGeometryTag(const TiXmlElement* elem);

protected:
  unsigned char nf[2]; //!< Number of scalar fields
};

#endif
