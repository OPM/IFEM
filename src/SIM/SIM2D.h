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

  //! \brief Read patches from given stream
  //! \param[in] isp The stream to read from
  void readPatches(std::istream& isp);

protected:
  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);

  //! \brief Preprocesses a user-defined Dirichlet boundary property.
  //! \param[in] patch 1-based index of the patch to receive the property
  //! \param[in] lndx Local index of the boundary item to receive the property
  //! \param[in] ldim Dimension of the boundary item to receive the property
  //! \param[in] dirs Which local DOFs to constrain
  //! \param[in] code In-homegeneous Dirichlet condition property code
  virtual bool addConstraint(int patch, int lndx, int ldim,
			     int dirs, int code = 0);

protected:
  unsigned char nf[2]; //!< Number of scalar fields
};

#endif
