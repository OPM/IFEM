// $Id: SIM1D.h,v 1.1 2010-09-02 15:28:01 kmo Exp $
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

#include "SIMbase.h"


/*!
  \brief Driver class for 1D NURBS-based FEM solver.
*/

class SIM1D : public SIMbase
{
public:
  //! \brief Default constructor.
  //! \param[in] n_f Number of components in the primary solution field
  SIM1D(unsigned char n_f = 1);
  //! \brief Empty destructor.
  virtual ~SIM1D() {}

  //! \brief Defines the spatial numerical integration scheme to use.
  //! \param[in] ng Number of Gauss points in each parameter direction
  virtual void setQuadratureRule(size_t ng);

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
  unsigned char nf; //!< Number of scalar fields
};

#endif
