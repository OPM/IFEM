// $Id$
//==============================================================================
//!
//! \file SIM3D.h
//!
//! \date Dec 08 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for 3D NURBS-based FEM analysis.
//!
//==============================================================================

#ifndef _SIM_3D_H
#define _SIM_3D_H

#include "SIMbase.h"


/*!
  \brief Driver class for 3D NURBS-based FEM solver.
  \details The class implements the parse method of the parent class, and can
  be used for any 3D continuum problem.
*/

class SIM3D : public SIMbase
{
public:
  //! \brief Default constructor.
  //! \param[in] check If \e true, ensure the model is in a right-hand system
  //! \param[in] n1 Dimension of the primary solution field
  //! \param[in] n2 Dimension of the second solution field (mixed method)
  SIM3D(bool check = false, unsigned char n1 = 3, unsigned char n2 = 0);
  //! \brief Empty destructor.
  virtual ~SIM3D() {}

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

  //! \brief Constrains a parametric line on a boundary face.
  //! \param[in] patch 1-based index of the patch to receive the property
  //! \param[in] lndx Local index of the boundary face to receive the property
  //! \param[in] line Local direction of the line on the face (1=I, 2=J)
  //! \param[in] xi Relative coordinate [0,1] defining the line placement
  //! \param[in] dirs Which local DOFs to constrain
  bool addConstraint(int patch, int lndx, int line, double xi, int dirs);

protected:
  unsigned char nf[2]; //!< Number of scalar fields

  bool checkRHSys;  //!< Check if all patches are in a right-hand system
};

#endif
