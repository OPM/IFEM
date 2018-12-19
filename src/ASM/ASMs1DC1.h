// $Id$
//==============================================================================
//!
//! \file ASMs1DC1.h
//!
//! \date Oct 25 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of C1-continuous structured 1D spline FE models.
//!
//==============================================================================

#ifndef _ASM_S1D_C1_H
#define _ASM_S1D_C1_H

#include "ASMs1D.h"


/*!
  \brief Driver for assembly of C1-continuous structured 1D spline FE models.
  \details This class extends the ASMs1D class to handle C1-continuity
  over patch interfaces, as well as boundary conditions on derivatives.
*/

class ASMs1DC1 : public ASMs1D
{
public:
  //! \brief Default constructor.
  ASMs1DC1(unsigned char n_s = 1, unsigned char n_f = 1) : ASMs1D(n_s,n_f) {}
  //! \brief Copy constructor.
  ASMs1DC1(const ASMs1DC1& patch, unsigned char n_f = 0) : ASMs1D(patch,n_f) {}
  //! \brief Empty destructor.
  virtual ~ASMs1DC1() {}


  // Methods for model generation
  // ============================

  //! \brief Generates the finite element topology data for the patch.
  //! \details This method is reimplemented to check that the patch has
  //! sufficient polynomial order.
  virtual bool generateFEMTopology();


  // Various methods for preprocessing of boundary conditions and patch topology
  // ===========================================================================

  //! \brief Constrains a node identified by two relative parameter values.
  //! \param[in] xi Parameter value along the curve
  //! \param[in] dof Which DOFs to constrain at the node
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \return 1-based index of the constrained node
  //!
  //! \details The parameter value has to be in the domain [0.0,1.0], where
  //! 0.0 means the beginning of the domain and 1.0 means the end. For values
  //! in between, the actual index is taken as the integer value closest to
  //! \a r*n, where \a r denotes the given relative parameter value,
  //! and \a n is the number of nodes along that parameter direction.
  virtual int constrainNode(double xi, int dof, int code);
  //! \brief Constrains all DOFs in local directions at a given end point.
  //! \param[in] dir Parameter direction defining the end to constrain
  //! \param[in] dof Which local DOFs to constrain at the end point
  //! \param[in] code Inhomogeneous dirichlet condition code
  virtual size_t constrainEndLocal(int dir, int dof, int code);
};

#endif
