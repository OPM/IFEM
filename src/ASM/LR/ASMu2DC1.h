// $Id$
//==============================================================================
//!
//! \file ASMu2DC1.h
//!
//! \date Oct 5 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of C1-continuous 2D LR-spline FE models.
//!
//==============================================================================

#ifndef _ASM_U2D_C1_H
#define _ASM_U2D_C1_H

#include "ASMu2D.h"


/*!
  \brief Driver for assembly of C1-continuous 2D LR-spline FE models.
  \details This class extends the ASMu2D class to handle C1-continuity
  over patch interfaces, as well as boundary conditions on derivatives.
*/

class ASMu2DC1 : public ASMu2D
{
public:
  //! \brief Default constructor.
  ASMu2DC1(unsigned char n_s = 2, unsigned char n_f = 1) : ASMu2D(n_s,n_f) {}
  //! \brief Copy constructor.
  ASMu2DC1(const ASMu2DC1& patch, unsigned char n_f = 0) : ASMu2D(patch,n_f) {}
  //! \brief Empty destructor.
  virtual ~ASMu2DC1() {}


  // Methods for model generation
  // ============================

  //! \brief Generates the finite element topology data for the patch.
  //! \details This method is reimplemented to check that the patch has
  //! sufficient polynomial order.
  virtual bool generateFEMTopology();


  // Various methods for preprocessing of boundary conditions and patch topology
  // ===========================================================================

  //! \brief Constrains all DOFs on a given boundary edge.
  //! \param[in] dir Parameter direction defining the edge to constrain
  //! \param[in] open If \e true, exclude the end points of the edge
  //! \param[in] dof Which DOFs to constrain at each node on the edge
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] basis Which basis to constrain edge for
  virtual void constrainEdge(int dir, bool open, int dof, int code, char basis);
};

#endif
