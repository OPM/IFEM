// $Id$
//==============================================================================
//!
//! \file ASM1D.h
//!
//! \date May 15 2014
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Abstract interface for 1D patches.
//!
//==============================================================================

#ifndef _ASM_1D_H
#define _ASM_1D_H

#include "ASMenums.h"
#include <vector>
#include <cstddef>

class ASMbase;
class Vec3;


/*!
  \brief Abstract interface for 1D spline patches.
  \details This class contains an interface to methods common for 1D patches.
*/

class ASM1D
{
protected:
  //! \brief The constructor is protected to allow objects of sub-classes only.
  ASM1D() {}

public:
  //! \brief Empty destructor.
  virtual ~ASM1D() {}

  //! \brief Creates a one-parametric patch of specified discretization type.
  //! \param[in] type The discretization method to use
  //! \param[in] nd Number of spatial dimensions
  //! \param[in] nf Number of unknowns per basis function in the patch
  static ASMbase* create(ASM::Discretization type,
                         unsigned char nd, unsigned char nf);
  //! \brief Creates a one-parametric patch of specified discretization type.
  //! \param[in] type The discretization method to use
  //! \param[in] nf Number of unknowns per basis function in the patch
  static ASMbase* create(ASM::Discretization type, unsigned char nf = 1);

  //! \brief Returns a copy of this patch with identical FE discretization.
  //! \param[in] nf Number of unknown per basis function in the patch
  //!
  //! \note The copied patch shares the FE data structures with the copy,
  //! in order to save memory. Thus, the copy cannot be read from file, refined,
  //! or changed in other ways that affect the FE geometry and/or topology.
  //! The other properties of the patch (boundary conditions, constraints,
  //! loads, etc.) are however not copied.
  ASMbase* clone(unsigned char* nf = nullptr) const;

  //! \brief Refines the parametrization by inserting extra knots uniformly.
  //! \param[in] nInsert Number of extra knots to insert in each knot-span
  virtual bool uniformRefine(int nInsert) = 0;
  //! \brief Raises the order of the spline object for this patch.
  //! \param[in] ru Number of times to raise the order in u-direction
  virtual bool raiseOrder(int ru) = 0;
  //! \brief Refines the parametrization by inserting extra knots.
  //! \param[in] xi Relative positions of added knots in each existing knot span
  virtual bool refine(const std::vector<double>& xi) = 0;

  //! \brief Constrains a node identified by two relative parameter values.
  //! \param[in] xi Parameter value along the curve
  //! \param[in] dof Which DOFs to constrain at the node
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] basis Basis to constrain node for
  //! \return 1-based index of the constrained node
  //!
  //! \details The parameter value has to be in the domain [0.0,1.0], where
  //! 0.0 means the beginning of the domain and 1.0 means the end. For values
  //! in between, the actual index is taken as the integer value closest to
  //! \a r*n, where \a r denotes the given relative parameter value,
  //! and \a n is the number of nodes along that parameter direction.
  virtual int constrainNode(double xi, int dof,
                            int code = 0, char basis = 1) = 0;

  //! \brief Calculates parameter values for visualization nodal points.
  //! \param[out] prm Parameter values for all points
  //! \param[in] nSegSpan Number of visualization segments over each knot-span
  virtual bool getGridParameters(std::vector<double>& prm,
                                 int nSegSpan) const = 0;

  //! \brief Returns characteristic element size based on end point coordinates.
  static double getElementSize(const std::vector<Vec3>& XC);
};

#endif
