// $Id$
//==============================================================================
//!
//! \file ASM3D.h
//!
//! \date January 2013
//!
//! \author Kjetil A. Johannessen / NTNU
//!
//! \brief Abstract interface for 3D patches.
//!
//==============================================================================

#ifndef _ASM_3D_H
#define _ASM_3D_H

#include "ASMenums.h"
#include <vector>
#include <cstddef>

class ASMbase;


/*!
  \brief Abstract interface for 3D spline patches.
  \details This class contains an interface to methods common for structured and
  unstructured 3D patches, such that these methods can be invoked without the
  need to type-cast the patch object to the actual class type.
*/

class ASM3D
{
protected:
  //! \brief The constructor is protected to allow objects of sub-classes only.
  ASM3D() {}

public:
  //! \brief Empty destructor.
  virtual ~ASM3D() {}

  typedef std::vector<unsigned char> CharVec; //!< Convenience type

  //! \brief Creates a three-parametric patch of specified discretization type.
  //! \param[in] type The discretization method to use
  //! \param[in] nf Number of unknowns per basis function in the patch
  //! \param[in] mixedFEM If \e true, force mixed formulation even if \a nf[1]=0
  static ASMbase* create(ASM::Discretization type,
                         const CharVec& nf, bool mixedFEM = false);

  //! \brief Creates a three-parametric patch of specified discretization type.
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
  ASMbase* clone(const CharVec& nf = CharVec()) const;

  //! \brief Checks that the patch is modelled in a right-hand-side system.
  virtual bool checkRightHandSystem() = 0;

  //! \brief Refines the parametrization by inserting extra knots uniformly.
  //! \param[in] dir Parameter direction to refine
  //! \param[in] nInsert Number of extra knots to insert in each knot-span
  virtual bool uniformRefine(int dir, int nInsert) = 0;
  //! \brief Raises the order of the spline object for this patch.
  //! \param[in] ru Number of times to raise the order in u-direction
  //! \param[in] rv Number of times to raise the order in v-direction
  //! \param[in] rw Number of times to raise the order in w-direction
  virtual bool raiseOrder(int ru, int rv, int rw) = 0;
  //! \brief Refines the parametrization by inserting extra knots.
  //! \param[in] dir Parameter direction to refine
  //! \param[in] xi Relative positions of added knots in each existing knot span
  virtual bool refine(int dir, const std::vector<double>& xi) = 0;

  //! \brief Constrains all DOFs on a given boundary face.
  //! \param[in] dir Parameter direction defining the face to constrain
  //! \param[in] open If \e true, exclude all points along the face boundary
  //! \param[in] dof Which DOFs to constrain at each node on the face
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] basis Which basis to constrain face for
  virtual void constrainFace(int dir, bool open, int dof,
                             int code = 0, char basis = 1) = 0;
  //! \brief Constrains all DOFs in local directions on a given boundary face.
  //! \param[in] dir Parameter direction defining the face to constrain
  //! \param[in] open If \e true, exclude all points along the face boundary
  //! \param[in] dof Which local DOFs to constrain at each node on the face
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] project If \e true, the local axis directions are projected
  //! \param[in] T1 Desired global direction of first local tangent direction
  //! \return Number of additional nodes added due to local axis constraints
  virtual size_t constrainFaceLocal(int dir, bool open, int dof, int code = 0,
                                    bool project = false, char T1 = '\0') = 0;

  //! \brief Constrains all DOFs on a given boundary edge.
  //! \param[in] lEdge Local index [1,12] of the edge to constrain
  //! \param[in] open If \e true, exclude the end points of the edge
  //! \param[in] dof Which DOFs to constrain at each node on the edge
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] basis Which basis to constrain edge for
  virtual void constrainEdge(int lEdge, bool open, int dof, int code = 0,
                             char basis = 1) = 0;
  //! \brief Constrains all DOFs along a line on a given boundary face.
  //! \param[in] fdir Parameter direction defining the face to constrain
  //! \param[in] ldir Parameter direction defining the line to constrain
  //! \param[in] xi Parameter value defining the line to constrain
  //! \param[in] dof Which DOFs to constrain at each node along the line
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] basis Which basis to constrain line for
  //!
  //! \details The parameter \a xi has to be in the domain [0.0,1.0], where
  //! 0.0 means the beginning of the domain and 1.0 means the end. The line to
  //! constrain goes along the parameter direction \a ldir in the face with
  //! with normal in parameter direction \a fdir, and positioned along the third
  //! parameter direction as indicated by \a xi. The actual value of \a xi
  //! is converted to the integer value closest to \a xi*n, where \a n is the
  //! number of nodes (control points) in that parameter direction.
  virtual void constrainLine(int fdir, int ldir, double xi, int dof = 123,
                             int code = 0, char basis = 1) = 0;

  //! \brief Constrains a corner node identified by the three parameter indices.
  //! \param[in] I Parameter index in u-direction
  //! \param[in] J Parameter index in v-direction
  //! \param[in] K Parameter index in w-direction
  //! \param[in] dof Which DOFs to constrain at the node
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] basis Which basis to constrain corner for
  //!
  //! \details The sign of the three indices is used to define whether we want
  //! the node at the beginning or the end of that parameter direction.
  //! The magnitude of the indices are not used.
  virtual void constrainCorner(int I, int J, int K, int dof = 123,
                               int code = 0, char basis = 1) = 0;
  //! \brief Constrains a node identified by three relative parameter values.
  //! \param[in] xi Parameter in u-direction
  //! \param[in] eta Parameter in v-direction
  //! \param[in] zeta Parameter in w-direction
  //! \param[in] dof Which DOFs to constrain at the node
  //! \param[in] code Inhomogeneous dirichlet condition code
  //! \param[in] basis Which basis to constrain node for
  //!
  //! \details The parameter values have to be in the domain [0.0,1.0], where
  //! 0.0 means the beginning of the domain and 1.0 means the end. For values
  //! in between, the actual index is taken as the integer value closest to
  //! \a r*n, where \a r denotes the given relative parameter value,
  //! and \a n is the number of nodes along that parameter direction.
  virtual void constrainNode(double xi, double eta, double zeta,
                             int dof = 123, int code = 0, char basis = 1) = 0;

  //! \brief Calculates parameter values for visualization nodal points.
  //! \param[out] prm Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1)
  //! \param[in] nSegSpan Number of visualization segments over each knot-span
  virtual bool getGridParameters(std::vector<double>& prm,
                                 int dir, int nSegSpan) const = 0;

  //! \brief Returns the node index for a given corner.
  virtual int getCorner(int I, int J, int K, int basis = 1) const = 0;

  //! \brief Returns the node indices for a given edge.
  virtual std::vector<int> getEdge(int lEdge, bool open, int basis = 1) const = 0;
};

#endif
