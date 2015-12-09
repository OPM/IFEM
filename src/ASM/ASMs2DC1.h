// $Id$
//==============================================================================
//!
//! \file ASMs2DC1.h
//!
//! \date Oct 25 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of C1-continuous structured 2D spline FE models.
//!
//==============================================================================

#ifndef _ASM_S2D_C1_H
#define _ASM_S2D_C1_H

#include "ASMs2D.h"


/*!
  \brief Driver for assembly of C1-continuous structured 2D spline FE models.
  \details This class extends the ASMs2D class to handle C1-continuity
  over patch interfaces, as well as boundary conditions on derivatives.
*/

class ASMs2DC1 : public ASMs2D
{
public:
  //! \brief Default constructor.
  ASMs2DC1(unsigned char n_s = 2, unsigned char n_f = 1) : ASMs2D(n_s,n_f) {}
  //! \brief Copy constructor.
  ASMs2DC1(const ASMs2DC1& patch, unsigned char n_f = 0) : ASMs2D(patch,n_f) {}
  //! \brief Empty destructor.
  virtual ~ASMs2DC1() {}


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
  //! \param[in] basis Basis to constrain
  virtual void constrainEdge(int dir, bool open, int dof = 12, int code = 0,
                             char basis = 1);

  using ASMs2D::constrainCorner;
  //! \brief Constrains a corner node identified by the two parameter indices.
  //! \param[in] I Parameter index in u-direction
  //! \param[in] J Parameter index in v-direction
  //! \param[in] dof Which DOFs to constrain at the node
  //! \param[in] code Inhomogeneous dirichlet condition code
  //!
  //! \details The sign of the two indices is used to define whether we want
  //! the node at the beginning or the end of that parameter direction.
  //! The magnitude of the indices are not used.
  virtual void constrainCorner(int I, int J, int dof = 12, int code = 0);
  //! \brief Constrains a node identified by two relative parameter values.
  //! \param[in] xi Parameter in u-direction
  //! \param[in] eta Parameter in v-direction
  //! \param[in] dof Which DOFs to constrain at the node
  //! \param[in] code Inhomogeneous dirichlet condition code
  //!
  //! \details The parameter values have to be in the domain [0.0,1.0], where
  //! 0.0 means the beginning of the domain and 1.0 means the end. For values
  //! in between, the actual index is taken as the integer value closest to
  //! \a r*n, where \a r denotes the given relative parameter value,
  //! and \a n is the number of nodes along that parameter direction.
  virtual void constrainNode(double xi, double eta, int dof = 12, int code = 0);

  //! \brief Connects all matching nodes on two adjacent boundary edges.
  //! \param[in] edge Local edge index of this patch, in range [1,4]
  //! \param neighbor The neighbor patch
  //! \param[in] nedge Local edge index of neighbor patch, in range [1,4]
  //! \param[in] revers Indicates whether the two edges have opposite directions
  bool connectC1(int edge, ASMs2DC1* neighbor, int nedge, bool revers = false);

  //! \brief Makes two opposite boundary edges periodic.
  //! \param[in] dir Parameter direction defining the periodic edges
  //! \param[in] basis Which basis to connect (mixed methods), 0 means both
  //! \param[in] master 1-based index of the first master node in this basis
  virtual void closeEdges(int dir, int basis = 0, int master = 1);

  //! \brief Renumbers the global node numbers in the \a neighbors map.
  //! \param[in] old2new Old-to-new node number mapping
  static void renumberNodes(const std::map<int,int>& old2new);

  //! \brief Initializes constraint equations enforcing C1-continuity.
  virtual bool initConstraints();

  //! \brief Updates the time-dependent in-homogeneous Dirichlet coefficients.
  //! \param[in] func Scalar property fields
  //! \param[in] vfunc Vector property fields
  //! \param[in] time Current time
  //! \param[in] g2l Pointer to global-to-local node number mapping
  virtual bool updateDirichlet(const std::map<int,RealFunc*>& func,
                               const std::map<int,VecFunc*>& vfunc, double time,
                               const std::map<int,int>* g2l = nullptr);

private:
  static std::map<int,ASMs2DC1*> neighbors; //!< Global node to patch mapping
};

#endif
