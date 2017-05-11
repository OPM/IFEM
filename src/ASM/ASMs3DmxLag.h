// $Id$
//==============================================================================
//!
//! \file ASMs3DmxLag.h
//!
//! \date Dec 28 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of structured 3D Lagrange mixed FE models.
//!
//==============================================================================

#ifndef _ASM_S3D_MX_LAG_H
#define _ASM_S3D_MX_LAG_H

#include "ASMs3DLag.h"
#include "ASMmxBase.h"
#include <array>


/*!
  \brief Driver for assembly of structured 3D Lagrange mixed FE models.
  \details This class implements a two-field mixed formulation with Lagrangian
  basis functions. The geometry and the first field are of equal order and
  is one order higher than the second field.
*/

class ASMs3DmxLag : public ASMs3DLag, private ASMmxBase
{
public:
  //! \brief The constructor initializes the dimension of each basis.
  ASMs3DmxLag(const CharVec& n_f);
  //! \brief Copy constructor.
  ASMs3DmxLag(const ASMs3DmxLag& patch, const CharVec& n_f = CharVec(2,0));
  //! \brief Empty destructor.
  virtual ~ASMs3DmxLag() {}


  // Methods for model generation
  // ============================

  //! \brief Generates the finite element topology data for the patch.
  //! \details The data generated are the element-to-node connectivity array,
  //! the nodal coordinate array, as well as global node and element numbers.
  virtual bool generateFEMTopology();

  //! \brief Clears the contents of the patch, making it empty.
  //! \param[in] retainGeometry If \e true, the spline geometry is not cleared.
  //! This is used to reinitialize the patch after it has been refined.
  virtual void clear(bool retainGeometry = false);

  //! \brief Returns the total number of nodes in this patch.
  virtual size_t getNoNodes(int basis = 0) const;
  //! \brief Returns the number of solution fields.
  virtual unsigned char getNoFields(int basis = 0) const;
  //! \brief Returns the number of DOFs per node.
  //! \param[in] inod 1-based node index local to current patch
  virtual unsigned char getNodalDOFs(size_t inod) const;
  //! \brief Returns the classification of a node.
  //! \param[in] inod 1-based node index local to current patch
  virtual char getNodeType(size_t inod) const;

  //! \brief Initializes the patch level MADOF array for mixed problems.
  virtual void initMADOF(const int* sysMadof);

  //! \brief Connects all matching nodes on two adjacent boundary faces.
  //! \param[in] face Local face index of this patch, in range [1,6]
  //! \param neighbor The neighbor patch
  //! \param[in] nface Local face index of neighbor patch, in range [1,6]
  //! \param[in] norient Relative face orientation flag (see class ASMs3D)
  //! \param[in] basis Which basis to connect, or 0 for all.
  //! \param[in] coordCheck False to disable coordinate checks (periodic connections)
  //! \param[in] thick Thickness of connection
  virtual bool connectPatch(int face, ASMs3D& neighbor, int nface, int norient,
                            int basis = 0, bool coordCheck = true, int thick = 1);

  //! \brief Makes two opposite boundary faces periodic.
  //! \param[in] dir Parameter direction defining the periodic faces
  virtual void closeFaces(int dir, int = 0, int = 1);


  // Methods for integration of finite element quantities.
  // These are the main computational methods of the ASM class hierarchy.
  // ====================================================================

  //! \brief Evaluates an integral over the interior patch domain.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrate(Integrand& integrand,
			 GlobalIntegral& glbInt, const TimeDomain& time);

  //! \brief Evaluates a boundary integral over a patch face.
  //! \param integrand Object with problem-specific data and methods
  //! \param[in] lIndex Local index [1,6] of the boundary face
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrate(Integrand& integrand, int lIndex,
			 GlobalIntegral& glbInt, const TimeDomain& time);


  // Post-processing methods
  // =======================

  //! \brief Extract the primary solution field at the specified nodes.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] nodes 1-based local node numbers to extract solution for
  virtual bool getSolution(Matrix& sField, const Vector& locSol,
			   const IntVec& nodes) const;

  using ASMs3DLag::evalSolution;
  //! \brief Evaluates the primary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] gpar Parameter values of the result sampling points
  //! \param[in] regular Flag indicating how the sampling points are defined
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const RealArray* gpar, bool regular = true,
                            int = 0) const;

  //! \brief Evaluates the secondary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] gpar Parameter values of the result sampling points
  //! \param[in] regular Flag indicating how the sampling points are defined
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integrand,
			    const RealArray* gpar, bool regular = true) const;

  //! \brief Extracts nodal results for this patch from the global vector.
  //! \param[in] globVec Global solution vector in DOF-order
  //! \param[out] nodeVec Nodal result vector for this patch
  //! \param[in] basis Which basis (or 0 for both) to extract nodal values for
  virtual void extractNodeVec(const Vector& globVec, Vector& nodeVec,
			      unsigned char = 0, int basis = 0) const;

protected:

  // Internal utility methods
  // ========================

  //! \brief Returns the number of nodal points in each parameter direction.
  //! \param[out] n1 Number of nodes in first (u) direction
  //! \param[out] n2 Number of nodes in second (v) direction
  //! \param[out] n3 Number of nodes in third (w) direction
  //! \param[in] basis Which basis to return size parameters for
  virtual bool getSize(int& n1, int& n2, int& n3, int basis = 0) const;

private:
  std::vector<size_t> nxx; //!< Basis dimensions in 1st parameter direction
  std::vector<size_t> nyx; //!< Basis dimensions in 2nd parameter direction
  std::vector<size_t> nzx; //!< Basis dimensions in 3rd parameter direction

  std::vector< std::array<size_t,3> > elem_sizes; //!< Size on each basis
};

#endif
