// $Id$
//==============================================================================
//!
//! \file ASMs3Dmx.h
//!
//! \date Dec 28 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of structured 3D spline mixed FE models.
//!
//==============================================================================

#ifndef _ASM_S3D_MX_H
#define _ASM_S3D_MX_H

#include "ASMs3D.h"
#include "ASMmxBase.h"
#include <memory>


/*!
  \brief Driver for assembly of structured 3D spline mixed FE models.
  \details This class implements a two-field mixed formulation with splines as
  basis functions. The first field is of one order higher than the second field,
  and its basis is obtained by order-elevating the input spline object once.
  By default, the geometry is represented by the second (lower order) basis,
  however, by negating the \a n_f1 argument to the constructor, the first
  basis is used instead.
*/

class ASMs3Dmx : public ASMs3D, private ASMmxBase
{
public:
  //! \brief The constructor initializes the dimension of each basis.
  ASMs3Dmx(const CharVec& n_f);
  //! \brief Copy constructor.
  ASMs3Dmx(const ASMs3Dmx& patch, const CharVec& n_f = CharVec(2,0));
  //! \brief Empty destructor.
  virtual ~ASMs3Dmx() {}

  //! \brief Returns the spline surface representing the basis of this patch.
  virtual Go::SplineVolume* getBasis(int basis = 1) const;
  //! \brief Returns the spline curve representing a boundary of this patch.
  //! \param[in] dir Parameter direction defining which boundary to return
  //! \param[in] basis The basis to get the boundary for
  virtual Go::SplineSurface* getBoundary(int dir, int basis = 1);


  // Methods for model generation
  // ============================

  //! \brief Generates the finite element topology data for the patch.
  //! \details The data generated are the element-to-node connectivity array,
  //! the node-to-IJ-index array, as well as global node and element numbers.
  virtual bool generateFEMTopology();

  //! \brief Clears the contents of the patch, making it empty.
  //! \param[in] retainGeometry If \e true, the spline geometry is not cleared.
  //! This is used to reinitialize the patch after it has been refined.
  virtual void clear(bool retainGeometry = false);

  //! \brief Returns a matrix with nodal coordinates for an element.
  //! \param[in] iel Element index
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in one element
  virtual bool getElementCoordinates(Matrix& X, int iel) const;

  //! \brief Returns the global coordinates for the given node.
  //! \param[in] inod 1-based node index local to current patch
  virtual Vec3 getCoord(size_t inod) const;

  //! \brief Writes the geometry/basis of the patch to given stream.
  virtual bool write(std::ostream& os, int basis = 0) const;

  //! \brief Returns the number of bases.
  virtual size_t getNoBasis() const { return m_basis.size(); }
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

  //! \brief Evaluates a boundary integral over a patch edge.
  //! \param integrand Object with problem-specific data and methods
  //! \param[in] lIndex Local index of the boundary edge
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrate(Integrand& integrand, int lIndex,
			 GlobalIntegral& glbInt, const TimeDomain& time);


  // Post-processing methods
  // =======================

  //! \brief Evaluates the geometry at a specified point.
  //! \param[in] xi Dimensionless parameters in range [0.0,1.0] of the point
  //! \param[out] param The (u,v,w) parameters of the point in knot-span domain
  //! \param[out] X The Cartesian coordinates of the point
  //! \return Local node number within the patch that matches the point, if any
  //! \return 0 if no node (control point) matches this point
  virtual int evalPoint(const double* xi, double* param, Vec3& X) const;

  //! \brief Extract the primary solution field at the specified nodes.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] nodes 1-based local node numbers to extract solution for
  virtual bool getSolution(Matrix& sField, const Vector& locSol,
			   const IntVec& nodes) const;

  //! \brief Evaluates the primary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] gpar Parameter values of the result sampling points
  //! \param[in] regular Flag indicating how the sampling points are defined
  //! \param[in] deriv Derivative order to return
  //!
  //! \details When \a regular is \e true, it is assumed that the parameter
  //! value array \a gpar forms a regular tensor-product point grid of dimension
  //! \a gpar[0].size() \a X \a gpar[1].size() \a X \a gpar[2].size().
  //! Otherwise, we assume that it contains the \a u, \a v and \a w parameters
  //! directly for each sampling point.
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const RealArray* gpar, bool regular = true,
                            int deriv = 0) const;

  //! \brief Evaluates the secondary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] gpar Parameter values of the result sampling points
  //! \param[in] regular Flag indicating how the sampling points are defined
  //!
  //! \details The secondary solution is derived from the primary solution,
  //! which is assumed to be stored within the \a integrand for current patch.
  //! When \a regular is \e true, it is assumed that the parameter value array
  //! \a gpar forms a regular tensor-product point grid of dimension
  //! \a gpar[0].size() \a X \a gpar[1].size() \a X \a gpar[2].size().
  //! Otherwise, we assume that it contains the \a u, \a v and \a w parameters
  //! directly for each sampling point.
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integrand,
			    const RealArray* gpar, bool regular = true) const;

  //! \brief Extracts nodal results for this patch from the global vector.
  //! \param[in] globVec Global solution vector in DOF-order
  //! \param[out] nodeVec Nodal result vector for this patch
  //! \param[in] basis Which basis (or 0 for both) to extract nodal values for
  virtual void extractNodeVec(const Vector& globVec, Vector& nodeVec,
			      unsigned char = 0, int basis = 0) const;

  //! \brief Injects nodal results for this patch into a global vector.
  //! \param[in] nodeVec Nodal result vector for this patch
  //! \param[out] globVec Global solution vector in DOF-order
  //! \param[in] basis Which basis (or 0 for both) to extract nodal values for
  virtual bool injectNodeVec(const Vector& nodeVec, Vector& globVec,
                             unsigned char = 0, int basis = 0) const;

  //! \brief Generates element groups for multi-threading of interior integrals.
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] silence If \e true, suppress threading group outprint
  //! \param[in] ignoreGlobalLM If \e true ignore global multipliers in sanity check
  virtual void generateThreadGroups(const Integrand& integrand, bool silence,
                                    bool ignoreGlobalLM);
  //! \brief Generates element groups for multi-threading of boundary integrals.
  //! \param[in] lIndex Local index [1,6] of the boundary face
  //! \param[in] silence If \e true, suppress threading group outprint
  virtual void generateThreadGroups(char lIndex, bool silence, bool);

  //! \brief Returns the number of nodal points in each parameter direction.
  //! \param[out] n1 Number of nodes in first (u) direction
  //! \param[out] n2 Number of nodes in second (v) direction
  //! \param[out] n3 Number of nodes in third (w) direction
  //! \param[in] basis Which basis to return size parameters for
  virtual bool getSize(int& n1, int& n2, int& n3, int basis = 0) const;

  //! \brief Manually set the bases for special purposes (e.g. subdivision).
  //! \param[in] inp The custom constructed basis
  void setBases(std::vector<std::shared_ptr<Go::SplineVolume>>& inp) { m_basis = inp; }
protected:
  //! \brief Returns the volume in the parameter space for an element.
  //! \param[in] iel Element index
  double getParametricVolume(int iel) const;
  //! \brief Returns boundary face area in the parameter space for an element.
  //! \param[in] iel Element index
  //! \param[in] dir Local face index of the boundary face
  double getParametricArea(int iel, int dir) const;

  //! \brief Finds the global (or patch-local) node numbers on a patch boundary.
  //! \param[in] lIndex Local index of the boundary edge
  //! \param nodes Array of node numbers
  //! \param basis Which basis to grab nodes for (0 for all)
  //! \param thick Thickness of connection
  //! \param local If \e true return patch-local node numbers
  virtual void getBoundaryNodes(int lIndex, IntVec& nodes, int basis = 0,
                                int thick = 1, bool local = false) const;

private:
  std::vector<std::shared_ptr<Go::SplineVolume>> m_basis; //!< Vector of bases
};

#endif
