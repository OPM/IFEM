// $Id: ASMs3Dmx.h,v 1.2 2010-12-30 15:02:02 kmo Exp $
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
  //! \brief Constructor creating an instance by reading the given file.
  ASMs3Dmx(const char* fileName, bool checkRHS = false,
	   char n_f1 = 3, unsigned char n_f2 = 1);
  //! \brief Constructor creating an instance by reading the given input stream.
  ASMs3Dmx(std::istream& is, bool checkRHS = false,
	   char n_f1 = 3, unsigned char n_f2 = 1);
  //! \brief Default constructor creating an empty patch.
  ASMs3Dmx(bool checkRHS = false, char n_f1 = 3, unsigned char n_f2 = 1);
  //! \brief Empty destructor.
  virtual ~ASMs3Dmx() {}


  // Methods for model generation
  // ============================

  //! \brief Generates the finite element topology data for the patch.
  //! \details The data generated are the element-to-node connectivity array,
  //! the node-to-IJ-index array, as well as global node and element numbers.
  virtual bool generateFEMTopology();

  //! \brief Clears the contents of the patch, making it empty.
  virtual void clear();

  //! \brief Returns the global coordinates for the given node.
  //! \param[in] inod 1-based node index local to current patch
  virtual Vec3 getCoord(size_t inod) const;

  //! \brief Returns the number of solution fields.
  virtual unsigned char getNoFields(int basis = 0) const;
  //! \brief Returns the number of DOFs per node.
  //! \param[in] inod 1-based node index local to current patch
  virtual unsigned char getNodalDOFs(size_t inod) const;

  //! \brief Initializes the patch level MADOF array for mixed problems.
  virtual void initMADOF(const int* sysMadof);


  // Methods for integration of finite element quantities.
  // These are the main computational methods of the ASM class hierarchy.
  // ====================================================================

  //! \brief Evaluates an integral over the interior patch domain.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param locInt Vector of element-wise contributions to \a glbInt
  virtual bool integrate(Integrand& integrand,
			 GlobalIntegral& glbInt, const TimeDomain& time,
			 const LintegralVec& locInt = LintegralVec());

  //! \brief Evaluates a boundary integral over a patch edge.
  //! \param integrand Object with problem-specific data and methods
  //! \param[in] lIndex Local index of the boundary edge
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param locInt Vector of element-wise contributions to \a glbInt
  virtual bool integrate(Integrand& integrand, int lIndex,
			 GlobalIntegral& glbInt, const TimeDomain& time,
			 const LintegralVec& locInt = LintegralVec());


  // Post-processing methods
  // =======================

  //! \brief Evaluates the primary solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector in DOF-order
  //! \param[in] npe Number of visualization nodes over each knot span
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
			    const int* npe) const;

  //! \brief Evaluates the secondary solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] npe Number of visualization nodes over each knot span
  virtual bool evalSolution(Matrix& sField, const Integrand& integrand,
			    const int* npe) const;

  //! \brief Extracts nodal results for this patch from the global vector.
  //! \param[in] globVec Global solution vector in DOF-order
  //! \param[out] nodeVec Nodal result vector for this patch
  virtual void extractNodeVec(const Vector& globVec, Vector& nodeVec,
			      unsigned char = 0) const;

protected:

  // Internal utility methods
  // ========================

  //! \brief Returns a matrix with nodal coordinates for an element.
  //! \param[in] iel Element index
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in one element
  virtual bool getElementCoordinates(Matrix& X, int iel) const;

  //! \brief Returns the number of nodal points in each parameter direction.
  //! \param[out] n1 Number of nodes in first (u) direction
  //! \param[out] n2 Number of nodes in second (v) direction
  //! \param[out] n3 Number of nodes in third (w) direction
  //! \param[in] basis Which basis to return size parameters for
  virtual bool getSize(int& n1, int& n2, int& n3, int basis = 0) const;

private:
  Go::SplineVolume* basis1; //!< Pointer to spline object for the first basis
  Go::SplineVolume* basis2; //!< Pointer to spline object for the second basis
};

#endif
