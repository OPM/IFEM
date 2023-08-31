// $Id$
//==============================================================================
//!
//! \file ASMs1DLag.h
//!
//! \date Apr 20 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for assembly of 1D Lagrange FE models.
//!
//==============================================================================

#ifndef _ASM_S1D_LAG_H
#define _ASM_S1D_LAG_H

#include "ASMs1D.h"
#include "Vec3.h"


/*!
  \brief Driver for assembly of 1D Lagrange FE models.
  \details This class contains methods for 1D Lagrange patches.
*/

class ASMs1DLag : public ASMs1D
{
public:
  //! \brief Default constructor.
  ASMs1DLag(unsigned char n_s = 1, unsigned char n_f = 1);
  //! \brief Special copy constructor for sharing of FE data.
  ASMs1DLag(const ASMs1DLag& patch, unsigned char n_f);
  //! \brief Default copy constructor copying everything.
  ASMs1DLag(const ASMs1DLag& patch);
  //! \brief Empty destructor.
  virtual ~ASMs1DLag() {}


  // Methods for model generation
  // ============================

  //! \brief Writes the FEM basis to given stream.
  virtual bool write(std::ostream& os, int = 0) const;

  //! \brief Generates a beam finite element model for the patch.
  //! \param[in] Zaxis Vector defining a point in the local XZ-plane
  virtual bool generateOrientedFEModel(const Vec3& Zaxis);

  //! \brief Clears the contents of the patch, making it empty.
  //! \param[in] retainGeometry If \e true, the spline geometry is not cleared.
  //! This is used to reinitialize the patch after it has been refined.
  virtual void clear(bool retainGeometry = false);

  //! \brief Returns a matrix with nodal coordinates for an element.
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in one element
  //! \param[in] iel Element index
  virtual bool getElementCoordinates(Matrix& X, int iel, bool = true) const;

  //! \brief Returns a matrix with all nodal coordinates within the patch.
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in the patch
  virtual void getNodalCoordinates(Matrix& X) const;

  //! \brief Returns the global coordinates for the given node.
  //! \param[in] inod 1-based node index local to current patch
  virtual Vec3 getCoord(size_t inod) const;

  //! \brief Updates the nodal coordinates for this patch.
  //! \param[in] displ Incremental displacements to update the coordinates with
  virtual bool updateCoords(const Vector& displ);


  // Methods for integration of finite element quantities.
  // These are the main computational methods of the ASM class hierarchy.
  // ====================================================================

  //! \brief Evaluates an integral over the interior patch domain.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrate(Integrand& integrand,
                         GlobalIntegral& glbInt, const TimeDomain& time);

  //! \brief Evaluates a boundary integral over a patch end.
  //! \param integrand Object with problem-specific data and methods
  //! \param[in] lIndex Local index of the end point
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrate(Integrand& integrand, int lIndex,
                         GlobalIntegral& glbInt, const TimeDomain& time);


  // Post-processing methods
  // =======================

  //! \brief Evaluates the geometry at a specified point.
  //! \param[in] xi Dimensionless parameters in range [0.0,1.0] of the point
  //! \param[out] param The parameters of the point in the knot-span domain
  //! \param[out] X The Cartesian coordinates of the point
  //! \return Local node number within the patch that is closest to the point
  virtual int evalPoint(const double* xi, double* param, Vec3& X) const;

  //! \brief Creates a line element model of this patch for visualization.
  //! \param[out] grid The generated line grid
  //! \param[in] npe Number of visualization nodes over each knot span
  //! \note The number of element nodes must be set in \a grid on input.
  virtual bool tesselate(ElementBlock& grid, const int* npe) const;

  using ASMs1D::evalSolution;
  //! \brief Evaluates the primary solution field at all visualization points.
  //! \details The number of visualization points is the same as the order of
  //! the Lagrange elements by default.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector in DOF-order
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const int*, int = 0) const;

  //! \brief Evaluates the primary solution field at the nodal points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const RealArray*, bool = false,
                            int = 0, int = 0) const;

  //! \brief Evaluates the secondary solution field at all visualization points.
  //! \details The number of visualization points is the same as the order of
  //! the Lagrange elements by default.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integrand,
                            const int*, char = 0) const;

  //! \brief Evaluates the secondary solution field at the nodal points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integrand,
                            const RealArray*, bool = false) const;

  using ASMs1D::getSize;
  //! \brief Returns the number of nodal points in the patch.
  virtual int getSize(int = 0) const { return nx; }

protected:
  size_t nx; //!< Number of nodes
  int    p1; //!< Polynomial order of the basis

private:
  const Vec3Vec& coord; //!< Nodal coordinates

protected:
  Vec3Vec myCoord; //!< The actual nodal coordinates
};

#endif
