// $Id: ASMs1DLag.h,v 1.2 2010-09-07 07:51:46 kmo Exp $
//==============================================================================
//!
//! \file ASMs1DLag.h
//!
//! \date Apr 20 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for assembly of structured 1D Lagrange FE models.
//!
//==============================================================================

#ifndef _ASM_S1D_LAG_H
#define _ASM_S1D_LAG_H

#include "ASMs1D.h"
#include "Vec3.h"


/*!
  \brief Driver for assembly of structured 2D Lagrange FE models.
  \details This class contains methods for structured 2D Lagrange patches.
*/

class ASMs1DLag : public ASMs1D
{
public:
  //! \brief Constructor creating an instance by reading the given file.
  ASMs1DLag(const char* fileName, unsigned char n_s = 1, unsigned char n_f = 1);
  //! \brief Constructor creating an instance by reading the given input stream.
  ASMs1DLag(std::istream& is, unsigned char n_s = 1, unsigned char n_f = 1);
  //! \brief Default constructor creating an empty patch.
  ASMs1DLag(unsigned char n_s = 1, unsigned char n_f = 1) : ASMs1D(n_s,n_f) {}
  //! \brief Empty destructor.
  virtual ~ASMs1DLag() {}


  // Methods for model generation
  // ============================

  //! \brief Generates the finite element topology data for the patch.
  //! \details The data generated are the element-to-node connectivity array,
  //! and the global node and element numbers.
  virtual bool generateFEMTopology();

  //! \brief Clears the contents of the patch, making it empty.
  virtual void clear();

  //! \brief Returns the global coordinates for the given node.
  //! \param[in] inod 1-based node index local to current patch
  virtual Vec3 getCoord(size_t inod) const;


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

  //! \brief Evaluates a boundary integral over a patch end.
  //! \param integrand Object with problem-specific data and methods
  //! \param[in] lIndex Local index of the end point
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param locInt Vector of element-wise contributions to \a glbInt
  virtual bool integrate(Integrand& integrand, int lIndex,
			 GlobalIntegral& glbInt, const TimeDomain& time,
			 const LintegralVec& locInt = LintegralVec());


  // Post-processing methods
  // =======================

  //! \brief Creates a quad element model of this patch for visualization.
  //! \param[out] grid The generated quadrilateral grid
  //! \param[in] npe Number of visualization nodes over each knot span
  //! \note The number of element nodes must be set in \a grid on input.
  virtual bool tesselate(ElementBlock& grid, const int* npe) const;

  //! \brief Evaluates the primary solution field at all visualization points.
  //! \details The number of visualization points is the same as the order of
  //! the Lagrange elements by default.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector in DOF-order
  //! \param[in] npe Number of visualization nodes over each knot span
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
			    const int* npe) const;

  //! \brief Evaluates the secondary solution field at all visualization points.
  //! \details The number of visualization points is the same as the order of
  //! the Lagrange elements by default.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] npe Number of visualization nodes over each knot span
  virtual bool evalSolution(Matrix& sField, const Integrand& integrand,
			    const int* npe) const;

protected:

  // Internal utility methods
  // ========================

  //! \brief Returns a matrix with nodal coordinates for an element.
  //! \param[in] iel Element index
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in one element
  virtual bool getElementCoordinates(Matrix& X, int iel) const;
  //! \brief Returns a matrix with all nodal coordinates within the patch.
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in the patch
  virtual void getNodalCoordinates(Matrix& X) const;

  //! \brief Returns the number of nodal points in the patch.
  virtual int getSize() const { return nx; }

private:
  //! \brief The number of nodes in each direction for the patch
  size_t            nx;    //!< Number of nodes in first parameter direction
  std::vector<Vec3> coord; //!< Nodal coordinates
};

#endif
