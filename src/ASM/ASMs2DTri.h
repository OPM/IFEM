// $Id$
//==============================================================================
//!
//! \file ASMs2DTri.h
//!
//! \date Feb 07 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of structured 2D triangle-based FE models.
//!
//==============================================================================

#ifndef _ASM_S2D_TRI_H
#define _ASM_S2D_TRI_H

#include "ASMs2DLag.h"


/*!
  \brief Driver for assembly of structured 2D triangle-based FE models.
  \details This class contains methods for structured 2D triangular patches.
*/

class ASMs2DTri : public ASMs2DLag
{
public:
  //! \brief Default constructor.
  ASMs2DTri(unsigned char n = 2, unsigned char n_f = 2) : ASMs2DLag(n,n_f) {}
  //! \brief Special copy constructor for sharing of FE data.
  ASMs2DTri(const ASMs2DTri& patch, unsigned char n_f) : ASMs2DLag(patch,n_f) {}
  //! \brief Default copy constructor copying everything.
  ASMs2DTri(const ASMs2DTri& patch) : ASMs2DLag(patch) {}
  //! \brief Empty destructor.
  virtual ~ASMs2DTri() {}


  // Methods for model generation
  // ============================

  //! \brief Writes the FEM basis to given stream.
  virtual bool write(std::ostream& os, int = 0) const;

  //! \brief Generates the finite element topology data for the patch.
  //! \details The data generated are the element-to-node connectivity array,
  //! the nodal coordinate array, as well as global node and element numbers.
  virtual bool generateFEMTopology();

  //! \brief Adds extraordinary elements associated with a patch boundary.
  //! \param[in] dim Dimension of the boundary (should be 1)
  //! \param[in] item Local index of the boundary edge
  //! \param[in] nXn Number of extraordinary nodes
  //! \param[out] nodes Global numbers assigned to the extraordinary nodes
  virtual bool addXElms(short int dim, short int item,
                        size_t nXn, IntVec& nodes);

  //! \brief Computes the number of boundary integration points in this patch.
  virtual void getNoBouPoints(size_t& nPt, char ldim, char lindx);


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

  //! \brief Creates a triangle element model of this patch for visualization.
  //! \param[out] grid The generated triangular grid
  //! \param[in] npe Number of visualization nodes over each knot span
  virtual bool tesselate(ElementBlock& grid, const int* npe) const;

  using ASMs2DLag::evalSolution;
  //! \brief Evaluates the secondary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integrand,
                            const int*, char = 0) const;

  //! \brief Evaluates the secondary solution field at the given points.
  virtual bool evalSolution(Matrix& , const IntegrandBase&,
                            const RealArray*, bool = false) const;

  using ASMs2DLag::generateThreadGroups;
  //! \brief Generates element groups for multi-threading of interior integrals.
  virtual void generateThreadGroups(const Integrand&, bool, bool);
};

#endif
