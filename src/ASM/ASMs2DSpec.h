// $Id$
//==============================================================================
//!
//! \file ASMs2DSpec.h
//!
//! \date Mar 22 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for assembly of structured 2D Spectral FE models.
//!
//==============================================================================

#ifndef _ASM_S2D_SPEC_H
#define _ASM_S2D_SPEC_H

#include "ASMs2DLag.h"


/*!
  \brief Driver for assembly of structured 2D Spectral FE models.
  \details This class contains methods for structured 2D Spectral patches.
*/

class ASMs2DSpec : public ASMs2DLag
{
public:
  //! \brief Default constructor.
  ASMs2DSpec(unsigned char n_s = 2, unsigned char n_f = 2)
    : ASMs2DLag(n_s,n_f) {}
  //! \brief Copy constructor.
  ASMs2DSpec(const ASMs2DSpec& patch, unsigned char n_f = 0)
    : ASMs2DLag(patch,n_f) {}
  //! \brief Empty destructor.
  virtual ~ASMs2DSpec() {}


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

  using ASMs2DLag::evalSolution;
  //! \brief Evaluates the secondary solution field at the given points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] gpar Parameter values of the result sampling points
  //! \param[in] regular Flag indicating how the sampling points are defined
  virtual bool evalSolution(Matrix& sField, const IntegrandBase& integrand,
			    const RealArray* gpar, bool regular = true) const;

protected:

  // Internal utility methods
  // ========================

  //! \brief Calculates parameter values for all visualization nodal points.
  //! \param[out] prm Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1)
  //! \param[in] nSegSpan Number of visualization segments over each knot-span
  //!
  //! This method is reimplemented to return the parameter values of the
  //! Gauss-Lobatto-Legendre points, which define the FE nodes when using
  //! spectral elements. This method is therefore also used when creating the
  //! FE model (the visualization nodes and the FE nodes are identical).
  virtual bool getGridParameters(RealArray& prm, int dir, int nSegSpan) const;
};

#endif
