// $Id$
//==============================================================================
//!
//! \file ASMs1DSpec.h
//!
//! \date Apr 20 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for assembly of structured 1D Spectral FE models.
//!
//==============================================================================

#ifndef _ASM_S1D_SPEC_H
#define _ASM_S1D_SPEC_H

#include "ASMs1DLag.h"


/*!
  \brief Driver for assembly of structured 1D Spectral FE models.
  \details This class contains methods for structured 1D Spectral patches.
*/

class ASMs1DSpec : public ASMs1DLag
{
public:
  //! \brief Constructor creating an instance by reading the given file.
  ASMs1DSpec(const char* fNam = 0, unsigned char n_s = 1, unsigned char n_f = 1)
    : ASMs1DLag(fNam,n_s,n_f) {}
  //! \brief Constructor creating an instance by reading the given input stream.
  ASMs1DSpec(std::istream& is, unsigned char n_s = 1, unsigned char n_f = 1)
    : ASMs1DLag(is,n_s,n_f) {}
  //! \brief Empty destructor.
  virtual ~ASMs1DSpec() {}


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

  //! \brief Evaluates the secondary solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] integrand Object with problem-specific data and methods
  //! \param[in] npe Number of visualization nodes over each knot span
  virtual bool evalSolution(Matrix& sField, const Integrand& integrand,
			    const int* npe) const;

protected:

  // Internal utility methods
  // ========================

  //! \brief Calculates parameter values for all visualization nodal points.
  //! \param[out] prm Parameter values for all points
  //! \param[in] nSegSpan Number of visualization segments over each knot-span
  //!
  //! This method is reimplemented to return the parameter values of the
  //! Gauss-Lobatto-Legendre points, which define the FE nodes when using
  //! spectral elements. This method is therefore also used when creating the
  //! FE model (the visualization nodes and the FE nodes are identical).
  virtual bool getGridParameters(RealArray& prm, int nSegSpan) const;
};

#endif
