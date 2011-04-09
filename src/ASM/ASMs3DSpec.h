// $Id$
//==============================================================================
//!
//! \file ASMs3DSpec.h
//!
//! \date Mar 22 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for assembly of structured 3D Spectral FE models.
//!
//==============================================================================

#ifndef _ASM_S3D_SPEC_H
#define _ASM_S3D_SPEC_H

#include "ASMs3DLag.h"


/*!
  \brief Driver for assembly of structured 3D Spectral FE models.
  \details This class contains methods for structured 3D Spectral patches.
*/

class ASMs3DSpec : public ASMs3DLag
{
public:
  //! \brief Constructor creating an instance by reading the given file.
  ASMs3DSpec(const char* fNam = 0, bool checkRHS = false, unsigned char n_f = 3)
    : ASMs3DLag(fNam,checkRHS,n_f) {}
  //! \brief Constructor creating an instance by reading the given input stream.
  ASMs3DSpec(std::istream& is, bool checkRHS = false, unsigned char n_f = 3)
    : ASMs3DLag(is,checkRHS,n_f) {}
  //! \brief Empty destructor.
  virtual ~ASMs3DSpec() {}


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

  //! \brief Evaluates a boundary integral over a patch face.
  //! \param integrand Object with problem-specific data and methods
  //! \param[in] lIndex Local index [1,6] of the boundary face
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  //! \param locInt Vector of element-wise contributions to \a glbInt
  virtual bool integrate(Integrand& integrand, int lIndex,
			 GlobalIntegral& glbInt, const TimeDomain& time,
			 const LintegralVec& locInt = LintegralVec());

  //! \brief Evaluates a boundary integral over a patch edge.
  //! \param integrand Object with problem-specific data and methods
  //! \param[in] lEdge Local index [1,12] of the patch edge
  //! \param glbInt The integrated quantity
  //! \param[in] time Parameters for nonlinear/time-dependent simulations
  virtual bool integrateEdge(Integrand& integrand, int lEdge,
  			     GlobalIntegral& glbInt, const TimeDomain& time);


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
  //! \param[out] prm Parameter values in given direction for all points
  //! \param[in] dir Parameter direction (0,1,2)
  //! \param[in] nSegSpan Number of visualization segments over each knot-span
  //!
  //! This method is reimplemented to return the parameter values of the
  //! Gauss-Lobatto-Legendre points, which define the FE nodes when using
  //! spectral elements. This method is therefore also used when creating the
  //! FE model (the visualization nodes and the FE nodes are identical).
  virtual bool getGridParameters(RealArray& prm, int dir, int nSegSpan) const;
};

#endif
