// $Id$
//==============================================================================
//!
//! \file ASMmxBase.h
//!
//! \date Dec 28 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for spline-based mixed finite element assembly drivers.
//!
//==============================================================================

#ifndef _ASM_MX_BASE_H
#define _ASM_MX_BASE_H

#include "MatVec.h"
#include <memory>

namespace Go {
  class SplineSurface;
  class SplineVolume;
}


/*!
  \brief Base class for spline-based mixed finite element assembly drivers.
*/

class ASMmxBase
{
protected:
  //! \brief The constructor sets the number of field variables.
  //! \param[in] n_f Number of nodal variables in each field
  explicit ASMmxBase(const std::vector<unsigned char>& n_f) : nfx(n_f) {}

  //! \brief Initializes the patch level MADOF array.
  //! \param[in] MLGN Matrix of local-to-global node numbers
  //! \param[out] sysMadof System-level matrix of accumulated DOFs
  void initMx(const std::vector<int>& MLGN, const int* sysMadof);

  //! \brief Extracts nodal results for this patch from the global vector.
  //! \param[in] glbVec Global solution vector in DOF-order
  //! \param[out] nodVec Nodal result vector for this patch
  //! \param[in] basis Which basis to extract the nodal values for
  void extractNodeVecMx(const Vector& glbVec, Vector& nodVec, int basis) const;

  //! \brief Injects nodal results for this patch into a global vector.
  //! \param[out] glbVec Global solution vector in DOF-order
  //! \param[in] nodVec Nodal result vector for this patch
  //! \param[in] basis Which basis to inject the nodal values for
  void injectNodeVecMx(Vector& glbVec, const Vector& nodVec, int basis) const;

  //! \brief Extracts the primary solution field at the specified nodes.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  //! \param[in] nodes 1-based local node numbers to extract solution for
  bool getSolutionMx(Matrix& sField, const Vector& locSol,
                     const std::vector<int>& nodes) const;

public:
  //! \brief Enum defining available mixed formulation types.
  enum MixedType {
    NONE = 0,                  //!< No Mixed formulation
    FULL_CONT_RAISE_BASIS1,    //!< Full continuity, raise order and use as basis 1
    REDUCED_CONT_RAISE_BASIS1, //!< Reduced continuity, raise order and use as basis 1
    FULL_CONT_RAISE_BASIS2,    //!< Full continuity, raise order and use as basis 2
    REDUCED_CONT_RAISE_BASIS2, //!< Reduced continuity, raise order and use as basis 2
    DIV_COMPATIBLE,            //!< Div-compatible space for incompressible problems
    SUBGRID,                   //!< Sub-grid spaces
  };

  static MixedType Type; //!< Type of mixed formulation used
  static char  geoBasis; //!< 1-based index of basis representing the geometry

protected:
  typedef std::shared_ptr<Go::SplineSurface> SurfacePtr; //!< Pointer to spline
  typedef std::shared_ptr<Go::SplineVolume>  VolumePtr;  //!< Pointer to spline
  typedef std::vector<SurfacePtr> SurfaceVec; //!< Convenience type
  typedef std::vector<VolumePtr>  VolumeVec;  //!< Convenience type

  //! \brief Establish mixed bases in 2D.
  //! \param[in] surf The base basis to use
  //! \param[in] type The type of bases to establish
  //! \return Vector with bases
  static SurfaceVec establishBases(Go::SplineSurface* surf, MixedType type);

  //! \brief Establish mixed bases in 3D.
  //! \param[in] svol The base basis to use
  //! \param[in] type The type of bases to establish
  //! \return Vector with bases
  static VolumeVec establishBases(Go::SplineVolume* svol, MixedType type);

  //! \brief Returns a C^p-1 basis of one degree higher than \a *surf.
  static Go::SplineSurface* raiseBasis(Go::SplineSurface* surf);
  //! \brief Returns a C^p-1 basis of one degree higher than \a *svol.
  static Go::SplineVolume*  raiseBasis(Go::SplineVolume* svol);

  //! \brief Mask off DOFs not residing on the specified basis
  //! \param[in] dofs Encoded DOF indices (like 123456 meaning dofs 1 to 6)
  //! \param[in] basis Which basis to return DOF indices for
  //! \return Encoded DOF indices related to the specified basis only
  int maskDOFs(int dofs, char basis) const;

private:
  std::vector<int> MADOF; //!< Matrix of accumulated DOFs for this patch

protected:
  //! Number of basis functions per element in each basis
  std::vector<size_t> elem_size;
  //! Total number of basis functions in each basis
  std::vector<size_t> nb;

  std::vector<unsigned char> nfx; //!< Number of fields on each basis
};

#endif
