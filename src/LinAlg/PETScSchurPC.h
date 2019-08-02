// $Id$
//==============================================================================
//!
//! \file PETScSchurPC.h
//!
//! \date Jun 4 2019
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Schur-complement preconditioner using PETSc.
//!
//==============================================================================

#include "PETScMatrix.h"
#include <memory>


/*!
  \brief Class implementing a Schur-complement preconditioner for 2-block systems.
 */


class PETScSchurPC {
public:
  //! \brief The constructor sets up the preconditioner.
  //! \param pc_init The PETSc PC to set up
  //! \param blocks The matrix blocks
  //! \param params Linear solver parameters
  //! \param adm Processes administrator
  PETScSchurPC(PC& pc_init, const std::vector<Mat>& blocks,
               const LinSolParams::BlockParams& params, const ProcessAdm& adm);

  //! \brief The destructor frees the PETSc structures.
  ~PETScSchurPC();

  //! \brief PETSc compatible function applying the Schur complement matrix.
  //! \param A Shell matrix to apply
  //! \param x Vector to apply matrix to
  //! \param y Result of matrix-vector product
  static PetscErrorCode Apply_Schur(Mat A, Vec x, Vec y);

  //! \brief PETSc compatible function applying the Schur complement preconditioner.
  //! \param pc Shell preconditioner to apply
  //! \param x Vector to apply preconditioner to
  //! \param y Result of preconditioner evaluation
  static PetscErrorCode Apply_Outer(PC pc, Vec x, Vec y);

  //! \brief Destroy a Schur preconditioner.
  //! \param pc Preconditioner to destroy
  static PetscErrorCode Destroy(PC pc);

protected:
  KSP inner_ksp; //!< The KSP for the approximation of inner matrix inverse
  KSP outer_ksp; //!< The KSP for the approximated Schur complement
  Mat outer_mat; //!< Matrix shell describing the Schur complement
  Vec tmp; //!< Temporary vector
  const std::vector<Mat>* m_blocks; //!< Matrix blocks in system
};
