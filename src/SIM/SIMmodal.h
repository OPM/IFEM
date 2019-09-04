// $Id$
//==============================================================================
//!
//! \file SIMmodal.h
//!
//! \date Aug 30 2019
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of modal FEM system.
//!
//==============================================================================

#ifndef _SIM_MODAL_H
#define _SIM_MODAL_H

#include "SIMbase.h"

class NewmarkMats;


/*!
  \brief Class with support for assembly of modal linear equation systems.
  \details This class contains a separate AlgEqSystem object and an associated
  SAM object, used for assembling the modal system of equations.
*/

class SIMmodal
{
protected:
  //! \brief The constructor initializes the reference to the eigenmodes.
  SIMmodal(std::vector<Mode>& modes);
  //! \brief The destructor deletes the dynamically allocated members.
  virtual ~SIMmodal();

  //! \brief Calculates the dynamic solution from the previous modal solution.
  //! \param[in] mSol Previous modal solution
  //! \param[out] pSol Dynamic solution vectors
  bool expandSolution(const Vectors& mSol, Vectors& pSol) const;

  //! \brief Administers assembly of the modal equation system.
  //! \param[in] time Parameters for time-dependent simulations
  //! \param[in] mSol Previous modal solution
  //! \param[in] Rhs Current right-hand-side load vector
  //! \param[in] beta Newmark time integration parameter
  //! \param[in] gamma Newmark time integration parameter
  bool assembleModalSystem(const TimeDomain& time,
                           const Vectors& mSol, const Vector& Rhs,
                           double beta, double gamma);

  //! \brief Swaps the modal equation system before/after load vector assembly.
  bool swapSystem(AlgEqSystem*& sys, SAM*& sam);

  std::vector<Mode>& myModes; //!< Array of eigenmodes

private:
  AlgEqSystem* modalSys; //!< The modal equation system
  SAM*         modalSam; //!< Auxiliary data for FE assembly management
  NewmarkMats* myElmMat; //!< Nodal (single-DOF) element matrices
};

#endif
