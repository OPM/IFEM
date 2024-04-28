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

class ElmMats;


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

  //! \brief Parses structural damping parameters from an XML element.
  bool parseParams(const tinyxml2::XMLElement* elem);

  //! \brief Calculates the dynamic solution from the previous modal solution.
  //! \param[in] mSol Previous modal solution
  const Vectors& expandSolution(const Vectors& mSol);

  //! \brief Administers assembly of the modal equation system.
  //! \param[in] time Parameters for time-dependent simulations
  //! \param[in] mSol Previous modal solution
  //! \param[in] beta Newmark time integration parameter
  //! \param[in] gamma Newmark time integration parameter
  bool assembleModalSystem(const TimeDomain& time, const Vectors& mSol,
                           double beta, double gamma);

  //! \brief Swaps the modal equation system before/after load vector assembly.
  bool swapSystem(AlgEqSystem*& sys, SAM*& sam);

  //! \brief Writes the eigenmodes to a serialization container.
  bool saveModes(std::map<std::string,std::string>& data) const;
  //! \brief Restores the eigenmodes from a serialization container.
  bool restoreModes(const std::map<std::string,std::string>& data);

public:
  //! \brief Expands and returns the current dynamic solution.
  virtual const Vectors& expandSolution(const Vectors&, bool = false) = 0;
  //! \brief Returns the current expanded dynamic solution.
  //! \param[in] idx Solution vector index
  const Vector& expandedSolution(int idx) const;
  //! \brief Returns the number of expanded dynamic solution vectors.
  size_t numExpSolution() const { return sol.size(); }

  //! \brief Projects the secondary solution associated with the eigenmodes.
  virtual bool projectModes(Matrices&, std::vector<std::string>&,
                            SIMoptions::ProjectionMethod) { return false; }

protected:
  std::vector<Mode>& myModes; //!< Array of eigenmodes

  Vector  Rhs; //!< Current right-hand-side load vector of the dynamic system
  Vectors sol; //!< Expanded solution vectors from the modal solution

  bool   parsed; //!< Set to \e true after the model has been initialized
  double alpha1; //!< Mass-proportional damping parameter
  double alpha2; //!< Stiffness-proportional damping parameter

private:
  AlgEqSystem* modalSys; //!< The modal equation system
  SAM*         modalSam; //!< Auxiliary data for FE assembly management
  ElmMats*     myElmMat; //!< Nodal (single-DOF) element matrices
};

#endif
