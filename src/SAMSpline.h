// $Id: SAMSpline.h,v 1.2 2009-05-27 12:52:34 kmo Exp $
//==============================================================================
//!
//! \file SAMSpline.h
//!
//! \date Dec 10 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of FE matrices into system matrices for SplineVolume models.
//!
//==============================================================================

#ifndef _SAM_SPLINE_H
#define _SAM_SPLINE_H

#include "SAM.h"

class VolumePatch;
class LinEqSystem;


/*!
  \brief This is a sub-class of SAM with added functionality for spline models.
  \details It contains some additional functions for initializing the SAM-data
  for an FE model based on SplineVolume patches.
*/

class SAMSpline : public SAM
{
public:
  //! \brief Default constructor which initializes an empty SAMSpline object.
  SAMSpline() : SAM() {}
  //! \brief Empty destructor.
  virtual ~SAMSpline() {}

  //! \brief Allocates the dynamic arrays and populates them with data
  //! for a single-patch spline model.
  bool init(const VolumePatch& model);
  //! \brief Allocates the dynamic arrays and populates them with data
  //! for a multi-patch spline model.
  bool init(const std::vector<VolumePatch*>& model, int numNod = 0);

  //! \brief Performs the necessary initialization of the system matrices
  //! prior to the element assembly.
  //! \details This method must be called once before the first call to
  //! \a assembleSystem for a given load case or time step.
  //! \param sys The system left-hand-side matrices and right-hand-side vector
  //! \param withRHS If \a true, initialize the right-hand-side vector too
  //! \return \e false if no free DOFs in the system, otherwise \e true
  bool initForAssembly(LinEqSystem& sys, bool withRHS = false) const;

protected:
  //! \brief Initializes the nodal arrays \a MINEX, \a MADOF and \a MSC.
  void initNodeDofs(const std::vector<VolumePatch*>& p);
  //! \brief Initializes the element topology arrays \a MPMNPC and \a MMNPC.
  void initElementConn(const std::vector<VolumePatch*>& p);
  //! \brief Initializes the multi-point constraint arrays
  //! \a MPMCEQ, \a MMCEQ and \a TTCC.
  void initConstraintEqs(const std::vector<VolumePatch*>& p);
};

#endif
