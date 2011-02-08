// $Id: PressureLoad.h,v 1.2 2009-05-05 09:25:09 kmo Exp $
//==============================================================================
//!
//! \file PressureLoad.h
//!
//! \date Jan 27 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of a constant pressure load.
//!
//==============================================================================

#ifndef _PRESSURELOAD_H
#define _PRESSURELOAD_H


/*!
  \brief Struct for representing a constant pressure load.
*/

struct PLoad
{
  int volp; //!< Volume spline patch index [0,nvol]
  int face; //!< Local face index on spline volume [-3,3]
  int pdir; //!< Pressure direction [0,3] (0=normal pressure)
  double p; //!< Actual pressure value

  //! \brief Constructor creating an initialized pressure instance.
  //! \param[in] s Patch index
  //! \param[in] f Local face index
  //! \param[in] d Pressure direction
  //! \param[in] v Pressure value
  PLoad(int s, int f, int d, double v) : volp(s), face(f), pdir(d), p(v) {}
};

#endif
