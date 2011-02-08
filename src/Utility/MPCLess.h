// $Id: MPCLess.h,v 1.1 2009-12-03 15:37:07 kmo Exp $
//==============================================================================
//!
//! \file MPCLess.h
//!
//! \date Oct 1 2007
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Comparison of two \a MPC objects in a \a set<MPC*,MPCLess> object.
//!
//==============================================================================

#ifndef _MPC_LESS_H
#define _MPC_LESS_H

class MPC;


/*!
  \brief Functor used to sort a set of MPC pointers
*/

class MPCLess
{
public:
  //! \brief Comparison operator used when inserting an MPC-pointer into a
  //! \a set<MPC*,MPCLess> object.
  bool operator()(const MPC* lhs, const MPC* rhs) const;
  //! \brief Indicates whether only the slave dof number should affect sorting.
  //! \details The default is to also compare the associated coefficients.
  static bool compareSlaveDofOnly;
};

#endif
