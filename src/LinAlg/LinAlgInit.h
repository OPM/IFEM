// $Id: LinAlgInit.h,v 1.2 2011-01-02 15:50:35 kmo Exp $
//==============================================================================
//!
//! \file LinAlgInit.h
//!
//! \date May 06 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Initializer for linear algebra packages.
//!
//==============================================================================

#ifndef _LINALG_INIT_H
#define _LINALG_INIT_H


/*!
  \brief Initializer for linear algebra packages.
*/

class LinAlgInit
{
public:
  //! \brief The constructor uses the command-line arguments to set up things.
  LinAlgInit(int argc, char** argv);
  //! \brief The destructor finalizes the linear algebra packages.
  ~LinAlgInit();

  int myPid; //!< Processor ID in parallel simulations
};

#endif
