// $Id: LinEqSystem.h,v 1.4 2010-05-06 17:31:18 kmo Exp $
//==============================================================================
//!
//! \file LinEqSystem.h
//!
//! \date Apr 17 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Storage of a standard linear system of equations for structural FEM.
//!
//==============================================================================

#ifndef _LIN_EQ_SYSTEM_H
#define _LIN_EQ_SYSTEM_H

#include "SystemMatrix.h"


/*!
  \brief Class for storage of a linear system of equations.
*/

class LinEqSystem
{
public:
  SystemMatrix* K;   //!< System stiffness matrix
  SystemMatrix* M;   //!< System mass matrix
  StdVector     RHS; //!< System right-hand-side vector

  //! \brief Default constructor.
  LinEqSystem() { K = 0; M = 0; }

  //! \brief The destructor frees the dynamically allocated objects.
  ~LinEqSystem() { if (K) delete K; if (M) delete M; }

  //! \brief Erases the system matrices and frees dynamically allocated storage.
  void clear();
};

#endif
