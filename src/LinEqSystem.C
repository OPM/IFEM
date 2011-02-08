// $Id: LinEqSystem.C,v 1.1 2009-06-29 09:55:17 kmo Exp $
//==============================================================================
//!
//! \file LinEqSystem.C
//!
//! \date Apr 17 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Storage of a standard linear system of equations for structural FEM.
//!
//==============================================================================

#include "LinEqSystem.h"


void LinEqSystem::clear()
{
  RHS.clear();

  if (K) delete K;
  if (M) delete M;

  K = 0;
  M = 0;
}
