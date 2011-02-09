//==============================================================================
//!
//! \file AnaSol.h
//!
//! \date Des 7 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Analytical solution fields (primary and secondary)
//!
//==============================================================================

#include "AnaSol.h"

AnaSol::AnaSol()
{
  scalarSol = vectorSol = false;
}


AnaSol::AnaSol(RealFunc* sp, VecFunc* ss, VecFunc* vp, TensorFunc* vs)
  : scalSol(sp), scalSecSol(ss), vecSol(vp), vecSecSol(vs) 
{
  scalarSol = vectorSol = false;

  if (scalSol || scalSecSol) scalarSol = true;
  if (vecSol  || vecSecSol)  vectorSol = true;
}


AnaSol::~AnaSol()
{
  if (scalSol)    delete scalSol;
  if (scalSecSol) delete scalSecSol;
  if (vecSol)     delete vecSol;
  if (vecSecSol)  delete vecSecSol;

  scalarSol = vectorSol = false;
}
