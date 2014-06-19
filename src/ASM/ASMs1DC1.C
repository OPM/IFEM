// $Id$
//==============================================================================
//!
//! \file ASMs1DC1.C
//!
//! \date May 15 2014
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of C1-continuous structured 1D spline FE models.
//!
//==============================================================================

#include "ASMs1DC1.h"


void ASMs1DC1::constrainNode (double xi, int dof, int code)
{
  if (xi < 0.0 || xi > 1.0) return;

  int n1 = this->getSize();

  int node = 1;
  if (xi > 0.0) node += int(0.5+(n1-1)*xi);

  this->prescribe(node,dof%1000,code);
  if (dof < 1000 || code > 0) return;

  // Also fix the (up to two) neighboring nodes
  if (node > 1)  this->prescribe(node-1,dof%1000,0);
  if (node < n1) this->prescribe(node+1,dof%1000,0);
}
