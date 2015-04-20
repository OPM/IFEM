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

#include "GoTools/geometry/SplineCurve.h"

#include "ASMs1DC1.h"


bool ASMs1DC1::generateFEMTopology ()
{
  if (curv->order() > 2)
    return this->ASMs1D::generateFEMTopology();

  std::cerr <<" *** ASMs1DC1::generateFEMTopology: The polynomial order ("
            << curv->order() <<") is too low.\n    "
            <<" C1-continuity requires at least quadratic order."<< std::endl;
  return false;
}


void ASMs1DC1::constrainNode (double xi, int dof, int code, char basis)
{
  if (xi < 0.0 || xi > 1.0) return;

  int n1 = this->getSize(basis);
  int n=0;
  for (char i = 1; i < basis; ++i)
    n += this->getSize(i);

  int node = n+1;
  if (xi > 0.0) node += int(0.5+(n1-1)*xi);

  this->prescribe(node,dof%1000,code);
  if (dof < 1000 || code > 0) return;

  // Also fix the (up to two) neighboring nodes
  if (node > 1)  this->prescribe(node-1,dof%1000,0);
  if (node < n1) this->prescribe(node+1,dof%1000,0);
}
