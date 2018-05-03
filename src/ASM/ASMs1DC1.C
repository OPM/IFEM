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


int ASMs1DC1::constrainNode (double xi, int dof, int code, char basis)
{
  int node = this->ASMs1D::constrainNode(xi,dof%1000,code,basis);
  if (node > 0 && dof >= 1000 && code == 0)
  {
    // Also fix the (up to two) neighboring nodes
    if (node > 1)
      this->prescribe(node-1,dof/1000,code);
    if (node < this->getSize(basis))
      this->prescribe(node+1,dof/1000,code);
  }

  return node;
}
