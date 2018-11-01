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
#include "Utilities.h"


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
  if (node < 1 || dof < 1000) return node;

  // Also fix the (up to two) neighboring nodes
  for (int neighbour = node-1; neighbour <= node+1; neighbour += 2)
    if (neighbour > 0 && neighbour <= this->getSize(basis))
    {
      if (dof%1000 && code == 0)
        // The node is clamped, fix the neighboring node
        this->prescribe(neighbour,dof/1000,code);
      else for (int ldof : utl::getDigits(dof/1000))
        // The node has a prescribed rotation, add an MPC for that
        this->add2PC(MLGN[neighbour-1],ldof,MLGN[node-1],code);
    }

  return node;
}
