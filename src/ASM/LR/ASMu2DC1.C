// $Id$
//==============================================================================
//!
//! \file ASMu2DC1.C
//!
//! \date Oct 5 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of C1-continuous 2D LR-spline FE models.
//!
//==============================================================================

#include "LRSpline/LRSplineSurface.h"

#include "ASMu2DC1.h"


bool ASMu2DC1::generateFEMTopology ()
{
  if (!(geomB = this->createLRfromTensor()))
    return false;
  else if (geomB->order(0) > 2 || geomB->order(1) > 2)
    return this->ASMu2D::generateFEMTopology();

  std::cerr <<" *** ASMu2DC1::generateFEMTopology:"
            <<" The polynomial order "<< geomB->order(0) <<"x"<< geomB->order(1)
            <<" is too low.\n     C1-continuity requires"
            <<" at least quadratic order."<< std::endl;
  return false;
}


void ASMu2DC1::constrainEdge (int dir, bool open, int dof, int code, char basis)
{
  if (dof < 1000)
    this->ASMu2D::constrainEdge(dir,open,dof,code,basis);
  else if (basis > 1 || code != 0)
    std::cerr <<"  ** ASMu2DC1::constrainEdge: Not implemented for"
              <<" mixed basis and/or inhomogeneous Dirichlet conditions"
              <<" (ignored)."<< std::endl;
  else
  {
    // Figure out what edge we are at
    DirichletEdge de(lrspline.get(),dir,dof);

    // Get all basis functions on this boundary edge,
    // and the layer of functions next to the edge functions
    std::vector<LR::Basisfunction*> edgeFunc1, edgeFunc2;
    de.lr->getEdgeFunctions(edgeFunc1,de.edg,1);
    de.lr->getEdgeFunctions(edgeFunc2,de.edg,2);

    if (dof%1000) // Add constraints for all basis functions on the edge
      for (LR::Basisfunction* b : edgeFunc1)
        if (!open || !de.isCorner(1+b->getId())) // skip ends for open boundary
          this->prescribe(1+b->getId(),dof%1000,code);

    // Add constraints for the basis functions just inside the edge
    for (LR::Basisfunction* b : edgeFunc2)
      if (std::find(edgeFunc1.begin(),edgeFunc1.end(),b) == edgeFunc1.end())
        if (!open || !de.isCorner(1+b->getId())) // skip ends for open boundary
          this->prescribe(1+b->getId(),dof/1000,code);
  }
}
