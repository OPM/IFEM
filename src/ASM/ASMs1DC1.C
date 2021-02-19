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
#include "Tensor.h"
#include "Vec3.h"
#include "Vec3Oper.h"
#include "MPC.h"


bool ASMs1DC1::generateFEMTopology ()
{
  if (curv->order() > 2)
    return this->ASMs1D::generateFEMTopology();

  std::cerr <<" *** ASMs1DC1::generateFEMTopology: The polynomial order ("
            << curv->order() <<") is too low.\n    "
            <<" C1-continuity requires at least quadratic order."<< std::endl;
  return false;
}


int ASMs1DC1::constrainNode (double xi, int dof, int code)
{
  int node = this->ASMs1D::constrainNode(xi,dof%1000,code);
  if (node < 1 || dof < 1000) return node;

  // Also fix the (up to two) neighboring nodes
  for (int neighbour = node-1; neighbour <= node+1; neighbour += 2)
    if (neighbour > 0 && neighbour <= this->getSize())
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


size_t ASMs1DC1::constrainEndLocal (int dir, int dof, int code)
{
  if (shareFE == 'F')
  {
    std::cerr <<"\n *** ASMs1DC1::constrainEndLocal: Logic error, can not have"
              <<" constraints in local axes for shared patches."<< std::endl;
    return 0;
  }
  else if (dof < 1000) // No rotation constraint
    return this->ASMs1D::constrainEndLocal(dir,dof,code);
  else if (code != 0)
  {
    std::cerr <<"\n *** ASMs1DC1::constrainEndLocal: Prescribed DOFs in local"
              <<" axes not implemented."<< std::endl;
    return 0;
  }

  // Find the local-to-global transformation matrix at this end point,
  // where the X-axis corresponds to the tangent direction of the curve
  double uEnd = dir > 0 ? curv->endparam() : curv->startparam();
  Tensor Tlg(this->getLocal2Global(uEnd));

  // We need extra nodes representing the local (master) DOFs at this point
  int added = 0;
  int iSnod = dir > 0 ? this->getSize()-1 : 0;
  for (int i = 0; i < 2; i++, iSnod += (dir > 0 ? -1 : 1))
  {
    int dirs = i == 0 ? dof%1000 : dof/1000;
    if (this->allDofs(dirs))   // all DOFs at this node are fixed,
      this->fix(1+iSnod,dirs); // local axes not needed
    else
    {
      // Create an extra node for the local DOFs. The new node, for which
      // the Dirichlet boundary conditions will be defined, then inherits
      // the global node number of the original node. The original node, which
      // do not enter the equation system, receives a new global node number.
      std::map<int,int>::const_iterator xit = xNode.find(MLGN[iSnod]);
      if (xit != xNode.end()) // this node has been processed by another patch
        myMLGN.push_back(xit->second);
      else
      {
        // Store the original-to-extra node number mapping in ASMstruct::xNode
        int iMnod = myMLGN.size();
        myMLGN.push_back(++gNod);
        xNode[MLGN[iSnod]] = gNod;
        std::swap(myMLGN[iMnod],myMLGN[iSnod]);

        xnMap[1+iMnod] = 1+iSnod; // Store nodal connection needed by getCoord
        nxMap[1+iSnod] = 1+iMnod; // Store nodal connection needed by getNodeID

        // Add Dirichlet condition on the local DOF(s) of the added node
        this->fix(1+iMnod,dirs);

        // Establish constraint equations relating the global and local DOFs
        this->addLocal2GlobalCpl(iSnod,MLGN[iMnod],Tlg);
        ++added;
      }
    }
  }

  return added;
}


bool ASMs1DC1::addRigidCpl (int lindx, int, int,
                            int& gMaster, const Vec3& Xmaster, bool extraPt)
{
  if (extraPt) // The master point is not a patch node, create an extra node
    extraPt = this->createRgdMasterNode(gMaster);

  for (int i = 0; i < 2; i++)
  {
    int iSn = lindx == 1 ? i : this->getSize()-1-i;
    Vec3 dX = this->getCoord(iSn+1) - Xmaster;
    if (nsd < 3 && nf == 1)
    {
      // Scalar beam problem (1 DOF per node -> 3 DOFs in master point)
      MPC* cons = new MPC(MLGN[iSn],1);
      if (this->addMPC(cons) && cons)
      {
        cons->addMaster(gMaster,1,1.0);
        if (nsd == 2)
          cons->addMaster(gMaster,2,dX.y);
        cons->addMaster(gMaster,3,-dX.x);
#if SP_DEBUG > 1
        std::cout <<"Added constraint: "<< *cons;
#endif
      }
    }
    else if (nsd == 3 && nf == 3)
      // Spatial beam problem (3 DOF per node -> 6 DOFs in master point)
      this->addRigidMPC(MLGN[iSn],gMaster,dX);
    else
    {
      std::cerr <<"ASMs1DC1::addRigidCpl: Invalid nsd,nf combination: "
                << (int)nsd <<","<< (int)nf << std::endl;
      return false;
    }
  }

  return extraPt;
}
