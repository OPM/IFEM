// $Id$
//==============================================================================
//!
//! \file ASMs2DC1.C
//!
//! \date Oct 25 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of C1-continuous structured 2D spline FE models.
//!
//==============================================================================

#include "ASMs2DC1.h"
#include "Vec3Oper.h"
#include "Vec3.h"
#include "MPC.h"


bool ASMs2DC1::connectPatch (int edge, ASMs2D& neighbor, int nedge, bool revers)
{
  ASMs2DC1* neigh = dynamic_cast<ASMs2DC1*>(&neighbor);
  if (!neigh) return false;

  if (shareFE && neigh->shareFE)
    return true;
  else if (shareFE || neigh->shareFE)
  {
    std::cerr <<" *** ASMs2DC1::connectPatch: Logic error, cannot"
	      <<" connect a sharedFE patch with an unshared one"<< std::endl;
    return false;
  }

  // Set up the slave node numbers for this surface patch

  int n1, n2;
  if (!this->getSize(n1,n2)) return false;
  int node = 1, i1 = 0;

  switch (edge)
    {
    case 2: // Positive I-direction
      node += n1-1;
    case 1: // Negative I-direction
      i1 = n1;
      n1 = n2;
      break;

    case 4: // Positive J-direction
      node += n1*(n2-1);
    case 3: // Negative J-direction
      i1 = 1;
      break;

    default:
      std::cerr <<" *** ASMs2DC1::connectPatch: Invalid slave edge "
		<< edge << std::endl;
      return false;
    }

  int i;
  IntVec slaveNodes(n1,0);
  for (i = 0; i < n1; i++, node += i1)
    slaveNodes[i] = node;

  // Set up the master node numbers for the neighboring surface patch

  if (!neigh->getSize(n1,n2)) return false;
  node = 1; i1 = 0;

  switch (nedge)
    {
    case 2: // Positive I-direction
      node += n1-1;
    case 1: // Negative I-direction
      i1 = n1;
      n1 = n2;
      break;

    case 4: // Positive J-direction
      node += n1*(n2-1);
    case 3: // Negative J-direction
      i1 = 1;
      break;

    default:
      std::cerr <<" *** ASMs2DC1::connectPatch: Invalid master edge "
		<< nedge << std::endl;
      return false;
    }

  if (n1 != (int)slaveNodes.size())
  {
    std::cerr <<" *** ASMs2DC1::connectPatch: Non-matching edges, sizes "
	      << n1 <<" and "<< slaveNodes.size() << std::endl;
    return false;
  }

  const double xtol = 1.0e-4;
  for (i = 0; i < n1; i++, node += i1)
  {
    int k = revers ? n1-i-1 : i;
    if (!neigh->getCoord(node).equal(this->getCoord(slaveNodes[k]),xtol))
    {
      std::cerr <<" *** ASMs2DC1::connectPatch: Non-matching nodes "
		<< node <<": "<< neigh->getCoord(node)
		<<"\n                                          and "
		<< slaveNodes[k] <<": "<< this->getCoord(slaveNodes[k])
		<< std::endl;
      return false;
    }
    else
      ASMbase::collapseNodes(neigh->myMLGN[node-1],myMLGN[slaveNodes[k]-1]);
  }

  return true;
}


void ASMs2DC1::closeEdges (int dir, int, int)
{
  int n1, n2;
  if (!this->getSize(n1,n2)) return;

  int master = 1;
  switch (dir)
    {
    case 1: // Edges are closed in I-direction
      for (int i2 = 1; i2 <= n2; i2++, master += n1)
	this->makePeriodic(master,master+n1-1);
      break;

    case 2: // Edges are closed in J-direction
      for (int i1 = 1; i1 <= n1; i1++, master++)
	this->makePeriodic(master,master+n1*(n2-1));
      break;
    }
}


void ASMs2DC1::constrainEdge (int dir, int dof, int code)
{
  int n1, n2, node = 1;
  if (!this->getSize(n1,n2)) return;

  switch (dir)
    {
    case  1: // Right edge (positive I-direction)
      node += n1-1;
    case -1: // Left edge (negative I-direction)
      for (int i2 = 1; i2 <= n2; i2++, node += n1)
      {
	if (dof%100)
	  this->prescribe(node,dof%100,code);
	if (dof >= 100)
	{
	  if (dof%100 && code == 0)
	    // The edge is clamped, fix the neighboring node line
	    this->prescribe(node-dir,dof/100,0);
	  else
	    // The edge has a prescribed rotation, add an MPC for that
	    this->addMPC(MLGN[node-dir-1],dof/100,MLGN[node-1],code);
	}
      }
      break;

    case  2: // Back edge (positive J-direction)
      node += n1*(n2-1);
    case -2: // Front edge (negative J-direction)
      for (int i1 = 1; i1 <= n1; i1++, node++)
      {
	if (dof%100)
	  this->prescribe(node,dof%100,code);
	if (dof >= 100)
	{
	  if (dof%100 && code == 0)
	    // Edge is clamped, fix the neighboring node line
	    this->prescribe(node-n1*dir/2,dof/100,0);
	  else
	    // The edge has a prescribed rotation, add an MPC for that
	    this->addMPC(MLGN[node-n1*dir/2-1],MLGN[node-1],dof/100,code);
	}
      }
      break;
    }
}


void ASMs2DC1::constrainCorner (int I, int J, int dof, int code)
{
  int n1, n2;
  if (!this->getSize(n1,n2)) return;

  int node = 1;
  if (I > 0) node += n1-1;
  if (J > 0) node += n1*(n2-1);

  this->prescribe(node,dof%100,code);
  if (dof < 100 || code > 0) return;

  // Also fix the two neighboring nodes
  if (I > 0)
    this->prescribe(node-1,dof%100,0);
  else
    this->prescribe(node+1,dof%100,0);

  if (J > 0)
    this->prescribe(node-n1,dof%100,0);
  else
    this->prescribe(node+n1,dof%100,0);
}


void ASMs2DC1::constrainNode (double xi, double eta, int dof, int code)
{
  if (xi  < 0.0 || xi  > 1.0) return;
  if (eta < 0.0 || eta > 1.0) return;

  int n1, n2;
  if (!this->getSize(n1,n2)) return;

  int I = int(0.5+(n1-1)*xi);
  int J = int(0.5+(n2-1)*eta);
  this->prescribe(n1*J+I+1,dof%100,code);
  if (dof < 100 || code > 0) return;

  // Also fix the (up to four) neighboring nodes
  if (I > 0)    this->prescribe(n1*J+I      ,dof%100,0);
  if (I < n1-1) this->prescribe(n1*J+I+2    ,dof%100,0);
  if (J > 0)    this->prescribe(n1*(J-1)+I+1,dof%100,0);
  if (J < n2-1) this->prescribe(n1*(J+1)+I+1,dof%100,0);
}


bool ASMs2DC1::updateDirichlet (const std::map<int,RealFunc*>& func,
				const std::map<int,VecFunc*>& vfunc,
				double time)
{
  std::map<int,RealFunc*>::const_iterator fit;
  std::map<int,VecFunc*>::const_iterator vfit;
  for (MPCMap::iterator cit = dCode.begin(); cit != dCode.end(); cit++)
    if (cit->first->getNoMaster() == 1)
    {
      size_t inod = this->getNodeIndex(cit->first->getSlave().node);
      size_t jnod = this->getNodeIndex(cit->first->getMaster(0).node);
      if (inod < 1 || jnod < 1) return false;

      double theta = 0.0;
      Vec4 X(this->getCoord(jnod),time);
      if ((fit = func.find(cit->second)) != func.end())
	theta = (*fit->second)(X);
      else if ((vfit = vfunc.find(cit->second)) != vfunc.end())
	theta = (*vfit->second)(X)[cit->first->getSlave().dof-1];
      else
      {
	std::cerr <<" *** ASMs2DC1::updateDirichlet: Code "<< cit->second
		  <<" is not associated with any function."<< std::endl;
	return false;
      }
      cit->first->setSlaveCoeff((this->getCoord(inod) - X).length()*tan(theta));
    }

  return this->ASMs2D::updateDirichlet(func,vfunc,time);
}
