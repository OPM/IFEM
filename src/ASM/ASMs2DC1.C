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

#include "GoTools/geometry/SplineSurface.h"

#include "ASMs2DC1.h"
#include "Lagrange.h"
#include "Vec3Oper.h"
#include "Vec3.h"
#include "MPC.h"

//! \brief Solves A1j*xi*eta + A2j*ci + A3j*eta = A4j, j=1,2 for xi,eta.
extern "C" void dslbln_(const int& ipsw, const int& iwr, const double& eps,
                        const double* A, const double* B,
                        int& nsol, double* xi, double* eta);

std::map<int,ASMs2DC1*> ASMs2DC1::neighbors;


bool ASMs2DC1::generateFEMTopology ()
{
  if (surf->order_u() > 2 && surf->order_v() > 2)
    return this->ASMs2D::generateFEMTopology();

  std::cerr <<" *** ASMs2DC1::generateFEMTopology: The polynomial order ("
            << surf->order_u() <<","<< surf->order_v() <<") is too low.\n    "
            <<" C1-continuity requires at least quadratic order."<< std::endl;
  return false;
}


bool ASMs2DC1::connectC1 (int edge, ASMs2DC1* neighbor, int nedge, bool revers)
{
  if (shareFE && neighbor->shareFE)
    return true;
  else if (shareFE || neighbor->shareFE)
  {
    std::cerr <<" *** ASMs2DC1::connectPatch: Logic error, cannot"
              <<" connect a sharedFE patch with an unshared one"<< std::endl;
    return false;
  }

  // Set up the slave node numbers for this surface patch

  int n1, n2;
  if (!this->getSize(n1,n2)) return false;
  int node = 1, i1 = 0, i2 = 0;

  switch (edge)
    {
    case 2: // Positive I-direction
      node += n1-1;
    case 1: // Negative I-direction
      i1 = n1;
      n1 = n2;
      i2 = 3-edge*2;
      break;

    case 4: // Positive J-direction
      node += n1*(n2-1);
    case 3: // Negative J-direction
      i1 = 1;
      i2 = n1*(7-edge*2);
      break;

    default:
      std::cerr <<" *** ASMs2DC1::connectPatch: Invalid slave edge "
                << edge << std::endl;
      return false;
    }

  int i;
  IntMat slaveNodes(2,IntVec(n1,0));
  for (i = 0; i < n1; i++, node += i1)
  {
    slaveNodes.front()[i] = node;
    slaveNodes.back()[i] = node+i2;
  }

  // Set up the master node numbers for the neighboring surface patch

  if (!neighbor->getSize(n1,n2)) return false;
  node = 1; i1 = i2 = 0;

  switch (nedge)
    {
    case 2: // Positive I-direction
      node += n1-1;
    case 1: // Negative I-direction
      i1 = n1;
      n1 = n2;
      i2 = 3-nedge*2;
      break;

    case 4: // Positive J-direction
      node += n1*(n2-1);
    case 3: // Negative J-direction
      i1 = 1;
      i2 = n1*(7-nedge*2);
      break;

    default:
      std::cerr <<" *** ASMs2DC1::connectPatch: Invalid master edge "
                << nedge << std::endl;
      return false;
    }

  if (n1 != (int)slaveNodes.front().size())
  {
    std::cerr <<" *** ASMs2DC1::connectPatch: Non-matching edges, sizes "
              << n1 <<" and "<< slaveNodes.front().size() << std::endl;
    return false;
  }

  for (i = 0; i < n1; i++, node += i1)
  {
    int slave   = slaveNodes.front()[revers ? n1-i-1 : i];
    int master1 = node + i2;
    int master2 = slaveNodes.back()[revers ? n1-i-1 : i];
    if (neighbor->myMLGN[node-1] == myMLGN[slave-1])
    {
      for (unsigned char d = 1; d <= nf; d++)
        this->add3PC(MLGN[slave-1],d,neighbor->MLGN[master1-1],MLGN[master2-1]);
      neighbors[neighbor->MLGN[master1-1]] = neighbor;
      neighbors[MLGN[master2-1]] = this;
    }
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
      {
        this->makePeriodic(master,master+n1-1);
        for (unsigned char d = 1; d <= nf && n1 > 3; d++)
          this->add3PC(MLGN[master-1],d,MLGN[master],MLGN[master+n1-3]);
      }
      break;

    case 2: // Edges are closed in J-direction
      for (int i1 = 1; i1 <= n1; i1++, master++)
      {
        this->makePeriodic(master,master+n1*(n2-1));
        for (unsigned char d = 1; d <= nf && n2 > 3; d++)
          this->add3PC(MLGN[master-1],d,MLGN[master+n1-1],
                       MLGN[master+n1*(n2-2)-1]);
      }
      break;
    }
}


void ASMs2DC1::constrainEdge (int dir, bool open, int dof, int code, char)
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
        if (open && (i2 == 1 || i2 == n2))
          continue; // Skip the end points
        else if (dof%100)
          this->prescribe(node,dof%100,code);
        if (dof >= 100)
        {
          if (dof%100 && code == 0)
            // The edge is clamped, fix the neighboring node line
            this->prescribe(node-dir,dof/100,0);
          else
            // The edge has a prescribed rotation, add an MPC for that
            this->add2PC(MLGN[node-dir-1],dof/100,MLGN[node-1],code);
        }
      }
      break;

    case  2: // Back edge (positive J-direction)
      node += n1*(n2-1);
    case -2: // Front edge (negative J-direction)
      for (int i1 = 1; i1 <= n1; i1++, node++)
      {
        if (open && (i1 == 1 || i1 == n1))
          continue; // Skip the end points
        else if (dof%100)
          this->prescribe(node,dof%100,code);
        if (dof >= 100)
        {
          if (dof%100 && code == 0)
            // Edge is clamped, fix the neighboring node line
            this->prescribe(node-n1*dir/2,dof/100,0);
          else
            // The edge has a prescribed rotation, add an MPC for that
            this->add2PC(MLGN[node-n1*dir/2-1],MLGN[node-1],dof/100,code);
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


void ASMs2DC1::renumberNodes (const std::map<int,int>& old2new)
{
  std::map<int,int>::const_iterator it;
  std::map<int,ASMs2DC1*>::iterator nit;
  for (it = old2new.begin(); it != old2new.end(); it++)
    if (it->first != it->second)
      if ((nit = neighbors.find(it->first)) != neighbors.end())
      {
        neighbors[it->second] = nit->second;
        neighbors.erase(nit);
      }
}


/*!
  \brief Computes coupling coefficients for a 2-master constraint in a C1-patch.
*/

static void initMPC2 (MPC* mpc, std::vector<Vec3>& X)
{
  double s1 = (X[1]-X[0]).length();
  double s2 = (X[2]-X[0]).length();
  mpc->updateMaster(0,s1/(s1+s2));
  mpc->updateMaster(1,s2/(s1+s2));
}


/*!
  \brief Computes coupling coefficients for a 3-master constraint in a C1-patch.
*/

static void initMPC3 (MPC* mpc, std::vector<Vec3>& X)
{
  Vec3 V10(X[0]-X[1]);
  Vec3 V12(X[2]-X[1]);
  Vec3 V13(X[3]-X[1]);

  // Normal vector
  Vec3 VN(V12,V13);
  double Area = VN.length();
  VN /= Area;

  // Compute area coordinates
  V10 -= (VN*V10)*VN;
  Vec3 Aux1(V12,V10);
  Vec3 Aux2(V10,V13);
  double A3 = copysign(Aux1.length(),Aux1*VN);
  double A2 = copysign(Aux2.length(),Aux2*VN);
  mpc->updateMaster(0,(Area-A2-A3)/Area);
  mpc->updateMaster(1,A2/Area);
  mpc->updateMaster(2,A3/Area);
}


/*!
  \brief Computes coupling coefficients for a 4-master constraint in a C1-patch.
*/

static bool initMPC4flat (MPC* mpc, std::vector<Vec3>& X)
{
  // Check that all points are in the XY-plane
  double Zmin = X.front().z;
  double Zmax = Zmin;

  int i;
  for (i = 1; i <= 4; i++)
    if (X[i].z > Zmax)
      Zmax = X[i].z;
    else if (X[i].z < Zmin)
      Zmin = X[i].z;

  const double epsFlat = 1.0e-6;
  if (Zmax-Zmin > epsFlat) return false;

  // Order the master points counter-clockwise around the slave
  std::map<double,int> points;
  for (i = 1; i <= 4; i++)
  {
    double theta = atan2(X[i].y-X.front().y, X[i].x-X.front().x);
    if (theta < 0.0) theta += 2.0*M_PI;
    points[theta] = i;
  }

  double x[4], y[4], A[2][4];
  std::map<double,int>::const_iterator it;
  for (i = 0, it = points.begin(); it != points.end(); it++, i++)
  {
    x[i] = X[it->second].x - X.front().x;
    y[i] = X[it->second].y - X.front().y;
  }

  A[0][0] = 0.25*( x[0]-x[1]+x[2]-x[3]);
  A[0][1] = 0.25*(-x[0]+x[1]+x[2]-x[3]);
  A[0][2] = 0.25*(-x[0]-x[1]+x[2]+x[3]);
  A[0][3] = 0.25*(-x[0]-x[1]-x[2]-x[3]);

  A[1][0] = 0.25*( y[0]-y[1]+y[2]-y[3]);
  A[1][1] = 0.25*(-y[0]+y[1]+y[2]-y[3]);
  A[1][2] = 0.25*(-y[0]-y[1]+y[2]+y[3]);
  A[1][3] = 0.25*(-y[0]-y[1]-y[2]-y[3]);

  // We have to solve the following set of equations (j=1,2):
  // A1j*XI*ETA + A2j*XI + A3j*ETA = A4j
  const int IPSW = 0;
  const int IWR  = 6;
  dslbln_(IPSW,IWR,epsFlat,A[0],A[1],i,x,y);
  if (i != 1)
  {
    std::cerr <<" *** initMPC4flat: Detected "<< i <<" solutions."<< std::endl;
    return false;
  }

  // Evaluate the bi-linear interpolation polynomial
  RealArray N;
  Lagrange::computeBasis(N,2,x[0],2,y[0]);
  std::swap(N[2],N[3]);

  for (i = 0, it = points.begin(); it != points.end(); it++, i++)
    mpc->updateMaster(it->second-1,N[i]);

  return true;
}


bool ASMs2DC1::initConstraints ()
{
  // Compute coupling coefficients for the constraint equations
  // enforcing C1-continuity in the solution
  std::map<int,ASMs2DC1*>::const_iterator npit;
  for (MPCIter sit = mpcs.begin(); sit != mpcs.end(); sit++)
    if (dCode.find(*sit) == dCode.end())
    {
      // Extract coordinates of the control points involved, first the slave.
      // Note that some of the points are in a neighboring patch.
      ASMs2DC1* mpch = this;
      size_t nMaster = (*sit)->getNoMaster();
      std::vector<Vec3> X(1+nMaster);
      size_t inod = this->getNodeIndex((*sit)->getSlave().node);
      for (size_t m = 0; m < nMaster && inod > 0; m++)
      {
        X[m] = mpch->getCoord(inod);
        inod = (*sit)->getMaster(m).node;
        npit = neighbors.find(inod);
        mpch = npit == neighbors.end() ? this : npit->second;
        inod = mpch->getNodeIndex(inod);
      }
      if (inod > 0)
        X.back() = mpch->getCoord(inod);
      else
      {
        std::cerr <<" *** ASMs2DC1::initConstraints: Failed to initialize "
                  << **sit <<"     MPC masters:";
        for (npit = neighbors.begin(); npit != neighbors.end(); npit++)
          std::cerr <<" "<< npit->first;
        std::cerr << std::endl;
        return false;
      }

      switch (nMaster)
        {
        case 1:
          break;
        case 2:
          initMPC2(*sit,X);
          break;
        case 3:
          initMPC3(*sit,X);
          break;
        case 4:
          if (initMPC4flat(*sit,X))
            break;
        default:
          std::cerr <<" *** ASMs2DC1::initConstraints: Corner point with "
                    << nMaster << " connections not supported."<< std::endl;
          return false;
        }
    }

  return true;
}


bool ASMs2DC1::updateDirichlet (const std::map<int,RealFunc*>& func,
                                const std::map<int,VecFunc*>& vfunc,
                                double time, const std::map<int,int>* g2l)
{
  neighbors.clear();

  // Update the constraint equations defining Dirichlet conditions
  // on 1st-derivatives (normal rotations). Note: We here assume that
  // all MPCs with one master DOF and non-zero slave coefficient denote
  // prescribed 1st derivatives.
  std::map<int,RealFunc*>::const_iterator fit;
  std::map<int,VecFunc*>::const_iterator vfit;
  for (MPCMap::iterator cit = dCode.begin(); cit != dCode.end(); cit++)
    if (cit->first->getNoMaster() == 1)
    {
      size_t inod = this->getNodeIndex(cit->first->getSlave().node);
      size_t jnod = this->getNodeIndex(cit->first->getMaster(0).node);
      if (inod < 1 || jnod < 1) return false;

      // Evaluate the prescribed rotation
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
      // Update the slave coefficient, s = |X_j-X_i|*tan(theta)
      cit->first->setSlaveCoeff((this->getCoord(inod) - X).length()*tan(theta));
    }

  return this->ASMs2D::updateDirichlet(func,vfunc,time,g2l);
}
