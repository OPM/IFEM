// $Id$
//==============================================================================
//!
//! \file ASMs3DmxLag.C
//!
//! \date Dec 28 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of structured 3D Lagrange mixed FE models.
//!
//==============================================================================

#include "ASMs3DmxLag.h"
#include "Lagrange.h"
#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include <numeric>


ASMs3DmxLag::ASMs3DmxLag (const CharVec& n_f)
  : ASMs3DLag(std::accumulate(n_f.begin(),n_f.end(),0)), ASMmxBase(n_f)
{
}


ASMs3DmxLag::ASMs3DmxLag (const ASMs3DmxLag& patch, const CharVec& n_f)
  : ASMs3DLag(patch), ASMmxBase(n_f),
    nxx(patch.nxx), nyx(patch.nyx), nzx(patch.nzx)
{
}


void ASMs3DmxLag::clear (bool retainGeometry)
{
  nxx.clear();
  nyx.clear();
  nzx.clear();
  this->ASMs3DLag::clear(retainGeometry);
}


size_t ASMs3DmxLag::getNoNodes (int basis) const
{
  if (basis > 0 && basis <= (int)nb.size())
    return nb[basis-1];
  else
    return MLGN.size();
}


unsigned char ASMs3DmxLag::getNoFields (int basis) const
{
  if (basis == 0)
    return std::accumulate(nfx.begin(), nfx.end(), 0);

  return nfx[basis-1];
}


unsigned char ASMs3DmxLag::getNodalDOFs (size_t inod) const
{
  if (this->isLMn(inod)) return nLag;
  size_t nbc=0;
  for (size_t i=0;i<nb.size();++i)
    if (inod <= (nbc+=nb[i]))
      return nfx[i];

  return nfx[0];
}


char ASMs3DmxLag::getNodeType (size_t inod) const
{
  if (this->isLMn(inod)) return 'L';
  size_t nbc=0;
  for (size_t i=0;i<nb.size();++i)
    if (inod <= (nbc+=nb[i]))
      return i == 0 ? 'D' : 'P';

  return 'X';
}


void ASMs3DmxLag::initMADOF (const int* sysMadof)
{
  this->initMx(MLGN,sysMadof);
}


void ASMs3DmxLag::extractNodeVec (const RealArray& globRes, RealArray& nodeVec,
				  unsigned char, int basis) const
{
  this->extractNodeVecMx(globRes,nodeVec,basis);
}


bool ASMs3DmxLag::getSolution (Matrix& sField, const Vector& locSol,
			       const IntVec& nodes) const
{
  return this->getSolutionMx(sField,locSol,nodes);
}


bool ASMs3DmxLag::generateFEMTopology ()
{
  // Generate/check FE data for the geometry/basis1
  bool haveFEdata = !MLGN.empty();
  bool basis1IsOK = this->ASMs3DLag::generateFEMTopology();
  elmBasis = 1;
  if ((haveFEdata && !shareFE) || !basis1IsOK) return basis1IsOK;

  // Order of 2nd basis in the three parametric directions (order = degree + 1)
  const size_t q1 = p1 - 1;
  const size_t q2 = p2 - 1;
  const size_t q3 = p3 - 1;
  if (q1 < 2 || q2 < 2 || q3 < 2)
  {
    std::cerr <<" *** ASMs3DmxLag::generateFEMTopology: Too low order "<< q1
	      <<","<< q2 <<","<< q3 <<" for the second basis."<< std::endl;
    return false;
  }

  // Evaluate the parametric values
  RealArray gpar1, gpar2, gpar3;
  if (!this->getGridParameters(gpar1,0,q1-1)) return false;
  if (!this->getGridParameters(gpar2,1,q2-1)) return false;
  if (!this->getGridParameters(gpar3,2,q3-1)) return false;

  // Number of nodes in each direction for each basis
  nxx.resize(2);
  nyx.resize(2);
  nzx.resize(2);
  nxx[0] = nx;
  nyx[0] = ny;
  nzx[0] = nz;
  nxx[1] = gpar1.size();
  nyx[1] = gpar2.size();
  nzx[1] = gpar3.size();

  // Number of nodes per element in each direction for each basis
  elem_sizes.resize(2);
  elem_sizes[0][0] = p1;
  elem_sizes[0][1] = p2;
  elem_sizes[0][2] = p3;
  elem_sizes[1][0] = q1;
  elem_sizes[1][1] = q2;
  elem_sizes[1][2] = q3;

  // Number of nodes per element for each basis
  elem_size.resize(2);
  elem_size[0] = p1*p2*p3;
  elem_size[1] = q1*q2*q3;

  // Total number of nodes for each basis
  nb.resize(2);
  nb[0] = MLGN.size();
  nb[1] = nxx[1]*nyx[1]*nzx[1];

  if (shareFE == 'F') return true;

  // Add nodes for second basis (coordinates are not needed)
  nnod += nb[1];
  myMLGN.reserve(nnod);
  for (size_t i3 = 0; i3 < nzx[1]; i3++)
    for (size_t i2 = 0; i2 < nyx[1]; i2++)
      for (size_t i1 = 0; i1 < nxx[1]; i1++)
	myMLGN.push_back(++gNod);

  // Number of elements in each direction
  const int nelx = (nxx[1]-1)/(q1-1);
  const int nely = (nyx[1]-1)/(q2-1);
  const int nelz = (nzx[1]-1)/(q3-1);

  // Add connectivity for second basis: local --> global node relation
  int i, j, k, iel;
  for (k = iel = 0; k < nelz; k++)
    for (j = 0; j < nely; j++)
      for (i = 0; i < nelx; i++, iel++)
      {
	size_t nen1 = myMNPC[iel].size();
	myMNPC[iel].resize(nen1+q1*q2*q3);

	// First node in current element
	int corner = nb[0] + (q3-1)*nxx[1]*nyx[1]*k + (q2-1)*nxx[1]*j + q1*i-i;

	for (size_t c = 0; c < q3; c++)
	{
	  int cornod = nen1 + q1*q2*c;
	  myMNPC[iel][cornod] = corner + nxx[1]*nyx[1]*c;
	  for (size_t b = 1; b < q2; b++)
	  {
	    int facenod = cornod + q1*b;
	    myMNPC[iel][facenod] = myMNPC[iel][cornod] + nxx[1]*b;
	    for (size_t a = 1; a < q1; a++)
	    {
	      myMNPC[iel][facenod+a] = myMNPC[iel][facenod] + a;
	      myMNPC[iel][cornod+a]  = myMNPC[iel][cornod] + a;
	    }
	  }
	}
      }

  return true;
}


bool ASMs3DmxLag::connectPatch (int face, ASM3D& neighbor, int nface,
                                int norient, int basis,
                                bool coordCheck, int thick)
{
  ASMs3DmxLag* neighMx = dynamic_cast<ASMs3DmxLag*>(&neighbor);
  if (!neighMx) return false;

  if (swapW && face > 4) // Account for swapped parameter direction
    face = 11-face;

  if (neighMx->swapW && face > 4) // Account for swapped parameter direction
    nface = 11-nface;

  size_t nb1 = 0, nb2 = 0;
  for (size_t i = 1; i <= nxx.size(); i++) {
    if (basis == 0 || i == (size_t)basis)
      if (!this->connectBasis(face,*neighMx,nface,norient,i,nb1,nb2,
                              coordCheck,thick))
        return false;

    nb1 += nb[i-1];
    nb2 += neighMx->nb[i-1];
  }

  this->addNeighbor(neighMx);
  return true;
}


void ASMs3DmxLag::closeBoundaries (int dir, int, int)
{
  size_t nbi = 1;
  for (size_t i = 0; i < nxx.size(); nbi += nb[i++])
    this->ASMs3D::closeBoundaries(dir,1+i,nbi);
}


bool ASMs3DmxLag::getSize (int& n1, int& n2, int& n3, int basis) const
{
  if (basis <= 1)
    return this->ASMs3DLag::getSize(n1,n2,n3,1);

  n1 = nxx[1];
  n2 = nyx[1];
  n3 = nzx[1];

  return true;
}


bool ASMs3DmxLag::integrate (Integrand& integrand,
			     GlobalIntegral& glInt,
			     const TimeDomain& time)
{
  if (this->empty()) return true; // silently ignore empty patches

  if (myCache.empty()) {
    myCache.emplace_back(std::make_unique<BasisFunctionCache>(*this, cachePolicy, 1));
    const BasisFunctionCache& front = static_cast<const BasisFunctionCache&>(*myCache.front());
    for (size_t b = 2; b <= this->getNoBasis(); ++b)
      myCache.emplace_back(std::make_unique<BasisFunctionCache>(front, b));
  }

  for (std::unique_ptr<ASMs3D::BasisFunctionCache>& cache : myCache) {
    cache->setIntegrand(&integrand);
    cache->init(1);
  }

  BasisFunctionCache& cache = static_cast<BasisFunctionCache&>(*myCache.front());

  // Get Gaussian quadrature points and weights
  const std::array<int,3>& ng = cache.nGauss();
  const std::array<const double*,3>& xg = cache.coord();
  const std::array<const double*,3>& wg = cache.weight();

  // Number of elements in each direction
  const int nelx = cache.noElms()[0];
  const int nely = cache.noElms()[1];


  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t g = 0; g < threadGroupsVol.size() && ok; g++)
  {
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < threadGroupsVol[g].size(); t++)
    {
      MxFiniteElement fe(elem_size);
      Matrices dNxdu;
      Matrix Xnod, Jac;
      Vec4   X(nullptr,time.t);
      for (size_t l = 0; l < threadGroupsVol[g][t].size() && ok; l++)
      {
        int iel = threadGroupsVol[g][t][l];
        int i1  =  iel % nelx;
        int i2  = (iel / nelx) % nely;
        int i3  =  iel / (nelx*nely);

        // Set up control point coordinates for current element
        if (!this->getElementCoordinates(Xnod,++iel))
	{
          ok = false;
          break;
        }

        // Initialize element quantities
        fe.iel = MLGE[iel-1];
        LocalIntegral* A = integrand.getLocalIntegral(elem_size,fe.iel,false);
        if (!integrand.initElement(MNPC[iel-1], elem_size, nb, *A))
        {
          A->destruct();
          ok = false;
          break;
        }


        // --- Integration loop over all Gauss points in each direction --------

        size_t ip = 0;
        int jp = ((i3*nely + i2)*nelx + i1)*ng[0]*ng[1]*ng[2];
        fe.iGP = firstIp + jp; // Global integration point counter

        for (int k = 0; k < ng[2]; k++)
          for (int j = 0; j < ng[1]; j++)
            for (int i = 0; i < ng[0]; i++, fe.iGP++)
            {
              // Parameter value of current integration point
              fe.u = cache.getParam(0,i1,i);
              fe.v = cache.getParam(1,i2,j);
              fe.w = cache.getParam(2,i3,k);

              // Local coordinate of current integration point
              fe.xi   = xg[0][i];
              fe.eta  = xg[1][j];
              fe.zeta = xg[2][k];

              // Compute basis function derivatives at current integration point
              std::vector<const BasisFunctionVals*> bfs(this->getNoBasis());
              for (size_t b = 0; b < this->getNoBasis(); ++b) {
                bfs[b] = &myCache[b]->getVals(iel-1,ip);
                fe.basis(b+1) = bfs[b]->N;
              }

              // Compute Jacobian inverse of coordinate mapping and derivatives
              if (!fe.Jacobian(Jac,Xnod,elmBasis,&bfs))
                continue; // skip singular points

              // Cartesian coordinates of current integration point
              X.assign(Xnod * fe.basis(elmBasis));

              // Evaluate the integrand and accumulate element contributions
              fe.detJxW *= wg[0][i]*wg[1][j]*wg[2][k];
              if (!integrand.evalIntMx(*A,fe,time,X))
                ok = false;
            }

        // Finalize the element quantities
        if (ok && !integrand.finalizeElement(*A,fe,time,firstIp+jp))
          ok = false;

        // Assembly of global system integral
        if (ok && !glInt.assemble(A->ref(),fe.iel))
          ok = false;

        A->destruct();
      }
    }
  }

  return ok;
}


bool ASMs3DmxLag::integrate (Integrand& integrand, int lIndex,
			     GlobalIntegral& glInt,
			     const TimeDomain& time)
{
  if (this->empty()) return true; // silently ignore empty patches

  std::map<char,ThreadGroups>::const_iterator tit;
  if ((tit = threadGroupsFace.find(lIndex%10)) == threadGroupsFace.end())
  {
    std::cerr <<" *** ASMs3DmxLag::integrate: No thread groups for face "
              << lIndex%10 << std::endl;
    return false;
  }
  const ThreadGroups& threadGrp = tit->second;

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Find the parametric direction of the face normal {-3,-2,-1, 1, 2, 3}
  const int faceDir = (lIndex%10+1)/((lIndex%2) ? -2 : 2);

  const int t0 = abs(faceDir); // unsigned normal direction of the face
  const int t1 = 1 + t0%3; // first tangent direction of the face
  const int t2 = 1 + t1%3; // second tangent direction of the face

  // Extract the Neumann order flag (1 or higher) for the integrand
  integrand.setNeumannOrder(1 + lIndex/10);

  // Number of elements in each direction
  const int nel1 = (nxx[elmBasis-1]-1)/(elem_sizes[elmBasis-1][0]-1);
  const int nel2 = (nyx[elmBasis-1]-1)/(elem_sizes[elmBasis-1][1]-1);

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex%10);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;

  // === Assembly loop over all elements on the patch face =====================

  bool ok = true;
  for (size_t g = 0; g < threadGrp.size() && ok; g++)
  {
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < threadGrp[g].size(); t++)
    {
      MxFiniteElement fe(elem_size);
      Matrices dNxdu;
      Matrix Xnod, Jac;
      Vec4   X(nullptr,time.t);
      Vec3   normal;
      double xi[3];

      for (size_t l = 0; l < threadGrp[g][t].size() && ok; l++)
      {
        int iel = threadGrp[g][t][l];
        int i1  =  iel % nel1;
        int i2  = (iel / nel1) % nel2;
        int i3  =  iel / (nel1*nel2);

	// Set up control point coordinates for current element
	if (!this->getElementCoordinates(Xnod,++iel))
        {
          ok = false;
          break;
        }

	// Initialize element quantities
	fe.iel = MLGE[iel-1];
        LocalIntegral* A = integrand.getLocalIntegral(elem_size,fe.iel,true);
	if (!integrand.initElementBou(MNPC[iel-1],elem_size,nb,*A))
        {
          A->destruct();
          ok = false;
          break;
        }

        // Define some loop control variables depending on which face we are on
        int nf1, j1, j2;
        switch (abs(faceDir))
        {
          case 1: nf1 = nel2; j2 = i3; j1 = i2; break;
          case 2: nf1 = nel1; j2 = i3; j1 = i1; break;
          case 3: nf1 = nel1; j2 = i2; j1 = i1; break;
          default: nf1 = j1 = j2 = 0;
        }


	// --- Integration loop over all Gauss points in each direction --------

        int jp = (j2*nf1 + j1)*nGauss*nGauss;
        fe.iGP = firstp + jp; // Global integration point counter

	for (int j = 0; j < nGauss; j++)
	  for (int i = 0; i < nGauss; i++, fe.iGP++)
	  {
	    // Gauss point coordinates on the face
	    xi[t0-1] = faceDir < 0 ? -1.0 : 1.0;
	    xi[t1-1] = xg[i];
	    xi[t2-1] = xg[j];

	    // Compute the basis functions and their derivatives, using
	    // tensor product of one-dimensional Lagrange polynomials
            for (size_t b = 0; b < nxx.size(); ++b)
              if (!Lagrange::computeBasis(fe.basis(b+1),dNxdu[b],
                                          elem_sizes[b][0],xi[0],
                                          elem_sizes[b][1],xi[1],
                                          elem_sizes[b][2],xi[2]))
                ok = false;

	    // Compute basis function derivatives and the edge normal
        fe.detJxW = utl::Jacobian(Jac,normal,fe.grad(elmBasis),Xnod,
                                      dNxdu[elmBasis-1],t1,t2);
	    if (fe.detJxW == 0.0) continue; // skip singular points

            for (size_t b = 0; b < nxx.size(); ++b)
              if (b != (size_t)elmBasis-1)
                fe.grad(b+1).multiply(dNxdu[b],Jac);

	    if (faceDir < 0) normal *= -1.0;

	    // Cartesian coordinates of current integration point
        X.assign(Xnod * fe.basis(elmBasis));

	    // Evaluate the integrand and accumulate element contributions
	    fe.detJxW *= wg[i]*wg[j];
	    if (!integrand.evalBouMx(*A,fe,time,X,normal))
              ok = false;
	  }

        // Finalize the element quantities
        if (ok && !integrand.finalizeElementBou(*A,fe,time))
          ok = false;

	// Assembly of global system integral
	if (ok && !glInt.assemble(A->ref(),fe.iel))
          ok = false;

        A->destruct();
      }
    }
  }

  return ok;
}


bool ASMs3DmxLag::evalSolution (Matrix& sField, const Vector& locSol,
                                const RealArray*, bool, int, int nf) const
{
  size_t nc1 = nf ? nf : nfx[0];
  size_t nc2 = 0;
  if (nc1*nb[0] < locSol.size())
    nc2 = (locSol.size() - nc1*nb[0])/nb[1];
  else
    nc1 = locSol.size()/nb[0];

  if (nc1*nb[0] + nc2*nb[1] != locSol.size())
    return false;

  // TODO: Add evaluation of the second field at the nodes of the first field
  size_t nPoints = nb[0];
  size_t nComp = nc1;
  size_t i, n, ip = 0;
  sField.resize(nComp,nPoints);
  for (n = 1; n <= nPoints; n++)
    for (i = 1; i <= nComp; i++)
      sField(i,n) = locSol(++ip);

  return true;
}


bool ASMs3DmxLag::evalSolution (Matrix& sField, const IntegrandBase& integrand,
				const RealArray*, bool) const
{
  sField.resize(0,0);

  double incx = 2.0/double(p1-1);
  double incy = 2.0/double(p2-1);
  double incz = 2.0/double(p3-1);

  size_t nPoints = nb[0];
  IntVec check(nPoints,0);

  MxFiniteElement fe(elem_size);
  Vector          solPt;
  Vectors         globSolPt(nPoints);
  Matrices        dNxdu(nxx.size());
  Matrix          Xnod, Jac;

  // Evaluate the secondary solution field at each point
  for (size_t iel = 1; iel <= nel; iel++)
  {
    IntVec::const_iterator f2start = elmBasis == 1? MNPC[iel-1].begin() :
                                     MNPC[iel-1].begin() +
                                     std::accumulate(elem_size.begin()+elmBasis-2,
                                                     elem_size.begin()+elmBasis-1, 0);
    IntVec::const_iterator f2end = f2start + elem_size[elmBasis-1];
    IntVec mnpc1(f2start,f2end);

    this->getElementCoordinates(Xnod,iel);

    int i, j, k, loc = 0;
    for (k = 0; k < p3; k++)
      for (j = 0; j < p2; j++)
	for (i = 0; i < p1; i++, loc++)
	{
	  fe.xi   = -1.0 + i*incx;
	  fe.eta  = -1.0 + j*incy;
	  fe.zeta = -1.0 + k*incz;
          for (size_t b = 0; b < nxx.size(); ++b)
            if (!Lagrange::computeBasis(fe.basis(b+1),dNxdu[b],
                                        elem_sizes[b][0],fe.xi,
                                        elem_sizes[b][1],fe.eta,
                                        elem_sizes[b][2],fe.zeta))
              return false;

          // Compute Jacobian inverse of the coordinate mapping and
          // basis function derivatives w.r.t. Cartesian coordinates
          if (!fe.Jacobian(Jac,Xnod,elmBasis,nullptr,&dNxdu))
            continue; // skip singular points

	  // Now evaluate the solution field
      if (!integrand.evalSol(solPt,fe,Xnod*fe.basis(elmBasis),
                                 MNPC[iel-1],elem_size,nb))
	    return false;
	  else if (sField.empty())
	    sField.resize(solPt.size(),nPoints,true);

	  if (++check[mnpc1[loc]] == 1)
	    globSolPt[mnpc1[loc]] = solPt;
	  else
	    globSolPt[mnpc1[loc]] += solPt;
	}
  }

  for (size_t i = 0; i < nPoints; i++)
    sField.fillColumn(1+i,globSolPt[i] /= check[i]);

  return true;
}


ASMs3DmxLag::BasisFunctionCache::BasisFunctionCache (const ASMs3DLag& pch,
                                                     ASM::CachePolicy plcy,
                                                     int b) :
  ASMs3DLag::BasisFunctionCache(pch,plcy,b)
{
}


ASMs3DmxLag::BasisFunctionCache::BasisFunctionCache (const BasisFunctionCache& cache,
                                                     int b) :
  ASMs3DLag::BasisFunctionCache(cache,b)
{
}


BasisFunctionVals ASMs3DmxLag::BasisFunctionCache::calculatePt (size_t el,
                                                                size_t gp,
                                                                bool reduced) const
{
  std::array<size_t,3> gpIdx = this->gpIndex(gp,reduced);
  const Quadrature& q = reduced ? *reducedQ : *mainQ;

  const ASMs3DmxLag& pch = static_cast<const ASMs3DmxLag&>(patch);

  BasisFunctionVals result;
  if (nderiv == 1)
    Lagrange::computeBasis(result.N,result.dNdu,
                           pch.elem_sizes[basis-1][0],q.xg[0][gpIdx[0]],
                           pch.elem_sizes[basis-1][1],q.xg[1][gpIdx[1]],
                           pch.elem_sizes[basis-1][2],q.xg[2][gpIdx[2]]);

  return result;
}
