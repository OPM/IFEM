// $Id$
//==============================================================================
//!
//! \file ASMs2DmxLag.C
//!
//! \date Dec 28 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of structured 2D Lagrange mixed FE models.
//!
//==============================================================================

#include "ASMs2DmxLag.h"
#include "Lagrange.h"
#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include <numeric>


ASMs2DmxLag::ASMs2DmxLag (unsigned char n_s, const CharVec& n_f)
  : ASMs2DLag(n_s), ASMmxBase(n_f)
{
}


ASMs2DmxLag::ASMs2DmxLag (const ASMs2DmxLag& patch, const CharVec& n_f)
  : ASMs2DLag(patch), ASMmxBase(n_f),
    nxx(patch.nxx), nyx(patch.nyx)
{
}


void ASMs2DmxLag::clear (bool retainGeometry)
{
  nxx.clear();
  nyx.clear();
  this->ASMs2DLag::clear(retainGeometry);
}


size_t ASMs2DmxLag::getNoNodes (int basis) const
{
  if (basis > 0 && basis <= (int)nb.size())
    return nb[basis-1];
  else
    return MLGN.size();
}


unsigned char ASMs2DmxLag::getNoFields (int basis) const
{
  if (basis == 0)
    return std::accumulate(nfx.begin(), nfx.end(), 0);

  return nfx[basis-1];
}


unsigned char ASMs2DmxLag::getNodalDOFs (size_t inod) const
{
  if (this->isLMn(inod)) return nLag;
  size_t nbc=0;
  for (size_t i=0;i<nb.size();++i)
    if (inod <= (nbc+=nb[i]))
      return nfx[i];

  return nfx[0];
}


char ASMs2DmxLag::getNodeType (size_t inod) const
{
  if (this->isLMn(inod)) return 'L';
  size_t nbc=0;
  for (size_t i=0;i<nb.size();++i)
    if (inod <= (nbc+=nb[i]))
      return i == 0 ? 'D' : 'P';

  return 'X';
}


void ASMs2DmxLag::initMADOF (const int* sysMadof)
{
  this->initMx(MLGN,sysMadof);
}


void ASMs2DmxLag::extractNodeVec (const RealArray& globRes, RealArray& nodeVec,
				  unsigned char, int basis) const
{
  this->extractNodeVecMx(globRes,nodeVec,basis);
}


bool ASMs2DmxLag::getSolution (Matrix& sField, const Vector& locSol,
			       const IntVec& nodes) const
{
  return this->getSolutionMx(sField,locSol,nodes);
}


bool ASMs2DmxLag::generateFEMTopology ()
{
  // Generate/check FE data for the geometry/basis1
  bool haveFEdata = !MLGN.empty();
  bool basis1IsOK = this->ASMs2DLag::generateFEMTopology();
  itgBasis = 1;
  if ((haveFEdata && !shareFE) || !basis1IsOK) return basis1IsOK;

  // Order of 2nd basis in the two parametric directions (order = degree + 1)
  const size_t q1 = p1 - 1;
  const size_t q2 = p2 - 1;
  if (q1 < 2 || q2 < 2)
  {
    std::cerr <<" *** ASMs2DmxLag::generateFEMTopology: Too low order "
              << q1 <<","<< q2 <<" for the second basis."<< std::endl;
    return false;
  }

  // Evaluate the parametric values
  RealArray gpar1, gpar2;
  if (!this->getGridParameters(gpar1,0,q1-1)) return false;
  if (!this->getGridParameters(gpar2,1,q2-1)) return false;

  // Number of nodes in each direction for each basis
  nxx.resize(2);
  nyx.resize(2);
  nxx[0] = nx;
  nyx[0] = ny;
  nxx[1] = gpar1.size();
  nyx[1] = gpar2.size();

  // Number of nodes per element in each direction for each basis
  elem_sizes.resize(2);
  elem_sizes[0][0] = p1;
  elem_sizes[0][1] = p2;
  elem_sizes[1][0] = q1;
  elem_sizes[1][1] = q2;

  // Number of nodes per element for each basis
  elem_size.resize(2);
  elem_size[0] = p1*p2;
  elem_size[1] = q1*q2;

  // Total number of nodes for each basis
  nb.resize(2);
  nb[0] = MLGN.size();
  nb[1] = nxx[1]*nyx[1];

  if (shareFE == 'F') return true;

  // Add nodes for second basis (coordinates are not needed)
  nnod += nb[1];
  myMLGN.reserve(nnod);
  for (size_t i2 = 0; i2 < nyx[1]; i2++)
    for (size_t i1 = 0; i1 < nxx[1]; i1++)
      myMLGN.push_back(++gNod);

  // Number of elements in each direction
  const int nelx = (nxx[1]-1)/(q1-1);
  const int nely = (nyx[1]-1)/(q2-1);

  // Add connectivity for second basis: local --> global node relation
  int i, j, iel;
  for (j = iel = 0; j < nely; j++)
    for (i = 0; i < nelx; i++, iel++)
    {
      size_t nen1 = myMNPC[iel].size();
      myMNPC[iel].resize(nen1+q1*q2);

      // First node in current element
      int corner = nb[0] + (q2-1)*nxx[1]*j + (q1-1)*i;

      for (size_t b = 0; b < q2; b++)
      {
	int facenod = nen1 + b*q1;
	myMNPC[iel][facenod] = corner + b*nxx[1];
	for (size_t a = 1; a < q1; a++)
	  myMNPC[iel][facenod+a] = myMNPC[iel][facenod] + a;
      }
    }

  return true;
}


bool ASMs2DmxLag::connectPatch (int edge, ASM2D& neighbor, int nedge, bool revers,
                                int basis, bool coordCheck, int thick)
{
  ASMs2DmxLag* neighMx = dynamic_cast<ASMs2DmxLag*>(&neighbor);
  if (!neighMx) return false;

  size_t nb1 = 0, nb2 = 0;
  for (size_t i = 1; i <= nxx.size(); i++) {
    if (basis == 0 || i == (size_t)basis)
      if (!this->connectBasis(edge,*neighMx,nedge,revers,i,nb1,nb2,
                              coordCheck,thick))
        return false;
    nb1 += nb[i-1];
    nb2 += neighMx->nb[i-1];
  }

  this->addNeighbor(neighMx);
  return true;
}


void ASMs2DmxLag::closeBoundaries (int dir, int, int)
{
  size_t nbi = 1;
  for (size_t i = 0; i < nxx.size(); nbi += nb[i++])
    this->ASMs2D::closeBoundaries(dir,1+i,nbi);
}


bool ASMs2DmxLag::getSize (int& n1, int& n2, int basis) const
{
  if (basis <= 1)
    return this->ASMs2DLag::getSize(n1,n2,1);

  n1 = nxx[basis-2];
  n2 = nyx[basis-2];

  return true;
}


bool ASMs2DmxLag::integrate (Integrand& integrand,
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

  for (std::unique_ptr<ASMs2D::BasisFunctionCache>& cache : myCache) {
    cache->setIntegrand(&integrand);
    cache->init(1);
  }

  BasisFunctionCache& cache = static_cast<BasisFunctionCache&>(*myCache.front());

  // Get Gaussian quadrature points and weights
  const std::array<int,2>& ng = cache.nGauss();
  const std::array<const double*,2>& xg = cache.coord();
  const std::array<const double*,2>& wg = cache.weight();

  const int nelx = cache.noElms()[0];


  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t g = 0; g < threadGroups.size() && ok; g++)
  {
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < threadGroups[g].size(); t++)
    {
      MxFiniteElement fe(elem_size);
      Matrix Xnod, Jac;
      Vec4   X(nullptr,time.t);
      for (size_t e = 0; e < threadGroups[g][t].size() && ok; ++e)
      {
        int iel = threadGroups[g][t][e];
        int i1  = iel % nelx;
        int i2  = iel / nelx;

        // Set up control point coordinates for current element
        if (!this->getElementCoordinates(Xnod,++iel))
        {
          ok = false;
          break;
        }

        // Initialize element quantities
        fe.iel = MLGE[iel-1];
        LocalIntegral* A = integrand.getLocalIntegral(elem_size,fe.iel,false);
        if (!integrand.initElement(MNPC[iel-1],elem_size,nb,*A))
        {
          A->destruct();
          ok = false;
          break;
        }

        // --- Integration loop over all Gauss points in each direction --------

        size_t ip = 0;
        int jp = (i2*nelx + i1)*ng[0]*ng[1];
        fe.iGP = firstIp + jp; // Global integration point counter

        for (int j = 0; j < ng[1]; j++)
          for (int i = 0; i < ng[0]; i++, fe.iGP++, ++ip)
          {
            // Parameter value of current integration point
            fe.u = cache.getParam(0,i1,i);
            fe.v = cache.getParam(1,i2,j);

            // Local coordinates of current integration point
            fe.xi  = xg[0][i];
            fe.eta = xg[1][j];

            // Fetch basis function derivatives at current integration point
            std::vector<const BasisFunctionVals*> bfs(this->getNoBasis());
            for (size_t b = 0; b < this->getNoBasis(); ++b) {
              bfs[b] = &myCache[b]->getVals(iel-1,ip);
              fe.basis(b+1) = bfs[b]->N;
            }

            // Compute Jacobian inverse of coordinate mapping and derivatives
            if (!fe.Jacobian(Jac,Xnod,itgBasis,&bfs))
              continue; // skip singular points

            // Cartesian coordinates of current integration point
            X.assign(Xnod * fe.basis(itgBasis));

            // Evaluate the integrand and accumulate element contributions
            fe.detJxW *= wg[0][i]*wg[1][j];
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


bool ASMs2DmxLag::integrate (Integrand& integrand, int lIndex,
			     GlobalIntegral& glInt,
			     const TimeDomain& time)
{
  if (this->empty()) return true; // silently ignore empty patches

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Find the parametric direction of the edge normal {-2,-1, 1, 2}
  const int edgeDir = (lIndex%10+1)/((lIndex%2) ? -2 : 2);

  const int t1 = abs(edgeDir); // tangent direction normal to the patch edge
  const int t2 = 3-t1;         // tangent direction along the patch edge

  // Extract the Neumann order flag (1 or higher) for the integrand
  integrand.setNeumannOrder(1 + lIndex/10);


  // Number of elements in each direction
  const int nelx = (nxx[itgBasis-1]-1)/(elem_sizes[itgBasis-1][0]-1);
  const int nely = (nyx[itgBasis-1]-1)/(elem_sizes[itgBasis-1][1]-1);

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex%10);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;

  MxFiniteElement fe(elem_size);
  Matrices dNxdu(nxx.size());
  Matrix Xnod, Jac;
  Vec4   X(nullptr,time.t);
  Vec3   normal;
  double xi[2];

  // === Assembly loop over all elements on the patch edge =====================

  int iel = 1;
  for (int i2 = 0; i2 < nely; i2++)
    for (int i1 = 0; i1 < nelx; i1++, iel++)
    {
      // Skip elements that are not on current boundary edge
      bool skipMe = false;
      switch (edgeDir)
	{
	case -1: if (i1 > 0)      skipMe = true; break;
	case  1: if (i1 < nelx-1) skipMe = true; break;
	case -2: if (i2 > 0)      skipMe = true; break;
	case  2: if (i2 < nely-1) skipMe = true; break;
	}
      if (skipMe) continue;

      // Set up control point coordinates for current element
      if (!this->getElementCoordinates(Xnod,iel)) return false;

      // Initialize element quantities
      fe.iel = MLGE[iel-1];
      LocalIntegral* A = integrand.getLocalIntegral(elem_size,fe.iel,true);
      bool ok = integrand.initElementBou(MNPC[iel-1],elem_size,nb,*A);

      // --- Integration loop over all Gauss points along the edge -------------

      int jp = (t1 == 1 ? i2 : i1)*nGauss;
      fe.iGP = firstp + jp; // Global integration point counter

      for (int i = 0; i < nGauss && ok; i++, fe.iGP++)
      {
	// Gauss point coordinates along the edge
	xi[t1-1] = edgeDir < 0 ? -1.0 : 1.0;
	xi[t2-1] = xg[i];

	// Compute the basis functions and their derivatives, using
	// tensor product of one-dimensional Lagrange polynomials
        for (size_t b = 0; b < nxx.size(); ++b)
          if (!Lagrange::computeBasis(fe.basis(b+1),dNxdu[b],
                                      elem_sizes[b][0],xi[0],
                                      elem_sizes[b][1],xi[1]))
            ok = false;

        // Compute basis function derivatives and the edge normal
        fe.detJxW = utl::Jacobian(Jac,normal,fe.grad(itgBasis),Xnod,
                                  dNxdu[itgBasis-1],t1,t2);
        if (fe.detJxW == 0.0) continue; // skip singular points

        for (size_t b = 0; b < nxx.size(); ++b)
          if (b != (size_t)itgBasis-1)
            fe.grad(b+1).multiply(dNxdu[b],Jac);

	if (edgeDir < 0) normal *= -1.0;

	// Cartesian coordinates of current integration point
    X.assign(Xnod * fe.basis(itgBasis));

	// Evaluate the integrand and accumulate element contributions
	fe.detJxW *= wg[i];
	if (ok && !integrand.evalBouMx(*A,fe,time,X,normal))
	  ok = false;
      }

      // Finalize the element quantities
      if (ok && !integrand.finalizeElementBou(*A,fe,time))
        ok = false;

      // Assembly of global system integral
      if (ok && !glInt.assemble(A->ref(),fe.iel))
	ok = false;

      A->destruct();

      if (!ok) return false;
    }

  return true;
}


bool ASMs2DmxLag::evalSolution (Matrix& sField, const Vector& locSol,
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


bool ASMs2DmxLag::evalSolution (Matrix& sField, const IntegrandBase& integrand,
				const RealArray*, bool) const
{
  sField.resize(0,0);

  double incx = 2.0/double(p1-1);
  double incy = 2.0/double(p2-1);

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
    IntVec::const_iterator f2start = itgBasis == 1? MNPC[iel-1].begin() :
                                     MNPC[iel-1].begin() +
                                     std::accumulate(elem_size.begin()+itgBasis-2,
                                                     elem_size.begin()+itgBasis-1, 0);
    IntVec::const_iterator f2end = f2start + elem_size[itgBasis-1];
    IntVec mnpc1(f2start,f2end);

    this->getElementCoordinates(Xnod,iel);

    int i, j, loc = 0;
    for (j = 0; j < p2; j++)
      for (i = 0; i < p1; i++, loc++)
      {
	double xi  = -1.0 + i*incx;
	double eta = -1.0 + j*incy;
        for (size_t b = 0; b < nxx.size(); ++b)
          if (!Lagrange::computeBasis(fe.basis(b+1),dNxdu[b],
                                      elem_sizes[b][0],xi,
                                      elem_sizes[b][1],eta))
	  return false;

        // Compute Jacobian inverse of the coordinate mapping and
        // basis function derivatives w.r.t. Cartesian coordinates
        if (!fe.Jacobian(Jac,Xnod,itgBasis,nullptr,&dNxdu))
          continue; // skip singular points

	// Now evaluate the solution field
    if (!integrand.evalSol(solPt,fe,Xnod*fe.basis(itgBasis),
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


ASMs2DmxLag::BasisFunctionCache::BasisFunctionCache (const ASMs2DLag& pch,
                                                     ASM::CachePolicy plcy,
                                                     int b) :
  ASMs2DLag::BasisFunctionCache(pch,plcy,b)
{
}


ASMs2DmxLag::BasisFunctionCache::BasisFunctionCache (const BasisFunctionCache& cache,
                                                     int b) :
  ASMs2DLag::BasisFunctionCache(cache,b)
{
}


BasisFunctionVals ASMs2DmxLag::BasisFunctionCache::calculatePt (size_t el,
                                                                size_t gp,
                                                                bool reduced) const
{
  std::array<size_t,2> gpIdx = this->gpIndex(gp,reduced);
  const Quadrature& q = reduced ? *reducedQ : *mainQ;

  const ASMs2DmxLag& pch = static_cast<const ASMs2DmxLag&>(patch);

  BasisFunctionVals result;
  if (nderiv == 1)
    Lagrange::computeBasis(result.N,result.dNdu,
                           pch.elem_sizes[basis-1][0],q.xg[0][gpIdx[0]],
                           pch.elem_sizes[basis-1][1],q.xg[1][gpIdx[1]]);

  return result;
}
