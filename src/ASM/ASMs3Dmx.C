// $Id$
//==============================================================================
//!
//! \file ASMs3Dmx.C
//!
//! \date Dec 28 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of structured 3D spline mixed FE models.
//!
//==============================================================================

#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/VolumeInterpolator.h"

#include "ASMs3Dmx.h"
#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "Utilities.h"
#include "Profiler.h"
#include "Vec3Oper.h"
#include "Vec3.h"


ASMs3Dmx::ASMs3Dmx (unsigned char n_f1, unsigned char n_f2)
  : ASMs3D(n_f1+n_f2), ASMmxBase(n_f1,n_f2)
{
  basis1 = basis2 = 0;
}


ASMs3Dmx::ASMs3Dmx (const ASMs3Dmx& patch, char n_f1, char n_f2)
  : ASMs3D(patch), ASMmxBase(patch.nf1,patch.nf2)
{
  basis1 = patch.basis1;
  basis2 = patch.basis2;
  nb1 = patch.nb1;
  nb2 = patch.nb2;
  if (n_f1 >= 0) nf1 = n_f1;
  if (n_f2 >= 0) nf2 = n_f2;
  nf = nf1 + nf2;
}


Go::SplineVolume* ASMs3Dmx::getBasis (int basis) const
{
  return basis == 2 ? basis2 : basis1;
}


Go::SplineSurface* ASMs3Dmx::getBoundary (int dir)
{
  if (dir < -3 || dir == 0 || dir > 3)
    return NULL;

  // The boundary surfaces are stored internally in the SplineVolume object
  int iface = dir > 0 ? 2*dir-1 : -2*dir-2;
  return basis1->getBoundarySurface(iface).get();
}


bool ASMs3Dmx::write (std::ostream& os, int basis) const
{
  if (basis1 && basis == 1)
    os <<"700 1 0 0\n" << *basis1;
  else if (basis2 && basis == 2)
    os <<"700 1 0 0\n" << *basis2;
  else if (svol)
    os <<"700 1 0 0\n" << *svol;
  else
    return false;

  return os.good();
}


void ASMs3Dmx::clear (bool retainGeometry)
{
  // Erase the spline data
  if (retainGeometry)
  {
    if (basis1 && basis1 != svol && !shareFE) delete basis1;
    if (basis2 && basis2 != svol && !shareFE) delete basis2;
  }
  else
  {
    if (basis1 && !shareFE) delete basis1;
    if (basis2 && !shareFE) delete basis2;
    svol = 0;
  }
  basis1 = basis2 = 0;

  // Erase the FE data
  this->ASMs3D::clear(retainGeometry);
}


size_t ASMs3Dmx::getNoNodes (int basis) const
{
  switch (basis)
    {
    case 1: return nb1;
    case 2: return nb2;
    }

  return this->ASMbase::getNoNodes(basis);
}


unsigned char ASMs3Dmx::getNoFields (int basis) const
{
  switch (basis)
    {
    case 1: return nf1;
    case 2: return nf2;
    }

  return nf;
}


unsigned char ASMs3Dmx::getNodalDOFs (size_t inod) const
{
  if (this->isLMn(inod)) return nLag;
  return inod <= nb1 || inod > nb1+nb2 ? nf1 : nf2;
}


char ASMs3Dmx::getNodeType (size_t inod) const
{
  if (this->isLMn(inod)) return 'L';
  return inod <= nb1 ? 'D' : (inod <= nb1+nb2 ? 'P' : 'X');
}


void ASMs3Dmx::initMADOF (const int* sysMadof)
{
  this->initMx(MLGN,sysMadof);
}


void ASMs3Dmx::extractNodeVec (const Vector& globRes, Vector& nodeVec,
			       unsigned char, int basis) const
{
  this->extractNodeVecMx(globRes,nodeVec,basis);
}


bool ASMs3Dmx::getSolution (Matrix& sField, const Vector& locSol,
			    const IntVec& nodes) const
{
  return this->getSolutionMx(sField,locSol,nodes);
}


bool ASMs3Dmx::generateFEMTopology ()
{
  if (!svol) return false;

  if (!basis1 && !basis2)
  {
    // With mixed methods we need two separate spline spaces
    if (useCpminus1)
    {
      // basis1 should be one degree higher than basis2 and C^p-1 continuous
      int ndim = svol->dimension();
      Go::BsplineBasis b1 = svol->basis(0).extendedBasis(svol->order(0)+1);
      Go::BsplineBasis b2 = svol->basis(1).extendedBasis(svol->order(1)+1);
      Go::BsplineBasis b3 = svol->basis(2).extendedBasis(svol->order(2)+1);
      /* To lower order and regularity this can be used instead
      std::vector<double>::const_iterator first =  ++svol->basis(0).begin();
      std::vector<double>::const_iterator last  =  --svol->basis(0).end();
      Go::BsplineBasis b1 = Go::BsplineBasis(svol->order(0)-1,first,last);
      first =  ++svol->basis(1).begin();
      last  =  --svol->basis(1).end();
      Go::BsplineBasis b2 = Go::BsplineBasis(svol->order(1)-1,first,last);
      first =  ++svol->basis(2).begin();
      last  =  --svol->basis(2).end();
      Go::BsplineBasis b3 = Go::BsplineBasis(svol->order(2)-1,first,last);
      */

      // Note: Currently this is implemented for non-rational splines only.
      // TODO: Ask the splines people how to fix this properly, that is, how
      // may we obtain the correct weights for basis1 when *svol is a NURBS?
      if (svol->rational())
	std::cout <<"WARNING: The geometry basis is rational (using NURBS).\n"
		  <<"         The basis for the unknown fields of one degree"
		  <<" higher will however be non-rational.\n"
		  <<"         This may affect accuracy.\n"<< std::endl;

      // Compute parameter values of the Greville points of the new basis
      size_t i;
      RealArray ug(b1.numCoefs()), vg(b2.numCoefs()), wg(b3.numCoefs());
      for (i = 0; i < ug.size(); i++)
	ug[i] = b1.grevilleParameter(i);
      for (i = 0; i < vg.size(); i++)
	vg[i] = b2.grevilleParameter(i);
      for (i = 0; i < wg.size(); i++)
	wg[i] = b3.grevilleParameter(i);

      // Evaluate the spline volume at all points
      RealArray XYZ(ndim*ug.size()*vg.size()*wg.size());
      svol->gridEvaluator(ug,vg,wg,XYZ);

      // Project the coordinates onto the new basis (the 2nd XYZ is dummy here)
      basis1 = Go::VolumeInterpolator::regularInterpolation(b1,b2,b3,
							    ug,vg,wg,XYZ,ndim,
							    false,XYZ);
    }
    else
    {
      // Order-elevate basis1 such that it is of one degree higher than basis2
      // but only C^p-2 continuous
      basis1 = new Go::SplineVolume(*svol);
      basis1->raiseOrder(1,1,1);
    }
    basis2 = svol;

    // Define which basis that should be used to represent the geometry
    if (geoUsesBasis1) svol = basis1;
  }

  const int n1 = basis1->numCoefs(0);
  const int n2 = basis1->numCoefs(1);
  const int n3 = basis1->numCoefs(2);
  const int m1 = basis2->numCoefs(0);
  const int m2 = basis2->numCoefs(1);
  const int m3 = basis2->numCoefs(2);
  if (!nodeInd.empty() && !shareFE)
  {
    if (nodeInd.size() == nb1 + nb2) return true;
    std::cerr <<" *** ASMs3Dmx::generateFEMTopology: Inconsistency between the"
	      <<" number of FE nodes "<< nodeInd.size()
	      <<"\n     and the number of spline coefficients "<< nb1 + nb2
	      <<" in the patch."<< std::endl;
    return false;
  }

  nb1 = n1*n2*n3; // Number of functions in first basis
  nb2 = m1*m2*m3; // Number of functions in second basis

  if (shareFE) return true;

  const int p1 = basis1->order(0);
  const int p2 = basis1->order(1);
  const int p3 = basis1->order(2);
  const int q1 = basis2->order(0);
  const int q2 = basis2->order(1);
  const int q3 = basis2->order(2);

#ifdef SP_DEBUG
  std::cout <<"numCoefs: "<< n1 <<" "<< n2 <<" "<< n3
	    <<", "<< m1 <<" "<< m2 <<" "<< m3;
  std::cout <<"\norder: "<< p1 <<" "<< p2 <<" "<< p3
	    <<", "<< q1 <<" "<< q2 <<" "<< q3;
  for (int d = 0; d < 3; d++)
  {
    std::cout <<"\nd"<< char('u'+d) <<':';
    for (int i = 0; i < basis1->numCoefs(d); i++)
      std::cout <<' '<< basis1->knotSpan(d,i);
    for (int j = 0; j < basis2->numCoefs(d); j++)
      std::cout <<' '<< basis2->knotSpan(d,j);
  }
  std::cout << std::endl;
#endif
  // Consistency checks, just to be fool-proof
  if (m1 <  2 || m2 <  2 || m3 <  2) return false;
  if (q1 <  1 || q2 <  1 || q3 <  1) return false;
  if (p1 > n1 || p2 > n2 || p3 > n3) return false;
  if (q1 > m1 || q2 > m2 || q3 > m3) return false;

  myMLGE.resize((n1-p1+1)*(n2-p2+1)*(n3-p3+1),0);
  myMLGN.resize(nb1 + nb2);
  myMNPC.resize(myMLGE.size());
  myNodeInd.resize(myMLGN.size());

  int i1, i2, i3, j1, j2, j3;
  size_t iel, inod = 0;
  for (i3 = 0; i3 < n3; i3++)
    for (i2 = 0; i2 < n2; i2++)
      for (i1 = 0; i1 < n1; i1++)
      {
	myNodeInd[inod].I = i1;
	myNodeInd[inod].J = i2;
	myNodeInd[inod].K = i3;
	myMLGN[inod++]    = ++gNod;
      }

  for (i3 = 0; i3 < m3; i3++)
    for (i2 = 0; i2 < m2; i2++)
      for (i1 = 0; i1 < m1; i1++)
      {
	myNodeInd[inod].I = i1;
	myNodeInd[inod].J = i2;
	myNodeInd[inod].K = i3;
	myMLGN[inod++]    = ++gNod;
      }

  if (geoUsesBasis1)
  {
    // Create nodal connectivities for basis 1
    iel = inod = 0;
    for (i3 = 1; i3 <= n3; i3++)
      for (i2 = 1; i2 <= n2; i2++)
	for (i1 = 1; i1 <= n1; i1++, inod++)
	  if (i1 >= p1 && i2 >= p2 && i3 >= p3)
	  {
	    if (basis1->knotSpan(0,i1-1) > 0.0)
	      if (basis1->knotSpan(1,i2-1) > 0.0)
		if (basis1->knotSpan(2,i3-1) > 0.0)
		{
		  myMLGE[iel] = ++gEl; // global element number over all patches
		  myMNPC[iel].resize(p1*p2*p3+q1*q2*q3,0);

		  int lnod = 0;
		  for (j3 = p3-1; j3 >= 0; j3--)
		    for (j2 = p2-1; j2 >= 0; j2--)
		      for (j1 = p1-1; j1 >= 0; j1--)
			myMNPC[iel][lnod++] = inod - n1*n2*j3 - n1*j2 - j1;
		}

	    iel++;
	  }

    // Create nodal connectivities for basis 2
    iel = 0;
    for (i3 = 1; i3 <= m3; i3++)
      for (i2 = 1; i2 <= m2; i2++)
	for (i1 = 1; i1 <= m1; i1++, inod++)
	  if (i1 >= q1 && i2 >= q2 && i3 >= q3)
	    if (basis2->knotSpan(0,i1-1) > 0.0)
	      if (basis2->knotSpan(1,i2-1) > 0.0)
		if (basis2->knotSpan(2,i3-1) > 0.0)
		{
		  while (iel < myMNPC.size() && myMNPC[iel].empty()) iel++;

		  int lnod = p1*p2*p3;
		  for (j3 = q3-1; j3 >= 0; j3--)
		    for (j2 = q2-1; j2 >= 0; j2--)
		      for (j1 = q1-1; j1 >= 0; j1--)
			myMNPC[iel][lnod++] = inod - m1*m2*j3 - m1*j2 - j1;

		  iel++;
		}
  }
  else
  {
    // Create nodal connectivities for basis 2
    iel = 0;
    inod = n1*n2*n3;
    for (i3 = 1; i3 <= m3; i3++)
      for (i2 = 1; i2 <= m2; i2++)
	for (i1 = 1; i1 <= m1; i1++, inod++)
	  if (i1 >= q1 && i2 >= q2 && i3 >= q3)
	  {
	    if (basis2->knotSpan(0,i1-1) > 0.0)
	      if (basis2->knotSpan(1,i2-1) > 0.0)
		if (basis2->knotSpan(2,i3-1) > 0.0)
		{
		  myMLGE[iel] = ++gEl; // global element number over all patches
		  myMNPC[iel].resize(p1*p2*p3+q1*q2*q3,0);

		  int lnod = p1*p2*p3;
		  for (j3 = q3-1; j3 >= 0; j3--)
		    for (j2 = q2-1; j2 >= 0; j2--)
		      for (j1 = q1-1; j1 >= 0; j1--)
			myMNPC[iel][lnod++] = inod - m1*m2*j3 - m1*j2 - j1;
		}

	    iel++;
	  }

    // Create nodal connectivities for basis 1
    iel = inod = 0;
    for (i3 = 1; i3 <= n3; i3++)
      for (i2 = 1; i2 <= n2; i2++)
	for (i1 = 1; i1 <= n1; i1++, inod++)
	  if (i1 >= p1 && i2 >= p2 && i3 >= p3)
	    if (basis1->knotSpan(0,i1-1) > 0.0)
	      if (basis1->knotSpan(1,i2-1) > 0.0)
		if (basis1->knotSpan(2,i3-1) > 0.0)
		{
		  while (iel < myMNPC.size() && myMNPC[iel].empty()) iel++;

		  int lnod = 0;
		  for (j3 = p3-1; j3 >= 0; j3--)
		    for (j2 = p2-1; j2 >= 0; j2--)
		      for (j1 = p1-1; j1 >= 0; j1--)
			myMNPC[iel][lnod++] = inod - n1*n2*j3 - n1*j2 - j1;

		  iel++;
		}
  }

#ifdef SP_DEBUG
  std::cout <<"NEL = "<< MLGE.size() <<" NNOD = "<< MLGN.size() << std::endl;
#endif
  return true;
}


bool ASMs3Dmx::connectPatch (int face, ASMs3D& neighbor, int nface, int norient)
{
  ASMs3Dmx* neighMx = dynamic_cast<ASMs3Dmx*>(&neighbor);
  if (!neighMx) return false;

  if (swapW && face > 4) // Account for swapped parameter direction
    face = 11-face;

  if (neighMx->swapW && face > 4) // Account for swapped parameter direction
    nface = 11-nface;

  return this->connectBasis(face,neighbor,nface,norient,1,0,0)
    &&   this->connectBasis(face,neighbor,nface,norient,2,nb1,neighMx->nb1);
}


void ASMs3Dmx::closeFaces (int dir, int, int)
{
  this->ASMs3D::closeFaces(dir,1,1);
  this->ASMs3D::closeFaces(dir,2,nb1+1);
}


bool ASMs3Dmx::getElementCoordinates (Matrix& X, int iel) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMs3Dmx::getElementCoordinates: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return false;
  }
#endif

  size_t nenod = svol->order(0)*svol->order(1)*svol->order(2);
  size_t lnod0 = basis1->order(0)*basis1->order(1)*basis1->order(2);
  if (geoUsesBasis1) lnod0 = 0;

  X.resize(3,nenod);
  const IntVec& mnpc = MNPC[iel-1];

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int nd = svol->dimension();
  RealArray::const_iterator cit = svol->coefs_begin();
  for (size_t n = 0; n < nenod; n++)
  {
    int iI = nodeInd[mnpc[lnod0+n]].I;
    int iJ = nodeInd[mnpc[lnod0+n]].J;
    int iK = nodeInd[mnpc[lnod0+n]].K;
    int ip = ((iK*n2 + iJ)*n1 + iI)*nd;
    for (size_t i = 0; i < 3; i++)
      X(i+1,n+1) = *(cit+(ip+i));
  }

#if SP_DEBUG > 2
  std::cout <<"\nCoordinates for element "<< iel << X << std::endl;
#endif
  return true;
}


Vec3 ASMs3Dmx::getCoord (size_t inod) const
{
  if (inod > nodeInd.size() && inod <= MLGN.size())
  {
    // This is a node added due to constraints in local directions.
    // Find the corresponding original node (see constrainFaceLocal)
    std::map<size_t,size_t>::const_iterator it = xnMap.find(inod);
    if (it != xnMap.end()) inod = it->second;
  }

#ifdef INDEX_CHECK
  if (inod < 1 || inod > nodeInd.size())
  {
    std::cerr <<" *** ASMs3Dmx::getCoord: Node index "<< inod
	      <<" out of range [1,"<< nodeInd.size() <<"]."<< std::endl;
    return Vec3();
  }
#endif

  RealArray::const_iterator cit;
  const int I = nodeInd[inod-1].I;
  const int J = nodeInd[inod-1].J;
  const int K = nodeInd[inod-1].K;
  if (inod <= nb1)
    cit = basis1->coefs_begin()
      + ((K*basis1->numCoefs(1)+J)*basis1->numCoefs(0)+I) * basis1->dimension();
  else
    cit = basis2->coefs_begin()
      + ((K*basis2->numCoefs(1)+J)*basis2->numCoefs(0)+I) * basis2->dimension();

  return Vec3(*cit,*(cit+1),*(cit+2));
}


bool ASMs3Dmx::getSize (int& n1, int& n2, int& n3, int basis) const
{
  switch (basis)
    {
    case 1:
      if (!basis1) return false;
      n1 = basis1->numCoefs(0);
      n2 = basis1->numCoefs(1);
      n3 = basis1->numCoefs(2);
      return true;

    case 2:
      if (!basis2) return false;
      n1 = basis2->numCoefs(0);
      n2 = basis2->numCoefs(1);
      n3 = basis2->numCoefs(2);
      return true;
    }

  return this->ASMs3D::getSize(n1,n2,n3);
}


bool ASMs3Dmx::integrate (Integrand& integrand,
			  GlobalIntegral& glInt,
			  const TimeDomain& time)
{
  if (!svol) return true; // silently ignore empty patches
  if (!basis1 || !basis2) return false;

  PROFILE2("ASMs3Dmx::integrate(I)");

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Compute parameter values of the Gauss points over the whole patch
  Matrix gpar[3];
  for (int d = 0; d < 3; d++)
    this->getGaussPointParameters(gpar[d],d,nGauss,xg);

  // Evaluate basis function derivatives at all integration points
  std::vector<Go::BasisDerivs> spline1, spline2;
  basis1->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline1);
  basis2->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline2);

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);

  const int nel1 = n1 - p1 + 1;
  const int nel2 = n2 - p2 + 1;
  const int nel3 = n3 - p3 + 1;


  // === Assembly loop over all elements in the patch ==========================

  bool ok=true;
  for (size_t g=0;g<threadGroupsVol.size() && ok;++g) {
#pragma omp parallel for schedule(static)
    for (size_t t=0;t<threadGroupsVol[g].size();++t) {
      MxFiniteElement fe(basis1->order(0)*basis1->order(1)*basis1->order(2),
                         basis2->order(0)*basis2->order(1)*basis2->order(2));
      Matrix dN1du, dN2du, Xnod, Jac;
      Vec4   X;
      for (size_t l = 0; l < threadGroupsVol[g][t].size() && ok; ++l)
      {
        int iel = threadGroupsVol[g][t][l];
        fe.iel = MLGE[iel];
        if (fe.iel < 1) continue; // zero-volume element

        int i1 = p1 + iel % nel1;
        int i2 = p2 + (iel / nel1) % nel3;
        int i3 = p3 + iel / (nel1*nel2);

        // Get element volume in the parameter space
        double dV = this->getParametricVolume(++iel);
        if (dV < 0.0) // topology error (probably logic error)
	{
          ok = false;
          break;
        }

        // Set up control point (nodal) coordinates for current element
        if (!this->getElementCoordinates(Xnod,iel))
        {
          ok = false;
          break;
        }

        // Initialize element quantities
        IntVec::const_iterator f2start = MNPC[iel-1].begin() + fe.N1.size();
        LocalIntegral* A = integrand.getLocalIntegral(fe.N1.size(),fe.N2.size(),
                                                      fe.iel,false);
        if (!integrand.initElement(IntVec(MNPC[iel-1].begin(),f2start),
                                   IntVec(f2start,MNPC[iel-1].end()),nb1,*A))
        {
          A->destruct();
          ok = false;
          break;
        }


        // --- Integration loop over all Gauss points in each direction --------

        int ip = (((i3-p3)*nGauss*nel2 + i2-p2)*nGauss*nel1 + i1-p1)*nGauss;
        int jp = (((i3-p3)*nel2 + i2-p2*nel1 + i1-p1))*nGauss*nGauss*nGauss;
        fe.iGP = firstIp + jp; // Global integration point counter

        for (int k = 0; k < nGauss; k++, ip += nGauss*(nel2-1)*nGauss*nel1)
          for (int j = 0; j < nGauss; j++, ip += nGauss*(nel1-1))
            for (int i = 0; i < nGauss; i++, ip++, fe.iGP++)
            {
              // Local element coordinates of current integration point
              fe.xi   = xg[i];
              fe.eta  = xg[j];
              fe.zeta = xg[k];

              // Parameter values of current integration point
              fe.u = gpar[0](i+1,i1-p1+1);
              fe.v = gpar[1](j+1,i2-p2+1);
              fe.w = gpar[2](k+1,i3-p3+1);

              // Fetch basis function derivatives at current integration point
              extractBasis(spline1[ip],fe.N1,dN1du);
              extractBasis(spline2[ip],fe.N2,dN2du);

              // Compute Jacobian inverse of the coordinate mapping and
              // basis function derivatives w.r.t. Cartesian coordinates
              if (geoUsesBasis1)
              {
                fe.detJxW = utl::Jacobian(Jac,fe.dN1dX,Xnod,dN1du);
                fe.dN2dX.multiply(dN2du,Jac); // dN2dX = dN2du * J^-1
              }
              else
              {
                fe.detJxW = utl::Jacobian(Jac,fe.dN2dX,Xnod,dN2du);
                fe.dN1dX.multiply(dN1du,Jac); // dN1dX = dN1du * J^-1
              }
              if (fe.detJxW == 0.0) continue; // skip singular points

              // Cartesian coordinates of current integration point
              X = Xnod * (geoUsesBasis1 ? fe.N1 : fe.N2);
              X.t = time.t;

              // Evaluate the integrand and accumulate element contributions
              fe.detJxW *= 0.125*dV*wg[i]*wg[j]*wg[k];
              if (!integrand.evalIntMx(*A,fe,time,X))
                ok = false;
            }

        // Finalize the element quantities
        if (ok && !integrand.finalizeElement(*A,time,firstIp+jp))
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


bool ASMs3Dmx::integrate (Integrand& integrand, int lIndex,
			  GlobalIntegral& glInt,
			  const TimeDomain& time)
{
  if (!svol) return true; // silently ignore empty patches
  if (!basis1 || !basis2) return false;

  PROFILE2("ASMs3Dmx::integrate(B)");

  std::map<char,ThreadGroups>::const_iterator tit;
  if ((tit = threadGroupsFace.find(lIndex)) == threadGroupsFace.end())
  {
    std::cerr <<" *** ASMs3D::integrate: No thread groups for face "<< lIndex
	      << std::endl;
    return false;
  }
  const ThreadGroups& threadGrp = tit->second;

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Find the parametric direction of the face normal {-3,-2,-1, 1, 2, 3}
  const int faceDir = (lIndex+1)/(lIndex%2 ? -2 : 2);

  const int t1 = 1 + abs(faceDir)%3; // first tangent direction
  const int t2 = 1 + t1%3;           // second tangent direction

  // Compute parameter values of the Gauss points over the whole patch face
  Matrix gpar[3];
  for (int d = 0; d < 3; d++)
    if (-1-d == faceDir)
    {
      gpar[d].resize(1,1);
      gpar[d].fill(svol->startparam(d));
    }
    else if (1+d == faceDir)
    {
      gpar[d].resize(1,1);
      gpar[d].fill(svol->endparam(d));
    }
    else
      this->getGaussPointParameters(gpar[d],d,nGauss,xg);

  // Evaluate basis function derivatives at all integration points
  std::vector<Go::BasisDerivs> spline1, spline2;
  basis1->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline1);
  basis2->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline2);

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  const int nel1 = n1 - p1 + 1;
  const int nel2 = n2 - p2 + 1;

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;


  // === Assembly loop over all elements on the patch face =====================

  bool ok = true;
  for (size_t g = 0; g < threadGrp.size() && ok; ++g) {
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < threadGrp[g].size(); ++t) {
      MxFiniteElement fe(basis1->order(0)*basis1->order(1)*basis1->order(2),
                         basis2->order(0)*basis2->order(1)*basis2->order(2));
      fe.xi = fe.eta = fe.zeta = faceDir < 0 ? -1.0 : 1.0;
      fe.u = gpar[0](1,1);
      fe.v = gpar[1](1,1);
      fe.w = gpar[2](1,1);

      Matrix dN1du, dN2du, Xnod, Jac;
      Vec4   X;
      Vec3   normal;
      for (size_t l = 0; l < threadGrp[g][t].size() && ok; ++l)
      {
        int iel = threadGrp[g][t][l];
        fe.iel = MLGE[iel];
        if (fe.iel < 1) continue; // zero-volume element

        int i1 = p1 + iel % nel1;
        int i2 = p2 + (iel / nel1) % nel2;
        int i3 = p3 + iel / (nel1*nel2);

	// Get element face area in the parameter space
	double dA = this->getParametricArea(++iel,abs(faceDir));
	if (dA < 0.0) // topology error (probably logic error)
	{
          ok = false;
          break;
        }

	// Set up control point coordinates for current element
	if (!this->getElementCoordinates(Xnod,iel))
	{
          ok = false;
          break;
        }

	// Initialize element quantities
	IntVec::const_iterator f2start = MNPC[iel-1].begin() + fe.N1.size();
        LocalIntegral* A = integrand.getLocalIntegral(fe.N1.size(),fe.N2.size(),
                                                      fe.iel,true);
	if (!integrand.initElementBou(IntVec(MNPC[iel-1].begin(),f2start),
				      IntVec(f2start,MNPC[iel-1].end()),nb1,*A))
        {
          A->destruct();
          ok = false;
          break;
        }

        // Define some loop control variables depending on which face we are on
        int nf1, j1, j2;
        switch (abs(faceDir))
        {
          case 1: nf1 = nel2; j2 = i3-p3; j1 = i2-p2; break;
          case 2: nf1 = nel1; j2 = i3-p3; j1 = i1-p1; break;
          case 3: nf1 = nel1; j2 = i2-p2; j1 = i1-p1; break;
          default: nf1 = j1 = j2 = 0;
        }


	// --- Integration loop over all Gauss points in each direction --------

        int k1, k2, k3;
        int ip = (j2*nGauss*nf1 + j1)*nGauss;
        int jp = (j2*nf1 + j1)*nGauss*nGauss;
        fe.iGP = firstp + jp; // Global integration point counter

        for (int j = 0; j < nGauss; j++, ip += nGauss*(nf1-1))
          for (int i = 0; i < nGauss; i++, ip++, fe.iGP++)
          {
            // Local element coordinates and parameter values
            // of current integration point
            switch (abs(faceDir))
            {
              case 1: k2 = i+1; k3 = j+1; k1 = 0; break;
              case 2: k1 = i+1; k3 = j+1; k2 = 0; break;
              case 3: k1 = i+1; k2 = j+1; k3 = 0; break;
              default: k1 = k2 = k3 = 0;
            }
            if (gpar[0].size() > 1)
            {
              fe.xi = xg[k1];
              fe.u = gpar[0](k1,i1-p1+1);
            }
            if (gpar[1].size() > 1)
            {
              fe.eta = xg[k2];
              fe.v = gpar[1](k2,i2-p2+1);
            }
            if (gpar[2].size() > 1)
            {
              fe.zeta = xg[k3];
              fe.w = gpar[2](k3,i3-p3+1);
            }

            // Fetch basis function derivatives at current integration point
            extractBasis(spline1[ip],fe.N1,dN1du);
            extractBasis(spline2[ip],fe.N2,dN2du);

            // Compute Jacobian inverse of the coordinate mapping and
            // basis function derivatives w.r.t. Cartesian coordinates
            if (geoUsesBasis1)
            {
              fe.detJxW = utl::Jacobian(Jac,normal,fe.dN1dX,Xnod,dN1du,t1,t2);
              fe.dN2dX.multiply(dN2du,Jac); // dN2dX = dN2du * J^-1
            }
            else
            {
              fe.detJxW = utl::Jacobian(Jac,normal,fe.dN2dX,Xnod,dN2du,t1,t2);
              fe.dN1dX.multiply(dN1du,Jac); // dN1dX = dN1du * J^-1
            }
            if (fe.detJxW == 0.0) continue; // skip singular points

            if (faceDir < 0) normal *= -1.0;

            // Cartesian coordinates of current integration point
            X = Xnod * (geoUsesBasis1 ? fe.N1 : fe.N2);
            X.t = time.t;

            // Evaluate the integrand and accumulate element contributions
            fe.detJxW *= 0.25*dA*wg[i]*wg[j];
            if (!integrand.evalBouMx(*A,fe,time,X,normal))
              ok = false;
          }

	// Assembly of global system integral
	if (ok && !glInt.assemble(A->ref(),fe.iel))
          ok = false;

        A->destruct();
      }
    }
  }

  return ok;
}


int ASMs3Dmx::evalPoint (const double* xi, double* param, Vec3& X) const
{
  if (!svol) return -3;

  int i;
  for (i = 0; i < 3; i++)
    param[i] = (1.0-xi[i])*svol->startparam(i) + xi[i]*svol->endparam(i);

  Go::Point X0;
  svol->point(X0,param[0],param[1],param[2]);
  for (i = 0; i < 3 && i < svol->dimension(); i++)
    X[i] = X0[i];

  // Check if this point matches any of the basis1 control points (nodes)
  Vec3 Xnod;
  size_t inod = 1;
  RealArray::const_iterator cit = basis1->coefs_begin();
  for (i = 0; cit != basis1->coefs_end(); cit++, i++)
  {
    if (i < 3) Xnod[i] = *cit;
    if (i+1 == basis1->dimension())
      if (X.equal(Xnod,0.001))
	return inod;
      else
      {
	inod++;
	i = -1;
      }
  }

  return 0;
}


bool ASMs3Dmx::evalSolution (Matrix& sField, const Vector& locSol,
			     const RealArray* gpar, bool regular) const
{
  if (!basis1 || !basis2) return false;

  // Evaluate the basis functions at all points
  std::vector<Go::BasisPts> spline1, spline2;
  if (regular)
  {
    basis1->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline1);
    basis2->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline2);
  }
  else if (gpar[0].size() == gpar[1].size() && gpar[0].size() == gpar[2].size())
  {
    spline1.resize(gpar[0].size());
    spline2.resize(gpar[0].size());
    for (size_t i = 0; i < spline1.size(); i++)
    {
      basis1->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline1[i]);
      basis2->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline2[i]);
    }
  }
  else
    return false;

  const int p1 = basis1->order(0);
  const int p2 = basis1->order(1);
  const int p3 = basis1->order(2);
  const int n1 = basis1->numCoefs(0);
  const int n2 = basis1->numCoefs(1);
  const int n3 = basis1->numCoefs(2);

  const int q1 = basis2->order(0);
  const int q2 = basis2->order(1);
  const int q3 = basis2->order(2);
  const int m1 = basis2->numCoefs(0);
  const int m2 = basis2->numCoefs(1);
  const int m3 = basis2->numCoefs(2);

  size_t nc1 = nf1;
  size_t nc2 = 0;
  if (nc1*nb1 < locSol.size())
    nc2 = (locSol.size() - nc1*nb1)/nb2;
  else
    nc1 = locSol.size()/nb1;

  if (nc1*nb1 + nc2*nb2 != locSol.size())
    return false;

  Matrix Xtmp;
  Vector Ytmp, Ztmp;

  // Evaluate the primary solution field at each point
  size_t nPoints = spline1.size();
  sField.resize(nc1+nc2,nPoints);
  for (size_t i = 0; i < nPoints; i++)
  {
    IntVec ip;
    scatterInd(n1,n2,n3,p1,p2,p3,spline1[i].left_idx,ip);

    utl::gather(ip,nc1,locSol,Xtmp);
    Xtmp.multiply(spline1[i].basisValues,Ytmp);

    if (nc2 > 0)
    {
      ip.clear();
      scatterInd(m1,m2,m3,q1,q2,q3,spline2[i].left_idx,ip);

      utl::gather(ip,nc2,locSol,Xtmp,nc1*nb1);
      Xtmp.multiply(spline2[i].basisValues,Ztmp);

      Ytmp.insert(Ytmp.end(),Ztmp.begin(),Ztmp.end());
    }
    sField.fillColumn(1+i,Ytmp);
  }

  return true;
}


bool ASMs3Dmx::evalSolution (Matrix& sField, const IntegrandBase& integrand,
			     const RealArray* gpar, bool regular) const
{
  sField.resize(0,0);

  if (!basis1 || !basis2) return false;

  // Evaluate the basis functions and their derivatives at all points
  std::vector<Go::BasisDerivs> spline1, spline2;
  if (regular)
  {
    basis1->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline1);
    basis2->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline2);
  }
  else if (gpar[0].size() == gpar[1].size() && gpar[0].size() == gpar[2].size())
  {
    spline1.resize(gpar[0].size());
    spline2.resize(gpar[0].size());
    for (size_t i = 0; i < spline1.size(); i++)
    {
      basis1->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline1[i]);
      basis2->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline2[i]);
    }
  }

  const int p1 = basis1->order(0);
  const int p2 = basis1->order(1);
  const int p3 = basis1->order(2);
  const int n1 = basis1->numCoefs(0);
  const int n2 = basis1->numCoefs(1);
  const int n3 = basis1->numCoefs(2);

  const int q1 = basis2->order(0);
  const int q2 = basis2->order(1);
  const int q3 = basis2->order(2);
  const int m1 = basis2->numCoefs(0);
  const int m2 = basis2->numCoefs(1);
  const int m3 = basis2->numCoefs(2);

  // Fetch nodal (control point) coordinates
  Matrix Xnod, Xtmp;
  this->getNodalCoordinates(Xnod);

  Vector N1(p1*p2*p3), N2(q1*q2*q3), solPt;
  Matrix dN1du, dN1dX, dN2du, dN2dX, Jac;
  Vec3   X;

  // Evaluate the secondary solution field at each point
  size_t nPoints = spline1.size();
  for (size_t i = 0; i < nPoints; i++)
  {
    // Fetch indices of the non-zero basis functions at this point
    IntVec ip1, ip2;
    scatterInd(n1,n2,n3,p1,p2,p3,spline1[i].left_idx,ip1);
    scatterInd(m1,m2,m3,q1,q2,q3,spline2[i].left_idx,ip2);

    // Fetch associated control point coordinates
    utl::gather(geoUsesBasis1 ? ip1 : ip2, 3, Xnod, Xtmp);

    // Fetch basis function derivatives at current integration point
    extractBasis(spline1[i],N1,dN1du);
    extractBasis(spline2[i],N2,dN2du);

    // Compute Jacobian inverse of the coordinate mapping and
    // basis function derivatives w.r.t. Cartesian coordinates
    if (geoUsesBasis1)
      if (utl::Jacobian(Jac,dN1dX,Xtmp,dN1du) == 0.0) // Jac = (Xtmp * dN1du)^-1
	continue; // skip singular points
      else
	dN2dX.multiply(dN2du,Jac); // dN2dX = dN2du * J^-1
    else
      if (utl::Jacobian(Jac,dN2dX,Xtmp,dN2du) == 0.0) // Jac = (Xtmp * dN2du)^-1
	continue; // skip singular points
      else
	dN1dX.multiply(dN1du,Jac); // dN1dX = dN1du * J^-1

    // Cartesian coordinates of current integration point
    X = Xtmp * (geoUsesBasis1 ? N1 : N2);

    // Now evaluate the solution field
    if (!integrand.evalSol(solPt,N1,N2,dN1dX,dN2dX,X,ip1,ip2))
      return false;
    else if (sField.empty())
      sField.resize(solPt.size(),nPoints,true);

    sField.fillColumn(1+i,solPt);
  }

  return true;
}
