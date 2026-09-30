// $Id$
//==============================================================================
//!
//! \file ASMs3Drecovery.C
//!
//! \date Mar 06 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Recovery of secondary solutions for structured 3D spline FE models.
//!
//==============================================================================

#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/VolumeInterpolator.h"

#include "ASMs3D.h"
#include "ASMmxBase.h"
#include "Function.h"
#include "GoTools/geometry/SplineSurface.h"
#include "ItgPoint.h"
#include "Field.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "GlbL2projector.h"
#include "SparseMatrix.h"
#include "SplineUtils.h"
#include "Vec3Oper.h"
#include "Profiler.h"
#include <array>


bool ASMs3D::getGrevilleParameters (RealArray& prm, int dir, int basisNum) const
{
  if (dir < 0 || dir > 2) return false;

  const Go::BsplineBasis& basis = this->getBasis(basisNum)->basis(dir);

  prm.resize(basis.numCoefs());
  for (size_t i = 0; i < prm.size(); i++)
    prm[i] = basis.grevilleParameter(i);

  return true;
}


bool ASMs3D::getQuasiInterplParameters (RealArray& prm, int dir) const
{
  if (!svol || dir < 0 || dir > 2) return false;

  const Go::BsplineBasis& basis = svol->basis(dir);

  RealArray knots_simple;
  basis.knotsSimple(knots_simple);

  prm.clear();
  for (size_t i = 0; i+1 < knots_simple.size(); i++)
  {
    prm.push_back(knots_simple[i]);
    prm.push_back(0.5*(knots_simple[i]+knots_simple[i+1]));
  }
  prm.push_back(knots_simple.back());

  return true;
}


bool ASMs3D::evaluate (const ASMbase* basis, const Vector& locVec,
                       RealArray& vec, int basisNum) const
{
  const ASMs3D* pch = dynamic_cast<const ASMs3D*>(basis);
  if (!pch) return false;

  // Compute parameter values of the result sampling points (Greville points)
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir,basisNum))
      return false;

  // Evaluate the result field at all sampling points.
  // Note: it is here assumed that *basis and *this have spline bases
  // defined over the same parameter domain.
  Matrix sValues;
  if (!pch->evalSolution(sValues,locVec,gpar.data()))
    return false;

  const Go::SplineVolume* svol = this->getBasis(basisNum);

  // Project the results onto the spline basis to find control point
  // values based on the result values evaluated at the Greville points.
  // Note that we here implicitly assume that the number of Greville points
  // equals the number of control points such that we don't have to resize
  // the result array. Think that is always the case, but beware if trying
  // other projection schemes later.

  RealArray weights;
  if (svol->rational())
    svol->getWeights(weights);

  Go::SplineVolume* vol_new =
    Go::VolumeInterpolator::regularInterpolation(svol->basis(0),
                                                 svol->basis(1),
                                                 svol->basis(2),
                                                 gpar[0], gpar[1], gpar[2],
                                                 sValues, sValues.rows(),
                                                 svol->rational(), weights);

  vec.assign(vol_new->coefs_begin(),vol_new->coefs_end());
  delete vol_new;

  return true;
}


Go::SplineVolume* ASMs3D::projectSolution (const IntegrandBase& integrand) const
{
  PROFILE2("ASMs3D::projectSolution");

  const int basis = 1;

  // Compute parameter values of the result sampling points (Greville points)
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir,basis))
      return nullptr;

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar.data()) || sValues.rows() == 0)
    return nullptr;

  const Go::SplineVolume* pvol = this->getBasis(basis);

  // Project the results onto the spline basis to find control point
  // values based on the result values evaluated at the Greville points.
  // Note that we here implicitly assume that the number of Greville points
  // equals the number of control points such that we don't have to resize
  // the result array. Think that is always the case, but beware if trying
  // other projection schemes later.

  RealArray weights;
  if (pvol->rational())
    pvol->getWeights(weights);

  return Go::VolumeInterpolator::regularInterpolation(pvol->basis(0),
                                                      pvol->basis(1),
                                                      pvol->basis(2),
                                                      gpar[0], gpar[1], gpar[2],
                                                      sValues, sValues.rows(),
                                                      pvol->rational(),
                                                      weights);
}


Go::GeomObject* ASMs3D::evalSolution (const IntegrandBase& integrand) const
{
  return this->projectSolution(integrand);
}


bool ASMs3D::assembleL2matrices (SystemMatrix& A, SystemVector& B,
                                 const L2Integrand& integrand,
                                 bool continuous) const
{
  const Go::SplineVolume* geo  = this->getBasis(ASM::GEOMETRY_BASIS);
  const Go::SplineVolume* proj = this->getBasis(ASM::PROJECTION_BASIS);
  const bool separateProjBasis = proj != geo;
  const bool singleBasis = !separateProjBasis && this->getNoBasis() == 1;
  const bool checkAge = this->getElementActivator() != nullptr;
  if (separateProjBasis && checkAge)
  {
    std::cerr <<" *** ASMs3D::assembleL2matrices: Separate projection basis "
              <<" can not be used with element activation."<< std::endl;
    return false;
  }

  const int p1 = proj->order(0);
  const int p2 = proj->order(1);
  const int p3 = proj->order(2);
  const int n1 = proj->numCoefs(0);
  const int n2 = proj->numCoefs(1);
  const int nel1 = proj->numCoefs(0) - p1 + 1;
  const int nel2 = proj->numCoefs(1) - p2 + 1;

  int pmax = p1 > p2 ? p1 : p2;
  if (pmax < p3) pmax = p3;

  // Get Gaussian quadrature point coordinates (and weights if continuous)
  const int ng1 = continuous ? this->getNoGaussPt(pmax,true) : p1 - 1;
  const int ng2 = continuous ? ng1 : p2 - 1;
  const int ng3 = continuous ? ng2 : p3 - 1;
  const double* xg = GaussQuadrature::getCoord(ng1);
  const double* yg = GaussQuadrature::getCoord(ng2);
  const double* zg = GaussQuadrature::getCoord(ng3);
  const double* wg = continuous ? GaussQuadrature::getWeight(ng1) : 0;
  if (!xg || !yg || !zg) return false;
  if (continuous && !wg) return false;

  // Compute parameter values of the Gauss points over the whole patch
  std::array<RealArray,3> gpar;
  SplineUtils::getGaussParameters(gpar[0],ng1,xg,proj->basis(0));
  SplineUtils::getGaussParameters(gpar[1],ng2,yg,proj->basis(1));
  SplineUtils::getGaussParameters(gpar[2],ng3,zg,proj->basis(2));

  // Evaluate basis functions at all integration points
  std::vector<Go::BasisPts>    spl1;
  std::vector<Go::BasisDerivs> spl2;
  if (continuous)
    geo->computeBasisGrid(gpar[0],gpar[1],gpar[2],spl2);

  if (!continuous || separateProjBasis)
    proj->computeBasisGrid(gpar[0],gpar[1],gpar[2],spl1);

  // Evaluate the secondary solution at all integration points
  Matrix sField;
  if (!integrand.evaluate(sField,gpar.data()))
  {
    std::cerr <<" *** ASMs3D::assembleL2matrices: Failed for patch "<< idx+1
              <<" nPoints="<< gpar[0].size()*gpar[1].size()*gpar[2].size()
              << std::endl;
    return false;
  }

  const size_t npnod = this->getNoProjectionNodes();
  const IntMat gmnpc = this->getElmNodes(ASM::PROJECTION_BASIS);
  A.preAssemble(gmnpc, gmnpc.size());


  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t g = 0; g < projThreadGroups.size() && ok; g++)
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < projThreadGroups[g].size(); t++) {
      double dV = 1.0;
      Vector phi(p1*p2*p3);
      Matrix dNdu, Xnod, J;
      for (int iel : projThreadGroups[g][t])
      {
        if (checkAge && this->inActiveElement(iel,integrand.getTimeLevel()))
          continue; // inactive element

        const IntVec& mnpc = gmnpc[iel];
        if (mnpc.empty())
          continue;

        int i1 = iel % nel1;
        int i2 = (iel / nel1) % nel2;
        int i3 = iel / (nel1*nel2);
        int ip = ((i3*ng2*nel2 + i2)*ng1*nel1 + i1)*ng3;

        if (continuous)
        {
          // Set up control point (nodal) coordinates for current element
          bool skip = false;
          if (singleBasis)
          {
            if (!this->getElementCoordinates(Xnod,1+iel))
              skip = true;
            if ((dV = 0.125*this->getParametricVolume(1+iel)) < 0.0)
              skip = true; // topology error (probably logic error)
          }
          else
          {
            if (!this->getElementCoordinatesPrm(Xnod,
                                                gpar[0][i1*ng1],
                                                gpar[1][i2*ng2],
                                                gpar[2][i3*ng3]))
              skip = true;
            if ((dV = 0.125 * proj->knotSpan(0, mnpc.back() % n1)
                            * proj->knotSpan(1,(mnpc.back() / n1) % n2)
                            * proj->knotSpan(2, mnpc.back() / (n1*n2))) < 0.0)
              skip = true;
          }
          if (skip)
          {
            ok = false;
            continue;
          }
        }

        // --- Integration loop over all Gauss points in each direction --------

        Matrix eA(p1*p2*p3, p1*p2*p3);
        Vectors eB(sField.rows(), Vector(p1*p2*p3));
        for (int k = 0; k < ng3; k++, ip += ng2*(nel2-1)*ng1*nel1)
          for (int j = 0; j < ng2; j++, ip += ng1*(nel1-1))
            for (int i = 0; i < ng1; i++, ip++)
            {
              if (continuous)
                SplineUtils::extractBasis(spl2[ip],phi,dNdu);

              if (!continuous || separateProjBasis)
                phi = spl1[ip].basisValues;

              // Compute the Jacobian inverse and derivatives
              double dJw = dV;
              if (continuous)
              {
                dJw *= wg[i]*wg[j]*wg[k]*utl::Jacobian(J,dNdu,Xnod,dNdu,false);
                if (dJw == 0.0) continue; // skip singular points
              }

              // Integrate the mass matrix
              eA.outer_product(phi, phi, true, dJw);

              // Integrate the rhs vector B
              for (size_t r = 1; r <= sField.rows(); r++)
                eB[r-1].add(phi,sField(r,ip+1)*dJw);
            }

        A.assemble(eA, mnpc);
        B.assemble(eB, mnpc, npnod);
      }
    }

  return ok;
}


#include "ASMs3DInterpolate.C" // TODO: inline these methods instead...


/*!
  \note A Variation Diminishing Spline Approximation is used here as the
  regular interpolation method in GoTools only works with uniform knots.
*/

bool ASMs3D::evaluate (const Field* field, RealArray& vec, int basisNum) const
{
  // Compute parameter values of the result sampling points (Greville points)
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir,basisNum))
      return false;

  // Evaluate the result field at all sampling points.
  // Note: it is here assumed that *basis and *this have spline bases
  // defined over the same parameter domain.
  Vector sValues;
  sValues.reserve(gpar[0].size()*gpar[1].size()*gpar[2].size());
  for (double w : gpar[2])
    for (double v : gpar[1])
      for (double u : gpar[0])
        sValues.push_back(field->valueFE(ItgPoint(u,v,w)));

  // Project the results onto the spline basis to find control point
  // values based on the result values evaluated at the Greville points.
  // Note that we here implicitly assume that the number of Greville points
  // equals the number of control points such that we don't have to resize
  // the result array. Think that is always the case, but beware if trying
  // other projection schemes later.

  Go::SplineVolume* vol_new =
    VariationDiminishingSplineApproximation(this->getBasis(basisNum),sValues,1);

  vec.assign(vol_new->coefs_begin(),vol_new->coefs_end());
  delete vol_new;

  return true;
}


Go::SplineVolume* ASMs3D::projectSolutionLeastSquare (const IntegrandBase& integrand) const
{
  if (!svol) return nullptr;

  // Get Gaussian quadrature points and weights
  int p = svol->order(0);
  for (int d = 1; d < 3; d++)
    if (p < svol->order(d)) p = svol->order(d);
  const int ng = this->getNoGaussPt(p,true);
  const double* xg = GaussQuadrature::getCoord(ng);
  const double* wg = GaussQuadrature::getWeight(ng);
  if (!xg || !wg) return nullptr;

  std::array<RealArray,3> gpar, wgpar;
  for (int dir = 0; dir < 3; dir++)
  {
    // Get Gauss parameter values and associated weights
    const Go::BsplineBasis& basis = svol->basis(dir);
    SplineUtils::getGaussParameters(gpar[dir],ng,xg,basis);
    RealArray::const_iterator knotit = basis.begin();
    RealArray& tmp = wgpar[dir];
    tmp.reserve(ng*(basis.numCoefs()-basis.order()));
    for (int i = basis.order(); i <= basis.numCoefs(); i++)
    {
      double d = knotit[i] - knotit[i-1];
      for (int j = 0; j < ng; j++)
        tmp.push_back(d > 0.0 ? wg[j]*d*0.5 : 0.0);
    }
  }

  // Evaluate the secondary solution at all sampling points (Gauss points)
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar.data()))
    return nullptr;

  RealArray weights;
  if (svol->rational())
    svol->getWeights(weights);

  return leastsquare_approximation(svol->basis(0),
                                   svol->basis(1),
                                   svol->basis(2),
                                   gpar[0], gpar[1], gpar[2],
                                   wgpar[0], wgpar[1], wgpar[2],
                                   sValues,
                                   sValues.rows(),
                                   svol->rational(),
                                   weights);
}


Go::SplineVolume* ASMs3D::projectSolutionLocal (const IntegrandBase& integrand) const
{
  // Compute parameter values of the result sampling points (Greville points)
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 2; dir++)
    if (!this->getQuasiInterplParameters(gpar[dir],dir))
      return nullptr;

  for (int dir = 2; dir < 3; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir))
      return nullptr;

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar.data()))
    return nullptr;

  RealArray weights;
  if (svol->rational())
    svol->getWeights(weights);

  return quasiInterpolation(svol->basis(0),
                            svol->basis(1),
                            svol->basis(2),
                            gpar[0], gpar[1], gpar[2],
                            sValues,
                            sValues.rows(),
                            svol->rational(),
                            weights);
}


Go::SplineVolume* ASMs3D::projectSolutionLocalApprox(const IntegrandBase& integrand) const
{
  // Compute parameter values of the result sampling points (Greville points)
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir))
      return nullptr;

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar.data()))
    return nullptr;

  // Project onto the geometry basis
  return VariationDiminishingSplineApproximation(svol.get(),sValues,sValues.rows());
}


/*!
  On a Piola mapped basis the prescribed value lives in the physical frame
  whereas the constrained degree of freedom is a coefficient of the reference
  basis, so the value has to be pulled back before it is fitted over the
  boundary face. On a div-compatible basis a prescribed velocity vector is
  pulled back by the adjugate of the jacobian, and a prescribed normal
  velocity by the dilation of the face, since the contravariant Piola
  transform preserves the normal flux. On a curl-compatible basis a
  prescribed vector field is pulled back by the transposed jacobian, as the
  covariant Piola transform preserves the circulation.
*/

Go::SplineSurface* ASMs3D::projectPiolaDirichlet (const DirichletFace& dface,
                                                  const RealFunc* sf,
                                                  const VecFunc* vf,
                                                  double time,
                                                  bool tangent) const
{
  const Go::SplineVolume* geo = this->getBasis(ASM::GEOMETRY_BASIS);
  if (!geo || !dface.surf)
    return nullptr;

  const int ndir = abs(dface.dir); // parameter direction of the normal
  const bool covariant = ASMmxBase::Type == ASMmxBase::CURL_COMPATIBLE;
  // The degrees of freedom of a compatible basis are the components of the
  // reference field, one component per basis, so the basis being constrained
  // says which component the prescribed value has to be pulled back to.
  const int comp = dface.basis;
  if (sf && (covariant || comp != ndir))
  {
    std::cerr <<" *** ASMs3D::projectPiolaDirichlet: Only the normal velocity"
              <<" on a div-compatible basis can be prescribed as a scalar on a"
              <<" face, not on basis "<< comp <<". Give the whole field."
              << std::endl;
    return nullptr;
  }
  // The reference normal points along the parameter direction on the upper
  // boundary and against it on the lower one
  const double nsign = dface.dir > 0 ? 1.0 : -1.0;

  // The parameter directions of the face surface in the volume
  const int s1 = ndir == 1 ? 1 : 0;
  const int s2 = ndir == 3 ? 1 : 2;

  double par[3] = { 0.0, 0.0, 0.0 };
  par[ndir-1] = dface.dir < 0 ? geo->startparam(ndir-1)
                              : geo->endparam(ndir-1);

  // Fit it in the parametric least-squares sense over the face
  const Go::BsplineBasis& bu = dface.surf->basis_u();
  const Go::BsplineBasis& bv = dface.surf->basis_v();
  const int nu = bu.numCoefs(), nv = bv.numCoefs();
  const int pu = bu.order(), pv = bv.order();
  const RealArray ku(bu.begin(),bu.end()), kv(bv.begin(),bv.end());

  const int nGP = this->getNoGaussPt(std::max(pu,pv),true);
  const double* xg = GaussQuadrature::getCoord(nGP);
  const double* wg = GaussQuadrature::getWeight(nGP);
  if (!xg || !wg) return nullptr;

  const int ncoef = nu*nv;

  // The weights of a rational face, whose basis the fitted coefficients
  // are the degrees of freedom of
  RealArray weights;
  if (dface.surf->rational())
  {
    const int dim = dface.surf->dimension();
    weights.reserve(ncoef);
    for (auto it = dface.surf->rcoefs_begin(); it != dface.surf->rcoefs_end(); it += dim+1)
      weights.push_back(*(it+dim));
  }

  SparseMatrix A(SparseMatrix::SUPERLU);
  StdVector B(ncoef);
  A.resize(ncoef,ncoef);
  RealArray Nu(pu), Nv(pv), R(pu*pv);
  std::vector<Go::Point> pts(4);

  for (size_t js = pv-1; js+pv < kv.size(); js++)
  {
    if (kv[js+1] <= kv[js]) continue; // skip the repeated knots
    const double dtv = 0.5*(kv[js+1]-kv[js]), tmv = 0.5*(kv[js+1]+kv[js]);
    for (size_t is = pu-1; is+pu < ku.size(); is++)
    {
      if (ku[is+1] <= ku[is]) continue;
      const double dtu = 0.5*(ku[is+1]-ku[is]), tmu = 0.5*(ku[is+1]+ku[is]);
      for (int jg = 0; jg < nGP; jg++)
        for (int ig = 0; ig < nGP; ig++)
        {
          double tu = tmu + dtu*xg[ig]; // not const, knotIntervalFuzzy may snap
          double tv = tmv + dtv*xg[jg];
          par[s1] = tu;
          par[s2] = tv;
          geo->point(pts,par[0],par[1],par[2],1);

          const Vec3 X = SplineUtils::toVec3(pts[0],nsd);
          const std::array<Vec3,3> dX{ SplineUtils::toVec3(pts[1],nsd),
                                       SplineUtils::toVec3(pts[2],nsd),
                                       SplineUtils::toVec3(pts[3],nsd) };
          const Vec4 Xt(X,time,par);

          double value = 0.0;
          if (sf)
          {
            // The prescribed value is the normal velocity itself, and the
            // Piola transform preserves the normal flux, u*n*dS = uhat*nhat*dShat
            const double dS = Vec3(dX[s1],dX[s2]).length();
            if (dS <= 0.0)
            {
              std::cerr <<" *** ASMs3D::projectPiolaDirichlet: Degenerate"
                        <<" boundary at parameters "<< tu <<" "<< tv
                        << std::endl;
              return nullptr;
            }
            value = nsign * dS * (tangent ? sf->timeDerivative(Xt) : (*sf)(Xt));
          }
          else if (vf)
          {
            const Vec3 u = tangent ? vf->timeDerivative(Xt) : (*vf)(Xt);
            if (covariant)
              // The component along the parameter direction of this basis
              value = dX[comp-1]*u;
            else
              // A row of the adjugate of the jacobian is the cross product of
              // the other two columns
              value = Vec3(dX[comp%3],dX[(comp+1)%3])*u;
          }
          else
            return nullptr;

          const int kiu = bu.knotIntervalFuzzy(tu);
          const int kiv = bv.knotIntervalFuzzy(tv);
          bu.computeBasisValues(tu,Nu.data(),0);
          bv.computeBasisValues(tv,Nv.data(),0);

          // The basis values of the face, rational if it is
          double W = 0.0;
          for (int b = 0; b < pv; b++)
            for (int a = 0; a < pu; a++)
            {
              const int ia = kiu-pu+1+a + nu*(kiv-pv+1+b);
              double& Ra = R[a+pu*b];
              Ra = Nu[a]*Nv[b];
              if (!weights.empty() && ia >= 0 && ia < ncoef)
                Ra *= weights[ia];
              W += Ra;
            }
          if (!weights.empty() && W > 0.0)
            for (double& Ra : R)
              Ra /= W;

          const double w = dtu*dtv*wg[ig]*wg[jg];
          for (int b = 0; b < pv; b++)
            for (int a = 0; a < pu; a++)
            {
              const int ia = kiu-pu+1+a + nu*(kiv-pv+1+b);
              if (ia < 0 || ia >= ncoef) continue;
              const double Ra = R[a+pu*b];
              B(1+ia) += Ra*value*w;
              for (int d = 0; d < pv; d++)
                for (int c = 0; c < pu; c++)
                {
                  const int ib = kiu-pu+1+c + nu*(kiv-pv+1+d);
                  if (ib >= 0 && ib < ncoef)
                    A(1+ia,1+ib) += Ra*R[c+pu*d]*w;
                }
            }
        }
    }
  }

  if (!A.solve(B))
  {
    std::cerr <<" *** ASMs3D::projectPiolaDirichlet: Failed to solve the"
              <<" least-squares system."<< std::endl;
    return nullptr;
  }

  return new Go::SplineSurface(nu,nv,pu,pv,ku.begin(),kv.begin(),B.begin(),1);
}
