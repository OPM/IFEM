// $Id$
//==============================================================================
//!
//! \file ASMmxBase.C
//!
//! \date Dec 28 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for spline-based mixed finite element assembly drivers.
//!
//==============================================================================

#include "ASMmxBase.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SurfaceInterpolator.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/VolumeInterpolator.h"
#include "SplineUtils.h"
#include <array>
#include <numeric>


char ASMmxBase::itgBasis             = 2;
ASMmxBase::MixedType ASMmxBase::Type = ASMmxBase::FULL_CONT_RAISE_BASIS1;


void ASMmxBase::initMx (const std::vector<int>& MLGN, const int* sysMadof)
{
  MADOF.clear();
  MADOF.reserve(MLGN.size());
  for (int n : MLGN)
    MADOF.push_back(sysMadof[n-1]-1);
}


void ASMmxBase::extractNodeVecMx (const RealArray& glbVec, RealArray& nodVec,
                                  int basis) const
{
  int b0 = basis-1, b1 = basis;
  size_t ofs = 0, len = 0;
  if (basis < 1 || basis > (int)nfx.size())
  {
    b0 = 0;
    b1 = nfx.size();
    for (int i = b0; i < b1; i++)
      len += nfx[i]*nb[i];
  }
  else
  {
    len = nfx[basis-1]*nb[basis-1];
    if (basis > 1)
      ofs = std::accumulate(nb.begin(), nb.begin()+basis-1, 0);
  }

  nodVec.resize(len);

  int ldof = 0;
  for (int b = b0; b < b1; ofs += nb[b++])
    for (size_t i = ofs; i < nb[b]+ofs; i++)
    {
      int idof = MADOF[i];
      for (size_t j = 0; j < nfx[b]; j++)
        nodVec[ldof++] = glbVec[idof++];
    }
}


void ASMmxBase::injectNodeVecMx (RealArray& glbVec, const RealArray& nodVec,
                                 int basis) const
{
  int b0 = basis-1, b1 = basis;
  size_t ofs = 0;
  if (basis < 1 || basis > (int)nfx.size())
  {
    b0 = 0;
    b1 = nfx.size();
  }
  else if (basis > 1)
    ofs = std::accumulate(nb.begin(), nb.begin()+basis-1, 0);

  int ldof = 0;
  for (int b = b0; b < b1; b++)
    for (size_t i = ofs; i < nb[b]+ofs; i++)
    {
      int idof = MADOF[i];
      for (size_t j = 0; j < nfx[b]; j++)
        glbVec[idof++] = nodVec[ldof++];
    }
}


bool ASMmxBase::getSolutionMx (Matrix& sField, const Vector& locSol,
                               const std::vector<int>& nodes) const
{
  if (nodes.empty()) return true;

  int nvar, low = 1, high = nb.front();
  if (nodes.front() <= high)
    nvar = nfx[0];
  else
  {
    nvar = nfx[1];
    low += nb[0];
    high += nb[1];
  }

  sField.resize(nvar,nodes.size());
  for (size_t i = 0; i < nodes.size(); i++)
    if (nodes[i] < low || nodes[i] > high)
    {
      std::cerr <<" *** ASMmxBase::getSolutionMx: Node #"<< nodes[i]
                <<" is out of range ["<< low <<","<< high <<"]."<< std::endl;
      return false;
    }
    else
    {
      int idof = nvar*(nodes[i]-1);
      if (low > 1) idof += nfx.front()*nb.front();
      for (int j = 0; j < nvar; j++)
        sField(j+1,i+1) = locSol[idof++];
    }

  return true;
}


ASMmxBase::SurfaceVec ASMmxBase::establishBases (Go::SplineSurface* surf,
                                                 MixedType type)
{
  SurfaceVec result(2);
  // With mixed methods we need two separate spline spaces
  if (type == FULL_CONT_RAISE_BASIS1 || type == FULL_CONT_RAISE_BASIS2)
  {
    // basis1 should be one degree higher than basis2 and C^p-1 continuous
    result[0].reset(ASMmxBase::raiseBasis(surf));
    result[1].reset(new Go::SplineSurface(*surf));
    itgBasis = 2;
  }
  else if (type == REDUCED_CONT_RAISE_BASIS1 || type == REDUCED_CONT_RAISE_BASIS2)
  {
    // Order-elevate basis1 such that it is of one degree higher than basis2
    // but only C^p-2 continuous
    result[0].reset(new Go::SplineSurface(*surf));
    result[0]->raiseOrder(1,1);
    result[1].reset(new Go::SplineSurface(*surf));
    itgBasis = 2;
  }
  else if (ASMmxBase::Type == ASMmxBase::DIV_COMPATIBLE)
  {
    result.resize(3);

    // basis1 should be one degree higher than basis2 and C^p-1 continuous
    int ndim = surf->dimension();
    Go::BsplineBasis a1 = surf->basis(0);
    Go::BsplineBasis a2 = surf->basis(1);
    Go::BsplineBasis b1 = surf->basis(0).extendedBasis(surf->order_u()+1);
    Go::BsplineBasis b2 = surf->basis(1).extendedBasis(surf->order_v()+1);

    // Compute parameter values of the Greville points
    size_t i;
    RealArray u0(a1.numCoefs()), v0(a2.numCoefs());
    for (i = 0; i < u0.size(); i++)
      u0[i] = a1.grevilleParameter(i);
    for (i = 0; i < v0.size(); i++)
      v0[i] = a2.grevilleParameter(i);
    RealArray ug(b1.numCoefs()), vg(b2.numCoefs());
    for (i = 0; i < ug.size(); i++)
      ug[i] = b1.grevilleParameter(i);
    for (i = 0; i < vg.size(); i++)
      vg[i] = b2.grevilleParameter(i);

    // Evaluate the spline surface at all points
    // Project the coordinates onto the new basis (the 2nd XYZ is dummy here)
    RealArray XYZ0(ndim*ug.size()*v0.size()), XYZ1(ndim*u0.size()*vg.size());
    surf->gridEvaluator(XYZ0,ug,v0);
    surf->gridEvaluator(XYZ1,u0,vg);
    result[2].reset(new Go::SplineSurface(*surf));
    result[0].reset(Go::SurfaceInterpolator::regularInterpolation(b1,a2,
                                                                  ug,v0,XYZ0,ndim,
                                                                  false,XYZ0));
    result[1].reset(Go::SurfaceInterpolator::regularInterpolation(a1,b2,
                                                                  u0,vg,XYZ1,ndim,
                                                                  false,XYZ1));
    itgBasis = 3;
  } else if (type == SUBGRID) {
    // basis1 should be one degree higher than basis2 and C^p-1 continuous
    int ndim = surf->dimension();
    result[1].reset(new Go::SplineSurface(*surf));
    Go::SplineSurface tmp(*surf);
    for (size_t i = 0; i < 2; ++i) {
      RealArray extraKnots;
      RealArray::const_iterator uit = tmp.basis(i).begin();
      double uprev = *(uit++);
      while (uit != tmp.basis(i).end())
      {
        double ucurr = *(uit++);
        if (ucurr > uprev)
          extraKnots.push_back(ucurr*0.5 + uprev*0.5);
        uprev = ucurr;
      }
      if (i == 0)
        tmp.insertKnot_u(extraKnots);
      else
        tmp.insertKnot_v(extraKnots);
    }
    Go::BsplineBasis b1 = tmp.basis(0).extendedBasis(tmp.order_u()+1);
    Go::BsplineBasis b2 = tmp.basis(1).extendedBasis(tmp.order_v()+1);

    // Compute parameter values of the Greville points
    size_t i;
    RealArray ug(b1.numCoefs()), vg(b2.numCoefs());
    for (i = 0; i < ug.size(); i++)
      ug[i] = b1.grevilleParameter(i);
    for (i = 0; i < vg.size(); i++)
      vg[i] = b2.grevilleParameter(i);

    if (surf->rational()) {
      std::vector<double> rCoefs(tmp.rcoefs_begin(), tmp.rcoefs_end());

      // we normally would set coefs as (x*w, y*w, w)
      // however, gotools use this representation internally already.

      // instance a Bspline surface in ndim+1
      Go::SplineSurface surf2(tmp.basis(0), tmp.basis(1), rCoefs.begin(), ndim+1, false);

      // interpolate the Bspline surface onto new basis
      RealArray XYZ((ndim+1)*ug.size()*vg.size());
      surf2.gridEvaluator(XYZ,ug,vg);
      std::unique_ptr<Go::SplineSurface> surf3(Go::SurfaceInterpolator::regularInterpolation(b1,b2,ug,vg,XYZ,ndim+1,false,XYZ));

      // new rational coefs are (x/w', y/w', w')
      // apparently gotools will rescale coeffs on surface creation.
      result[0].reset(new Go::SplineSurface(surf3->basis(0), surf3->basis(1), surf3->coefs_begin(), ndim, true));
    } else {
      RealArray XYZ(ndim*ug.size()*vg.size());
      // Evaluate the spline surface at all points
      tmp.gridEvaluator(XYZ,ug,vg);

      // Project the coordinates onto the new basis (the 2nd XYZ is dummy here)
      result[0].reset(Go::SurfaceInterpolator::regularInterpolation(b1,b2,
                                                                    ug,vg,XYZ,ndim,
                                                                    false,XYZ));
    }
    itgBasis = 1;
  }

  if (type == FULL_CONT_RAISE_BASIS2 || type == REDUCED_CONT_RAISE_BASIS2)
    std::swap(result[0], result[1]), itgBasis = 1;

  return result;
}


ASMmxBase::VolumeVec ASMmxBase::establishBases(Go::SplineVolume* svol,
                                               MixedType type)
{
  VolumeVec result(2);
  // With mixed methods we need two separate spline spaces
  if (type == FULL_CONT_RAISE_BASIS1 || type == FULL_CONT_RAISE_BASIS2)
  {
    // basis1 should be one degree higher than basis2 and C^p-1 continuous
    result[0].reset(ASMmxBase::raiseBasis(svol));
    result[1].reset(new Go::SplineVolume(*svol));
    itgBasis = 2;
  }
  else if (type == REDUCED_CONT_RAISE_BASIS1 || type == REDUCED_CONT_RAISE_BASIS2)
  {
    // Order-elevate basis1 such that it is of one degree higher than basis2
    // but only C^p-2 continuous
    result[0].reset(new Go::SplineVolume(*svol));
    result[0]->raiseOrder(1,1,1);
    result[1].reset(new Go::SplineVolume(*svol));
    itgBasis = 2;
  }
  else if (ASMmxBase::Type == ASMmxBase::DIV_COMPATIBLE)
  {
    result.resize(4);

    // basis1 should be one degree higher than basis2 and C^p-1 continuous
    int ndim = svol->dimension();
    Go::BsplineBasis a1 = svol->basis(0);
    Go::BsplineBasis a2 = svol->basis(1);
    Go::BsplineBasis a3 = svol->basis(2);
    Go::BsplineBasis b1 = svol->basis(0).extendedBasis(svol->order(0)+1);
    Go::BsplineBasis b2 = svol->basis(1).extendedBasis(svol->order(1)+1);
    Go::BsplineBasis b3 = svol->basis(2).extendedBasis(svol->order(2)+1);

    // Compute parameter values of the Greville points
    size_t i;
    RealArray u0(a1.numCoefs()), v0(a2.numCoefs()), w0(a3.numCoefs());
    for (i = 0; i < u0.size(); i++)
      u0[i] = a1.grevilleParameter(i);
    for (i = 0; i < v0.size(); i++)
      v0[i] = a2.grevilleParameter(i);
    for (i = 0; i < w0.size(); i++)
      w0[i] = a3.grevilleParameter(i);
    RealArray ug(b1.numCoefs()), vg(b2.numCoefs()), wg(b3.numCoefs());
    for (i = 0; i < ug.size(); i++)
      ug[i] = b1.grevilleParameter(i);
    for (i = 0; i < vg.size(); i++)
      vg[i] = b2.grevilleParameter(i);
    for (i = 0; i < wg.size(); i++)
      wg[i] = b3.grevilleParameter(i);

    // Evaluate the spline surface at all points
    // Project the coordinates onto the new basis (the 2nd XYZ is dummy here)
    RealArray XYZ0(ndim*ug.size()*v0.size()*w0.size()), XYZ1(ndim*u0.size()*vg.size()*w0.size()), XYZ2(ndim*u0.size()*v0.size()*wg.size());
    svol->gridEvaluator(ug,v0,w0,XYZ0);
    svol->gridEvaluator(u0,vg,w0,XYZ1);
    svol->gridEvaluator(u0,v0,wg,XYZ2);
    result[0].reset(Go::VolumeInterpolator::regularInterpolation(b1,a2,a3,
                                                                 ug,v0,w0,XYZ0,ndim,
                                                                 false,XYZ0));
    result[1].reset(Go::VolumeInterpolator::regularInterpolation(a1,b2,a3,
                                                                 u0,vg,w0,XYZ1,ndim,
                                                                 false,XYZ1));
    result[2].reset(Go::VolumeInterpolator::regularInterpolation(a1,a2,b3,
                                                                 u0,v0,wg,XYZ2,ndim,
                                                                 false,XYZ2));
    result[3].reset(new Go::SplineVolume(*svol));
    itgBasis = 4;
  } else if (type == SUBGRID) {
    // basis1 should be one degree higher than basis2 and C^p-1 continuous
    int ndim = svol->dimension();
    result[1].reset(new Go::SplineVolume(*svol));
    Go::SplineVolume tmp(*svol);
    for (size_t dir = 0; dir < 3; ++dir) {
      RealArray extraKnots;
      RealArray::const_iterator uit = tmp.basis(dir).begin();
      double uprev = *(uit++);
      while (uit != tmp.basis(dir).end())
      {
        double ucurr = *(uit++);
        if (ucurr > uprev)
          extraKnots.push_back(ucurr*0.5 + uprev*0.5);
        uprev = ucurr;
      }
      tmp.insertKnot(dir, extraKnots);
    }
    Go::BsplineBasis b1 = tmp.basis(0).extendedBasis(tmp.order(0)+1);
    Go::BsplineBasis b2 = tmp.basis(1).extendedBasis(tmp.order(1)+1);
    Go::BsplineBasis b3 = tmp.basis(2).extendedBasis(tmp.order(2)+1);

    // Compute parameter values of the Greville points
    size_t i;
    RealArray ug(b1.numCoefs()), vg(b2.numCoefs()), wg(b3.numCoefs());
    for (i = 0; i < ug.size(); i++)
      ug[i] = b1.grevilleParameter(i);
    for (i = 0; i < vg.size(); i++)
      vg[i] = b2.grevilleParameter(i);
    for (i = 0; i < wg.size(); i++)
      wg[i] = b3.grevilleParameter(i);

    if (svol->rational()) {
      std::vector<double> rCoefs(tmp.rcoefs_begin(), tmp.rcoefs_end());

      // we normally would set coefs as (x*w, y*w, w)
      // however, gotools use this representation internally already.

      // instance a Bspline surface in ndim+1
      Go::SplineVolume svol2(tmp.basis(0), tmp.basis(1), tmp.basis(2), rCoefs.begin(), ndim+1, false);

      // interpolate the Bspline surface onto new basis
      RealArray XYZ((ndim+1)*ug.size()*vg.size()*wg.size());
      svol2.gridEvaluator(ug,vg,wg,XYZ);
      std::unique_ptr<Go::SplineVolume> svol3(Go::VolumeInterpolator::regularInterpolation(b1,b2,b3,ug,vg,wg,XYZ,ndim+1,false,XYZ));

      // new rational coefs are (x/w', y/w', w')
      // apparently gotools will rescale coeffs on surface creation.
      result[0].reset(new Go::SplineVolume(svol3->basis(0), svol3->basis(1), svol3->basis(2),svol3->coefs_begin(), ndim, true));
    } else {
      RealArray XYZ(ndim*ug.size()*vg.size()*wg.size());
      // Evaluate the spline surface at all points
      tmp.gridEvaluator(ug,vg,wg,XYZ);

      // Project the coordinates onto the new basis (the 2nd XYZ is dummy here)
      result[0].reset(Go::VolumeInterpolator::regularInterpolation(b1,b2,b3,
                                                                   ug,vg,wg,
                                                                   XYZ,ndim,
                                                                   false,XYZ));
    }
    itgBasis = 1;
  }

  if (type == FULL_CONT_RAISE_BASIS2 || type == REDUCED_CONT_RAISE_BASIS2)
    std::swap(result[0], result[1]), itgBasis = 1;

  return result;
}


Go::SplineSurface* ASMmxBase::raiseBasis (Go::SplineSurface* surf)
{
  // Create a C^p-1 basis of one degree higher than *surf
  // but keep lines of reduced continuity
  std::array<Go::BsplineBasis,2> basis;
  for (size_t idx = 0; idx < 2; ++idx) {
    std::vector<double> knots;
    surf->basis(idx).knotsSimple(knots);
    std::vector<int> mult;
    surf->basis(idx).knotMultiplicities(mult);
    std::vector<int> cont(knots.size());
    cont.front() = cont.back() = -1;
    int order = idx == 0 ? surf->order_u() : surf->order_v();
    for (size_t i = 1; i < knots.size()-1; ++i)
      cont[i]  = order - (mult[i] == 1 ? 1 : mult[i]+1);

    std::vector<double> newknot = SplineUtils::buildKnotVector(order, knots, cont);
    basis[idx] = Go::BsplineBasis(order+1, newknot.begin(), newknot.end());
  }

  // Compute parameter values of the Greville points
  std::array<RealArray,2> ug;
  for (size_t idx = 0; idx < 2; ++idx) {
    ug[idx].resize(basis[idx].numCoefs());
    for (size_t i = 0; i < ug[idx].size(); i++)
      ug[idx][i] = basis[idx].grevilleParameter(i);
  }

  int ndim = surf->dimension();
  if (surf->rational())
  {
    std::vector<double> rCoefs(surf->rcoefs_begin(), surf->rcoefs_end());

    // We normally would set coefs as (x*w, y*w, w).
    // However, GoTools uses this representation internally already.

    // Instance a Bspline surface in ndim+1
    Go::SplineSurface surf2(surf->basis(0), surf->basis(1), rCoefs.begin(), ndim+1, false);

    // Interpolate the Bspline surface onto new basis
    RealArray XYZ((ndim+1)*ug[0].size()*ug[1].size());
    surf2.gridEvaluator(XYZ,ug[0],ug[1]);
    std::unique_ptr<Go::SplineSurface> surf3(Go::SurfaceInterpolator::regularInterpolation(basis[0],basis[1],ug[0],ug[1],XYZ,ndim+1,false,XYZ));

    // New rational coefs are (x/w', y/w', w')
    // Apparently, GoTools will rescale coeffs on surface creation.
    return new Go::SplineSurface(surf3->basis(0), surf3->basis(1), surf3->coefs_begin(), ndim, true);
  }

  // Evaluate the spline surface at all points
  RealArray XYZ(ndim*ug[0].size()*ug[1].size());
  surf->gridEvaluator(XYZ,ug[0],ug[1]);

  // Project the coordinates onto the new basis (the 2nd XYZ is dummy here)
  return Go::SurfaceInterpolator::regularInterpolation(basis[0],basis[1],
                                                      ug[0],ug[1],XYZ,ndim,false,XYZ);
}


Go::SplineVolume* ASMmxBase::raiseBasis (Go::SplineVolume* svol)
{
  // Create a C^p-1 basis of one degree higher than *svol
  // but keep lines of reduced continuity
  std::array<Go::BsplineBasis,3> basis;
  for (size_t idx = 0; idx < 3; ++idx) {
    std::vector<double> knots;
    svol->basis(idx).knotsSimple(knots);
    std::vector<int> mult;
    svol->basis(idx).knotMultiplicities(mult);
    std::vector<int> cont(knots.size());
    cont.front() = cont.back() = -1;
    int order = svol->order(idx);
    for (size_t i = 1; i < knots.size()-1; ++i)
      cont[i]  = order - (mult[i] == 1 ? 1 : mult[i]+1);

    std::vector<double> newknot = SplineUtils::buildKnotVector(order, knots, cont);
    basis[idx] = Go::BsplineBasis(order+1, newknot.begin(), newknot.end());
  }

  // Compute parameter values of the Greville points
  std::array<RealArray,3> ug;
  for (size_t idx = 0; idx < 3; ++idx) {
    ug[idx].resize(basis[idx].numCoefs());
    for (size_t i = 0; i < ug[idx].size(); i++)
      ug[idx][i] = basis[idx].grevilleParameter(i);
  }

  int ndim = svol->dimension();
  if (svol->rational())
  {
    std::vector<double> rCoefs(svol->rcoefs_begin(), svol->rcoefs_end());

    // We normally would set coefs as (x*w, y*w, w).
    // However, GoTools use this representation internally already.

    // Instance a Bspline surface in ndim+1
    Go::SplineVolume vol2(svol->basis(0), svol->basis(1), svol->basis(2), rCoefs.begin(), ndim+1, false);

    // Interpolate the Bspline surface onto new basis
    RealArray XYZ((ndim+1)*ug[0].size()*ug[1].size()*ug[2].size());
    vol2.gridEvaluator(XYZ,ug[0],ug[1],ug[2]);
    std::unique_ptr<Go::SplineVolume> svol3(Go::VolumeInterpolator::regularInterpolation(basis[0],basis[1],basis[2],ug[0],ug[1],ug[2],XYZ,ndim+1,false,XYZ));

    // New rational coefs are (x/w', y/w', w')
    // Apparently, GoTools will rescale coeffs on surface creation.
    return new Go::SplineVolume(svol3->basis(0), svol3->basis(1), svol3->basis(2), svol3->coefs_begin(), ndim, true);
  }

  RealArray XYZ(ndim*ug[0].size()*ug[1].size()*ug[2].size());
  // Evaluate the spline surface at all points
  svol->gridEvaluator(ug[0],ug[1],ug[2],XYZ);
  // Project the coordinates onto the new basis (the 2nd XYZ is dummy here)
  return Go::VolumeInterpolator::regularInterpolation(basis[0],basis[1],basis[2],ug[0],ug[1],ug[2],XYZ,ndim,false,XYZ);
}
