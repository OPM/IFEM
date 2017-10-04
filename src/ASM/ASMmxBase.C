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
#include <numeric>

char ASMmxBase::geoBasis            = 2;
ASMmxBase::MixedType ASMmxBase::Type = ASMmxBase::FULL_CONT_RAISE_BASIS1;


ASMmxBase::ASMmxBase (const std::vector<unsigned char>& n_f)
{
  nfx = n_f;
}


void ASMmxBase::initMx (const std::vector<int>& MLGN, const int* sysMadof)
{
  MADOF.resize(MLGN.size());
  for (size_t i = 0; i < MADOF.size(); i++)
    MADOF[i] = sysMadof[MLGN[i]-1]-1;
}


void ASMmxBase::extractNodeVecMx (const Vector& globRes, Vector& nodeVec,
				  int basis) const
{
  if (basis > (int)nfx.size())
    basis = 0;

  size_t len=0;
  if (basis == 0)
    for (size_t i=0;i<nfx.size();++i)
      len += nfx[i]*nb[i];
  else
    len = nfx[basis-1]*nb[basis-1];

  nodeVec.resize(len);

  size_t i, j;
  int idof, ldof = 0;
  size_t k=basis==0?1:basis;
  size_t ofs = std::accumulate(nb.begin(), nb.begin()+k-1, 0);
  for (; k < (basis==0?nfx.size()+1:(size_t)basis+1); ++k) {
    for (i = ofs; i < nb[k-1]+ofs; i++)
    {
      idof = MADOF[i];
      for (j = 0; j < nfx[k-1]; j++, ldof++)
	nodeVec[ldof] = globRes[idof++];
    }
    ofs += nb[k-1];
  }
}


void ASMmxBase::injectNodeVecMx (Vector& globRes, const Vector& nodeVec,
                                 int basis) const
{
  if (basis > (int)nfx.size())
    basis = 0;

  size_t i, j;
  int idof, ldof = 0;
  size_t k=basis==0?1:basis;
  size_t ofs = std::accumulate(nb.begin(), nb.begin()+k-1, 0);
  for (; k < (basis==0?nfx.size()+1:(size_t)basis+1); ++k) {
    for (i = ofs; i < nb[k-1]+ofs; i++)
    {
      idof = MADOF[i];
      for (j = 0; j < nfx[k-1]; j++, ldof++)
	globRes[idof++] = nodeVec[ldof];
    }
  }
}


bool ASMmxBase::getSolutionMx (Matrix& sField, const Vector& locSol,
			       const std::vector<int>& nodes) const
{
  if (nodes.empty()) return true;

  int low, high, nvar;
  if ((size_t)nodes.front() <= nb[0])
  {
    nvar = nfx[0];
    low  = 1;
    high = nb[0];
  }
  else
  {
    nvar = nfx[1];
    low  = nb[0]+1;
    high = nb[0]+nb[1];
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
      if (low > 1) idof += nfx[0]*nb[0];
      for (int j = 0; j < nvar; j++)
	sField(j+1,i+1) = locSol[idof++];
    }

  return true;
}


ASMmxBase::SurfaceVec ASMmxBase::establishBases(Go::SplineSurface* surf,
                                                MixedType type)
{
  SurfaceVec result(2);
  // With mixed methods we need two separate spline spaces
  if (type == FULL_CONT_RAISE_BASIS1 || type == FULL_CONT_RAISE_BASIS2)
  {
    // basis1 should be one degree higher than basis2 and C^p-1 continuous
    int ndim = surf->dimension();
    Go::BsplineBasis b1 = surf->basis(0).extendedBasis(surf->order_u()+1);
    Go::BsplineBasis b2 = surf->basis(1).extendedBasis(surf->order_v()+1);
    /* To lower order and regularity this can be used instead
    std::vector<double>::const_iterator first = ++surf->basis(0).begin();
    std::vector<double>::const_iterator last  = --surf->basis(0).end();
    Go::BsplineBasis b1 = Go::BsplineBasis(surf->order_u()-1,first,last);
    first =  ++surf->basis(1).begin();
    last  =  --surf->basis(1).end();
    Go::BsplineBasis b2 = Go::BsplineBasis(surf->order_v()-1,first,last);
    */

    // Compute parameter values of the Greville points
    size_t i;
    RealArray ug(b1.numCoefs()), vg(b2.numCoefs());
    for (i = 0; i < ug.size(); i++)
      ug[i] = b1.grevilleParameter(i);
    for (i = 0; i < vg.size(); i++)
      vg[i] = b2.grevilleParameter(i);

    if (surf->rational()) {
      std::vector<double> rCoefs(surf->rcoefs_begin(), surf->rcoefs_end());

      // we normally would set coefs as (x*w, y*w, w)
      // however, gotools use this representation internally already.

      // instance a Bspline surface in ndim+1
      Go::SplineSurface surf2(surf->basis(0), surf->basis(1), rCoefs.begin(), ndim+1, false);

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
      surf->gridEvaluator(XYZ,ug,vg);

      // Project the coordinates onto the new basis (the 2nd XYZ is dummy here)
      result[0].reset(Go::SurfaceInterpolator::regularInterpolation(b1,b2,
                                                                    ug,vg,XYZ,ndim,
                                                                    false,XYZ));
    }
    result[1].reset(new Go::SplineSurface(*surf));
  }
  else if (type == REDUCED_CONT_RAISE_BASIS1 || type == REDUCED_CONT_RAISE_BASIS2)
  {
    // Order-elevate basis1 such that it is of one degree higher than basis2
    // but only C^p-2 continuous
    result[0].reset(new Go::SplineSurface(*surf));
    result[0]->raiseOrder(1,1);
    result[1].reset(new Go::SplineSurface(*surf));
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
    geoBasis = 3;
  } else if (type == SUBGRID) {
    // basis1 should be one degree higher than basis2 and C^p-1 continuous
    int ndim = surf->dimension();
    result[1].reset(new Go::SplineSurface(*surf));
    Go::SplineSurface tmp(*surf);
    for (size_t i = 0; i < 2; ++i) {
      RealArray extraKnots;
      RealArray::const_iterator uit = tmp.basis(i).begin();
      double ucurr, uprev = *(uit++);
      while (uit != tmp.basis(i).end())
      {
        ucurr = *(uit++);
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
    geoBasis = 1;
  }

  if (type == FULL_CONT_RAISE_BASIS2 || type == REDUCED_CONT_RAISE_BASIS2)
    std::swap(result[0], result[1]);

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
    int ndim = svol->dimension();
    Go::BsplineBasis b1 = svol->basis(0).extendedBasis(svol->order(0)+1);
    Go::BsplineBasis b2 = svol->basis(1).extendedBasis(svol->order(1)+1);
    Go::BsplineBasis b3 = svol->basis(2).extendedBasis(svol->order(2)+1);
    /* To lower order and regularity this can be used instead
    std::vector<double>::const_iterator first = ++surf->basis(0).begin();
    std::vector<double>::const_iterator last  = --surf->basis(0).end();
    Go::BsplineBasis b1 = Go::BsplineBasis(surf->order_u()-1,first,last);
    first =  ++surf->basis(1).begin();
    last  =  --surf->basis(1).end();
    Go::BsplineBasis b2 = Go::BsplineBasis(surf->order_v()-1,first,last);
    */

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
      std::vector<double> rCoefs(svol->rcoefs_begin(), svol->rcoefs_end());

      // we normally would set coefs as (x*w, y*w, w)
      // however, gotools use this representation internally already.

      // instance a Bspline surface in ndim+1
      Go::SplineVolume vol2(svol->basis(0), svol->basis(1), svol->basis(2), rCoefs.begin(), ndim+1, false);

      // interpolate the Bspline surface onto new basis
      RealArray XYZ((ndim+1)*ug.size()*vg.size()*wg.size());
      vol2.gridEvaluator(XYZ,ug,vg,wg);
      std::unique_ptr<Go::SplineVolume> svol3(Go::VolumeInterpolator::regularInterpolation(b1,b2,b3,ug,vg,wg,XYZ,ndim+1,false,XYZ));

      // new rational coefs are (x/w', y/w', w')
      // apparently gotools will rescale coeffs on surface creation.
      result[0].reset(new Go::SplineVolume(svol3->basis(0), svol3->basis(1), svol3->basis(2), svol3->coefs_begin(), ndim, true));
    } else {
      RealArray XYZ(ndim*ug.size()*vg.size()*wg.size());
      // Evaluate the spline surface at all points
      svol->gridEvaluator(ug,vg,wg,XYZ);
      // Project the coordinates onto the new basis (the 2nd XYZ is dummy here)
      result[0].reset(Go::VolumeInterpolator::regularInterpolation(b1,b2,b3,
                                                                   ug,vg,wg,XYZ,ndim,
                                                                   false,XYZ));
    }
    result[1].reset(new Go::SplineVolume(*svol));
  }
  else if (type == REDUCED_CONT_RAISE_BASIS1 || type == REDUCED_CONT_RAISE_BASIS2)
  {
    // Order-elevate basis1 such that it is of one degree higher than basis2
    // but only C^p-2 continuous
    result[0].reset(new Go::SplineVolume(*svol));
    result[0]->raiseOrder(1,1,1);
    result[1].reset(new Go::SplineVolume(*svol));
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
    geoBasis = 4;
  } else if (type == SUBGRID) {
    // basis1 should be one degree higher than basis2 and C^p-1 continuous
    int ndim = svol->dimension();
    result[1].reset(new Go::SplineVolume(*svol));
    Go::SplineVolume tmp(*svol);
    for (size_t dir = 0; dir < 3; ++dir) {
      RealArray extraKnots;
      RealArray::const_iterator uit = tmp.basis(dir).begin();
      double ucurr, uprev = *(uit++);
      while (uit != tmp.basis(dir).end())
      {
        ucurr = *(uit++);
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
    geoBasis = 1;
  }

  if (type == FULL_CONT_RAISE_BASIS2 || type == REDUCED_CONT_RAISE_BASIS2)
    std::swap(result[0], result[1]);

  return result;
}
