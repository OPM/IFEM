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
  }
  else if (type == REDUCED_CONT_RAISE_BASIS1 || type == REDUCED_CONT_RAISE_BASIS2)
  {
    // Order-elevate basis1 such that it is of one degree higher than basis2
    // but only C^p-2 continuous
    result[0].reset(new Go::SplineSurface(*surf));
    result[0]->raiseOrder(1,1);
  }
  result[1].reset(new Go::SplineSurface(*surf));

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
      svol->gridEvaluator(XYZ,ug,vg,wg);
      // Project the coordinates onto the new basis (the 2nd XYZ is dummy here)
      result[0].reset(Go::VolumeInterpolator::regularInterpolation(b1,b2,b3,
                                                                   ug,vg,wg,XYZ,ndim,
                                                                   false,XYZ));
    }
  }
  else if (type == REDUCED_CONT_RAISE_BASIS1 || type == REDUCED_CONT_RAISE_BASIS2)
  {
    // Order-elevate basis1 such that it is of one degree higher than basis2
    // but only C^p-2 continuous
    result[0].reset(new Go::SplineVolume(*svol));
    result[0]->raiseOrder(1,1,1);
  }
  result[1].reset(new Go::SplineVolume(*svol));

  if (type == FULL_CONT_RAISE_BASIS2 || type == REDUCED_CONT_RAISE_BASIS2)
    std::swap(result[0], result[1]);

  return result;
}
