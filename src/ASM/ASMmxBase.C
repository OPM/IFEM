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
bool ASMmxBase::includeExtra         = false;


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
    result[0].reset(ASMmxBase::adjustBasis(*surf,{SplineUtils::AdjustOp::Raise,
                                                  SplineUtils::AdjustOp::Raise}));
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
    result[0].reset(ASMmxBase::adjustBasis(*surf,{SplineUtils::AdjustOp::Original,
                                                  SplineUtils::AdjustOp::Lower}));
    result[1].reset(ASMmxBase::adjustBasis(*surf,{SplineUtils::AdjustOp::Lower,
                                                  SplineUtils::AdjustOp::Original}));
    result[2].reset(ASMmxBase::adjustBasis(*surf,{SplineUtils::AdjustOp::Lower,
                                                  SplineUtils::AdjustOp::Lower}));
    itgBasis = 3;
  } else if (type == SUBGRID) {
    // basis1 should be one degree higher than basis2 and C^p-1 continuous
    result[1].reset(new Go::SplineSurface(*surf));
    result[0].reset(ASMmxBase::adjustBasis(*surf,{SplineUtils::AdjustOp::Raise,
                                                  SplineUtils::AdjustOp::Raise}));
    for (size_t i = 0; i < 2; ++i) {
      RealArray extraKnots;
      RealArray::const_iterator uit = result[0]->basis(i).begin();
      double uprev = *(uit++);
      while (uit != result[0]->basis(i).end())
      {
        double ucurr = *(uit++);
        if (ucurr > uprev)
          extraKnots.push_back(ucurr*0.5 + uprev*0.5);
        uprev = ucurr;
      }
      if (i == 0)
        result[0]->insertKnot_u(extraKnots);
      else
        result[0]->insertKnot_v(extraKnots);
    }
    itgBasis = 1;
  }

  if (type == FULL_CONT_RAISE_BASIS2 || type == REDUCED_CONT_RAISE_BASIS2)
    std::swap(result[0], result[1]), itgBasis = 1;

  return result;
}


ASMmxBase::VolumeVec ASMmxBase::establishBases (Go::SplineVolume* svol,
                                                MixedType type)
{
  VolumeVec result(2);
  // With mixed methods we need two separate spline spaces
  if (type == FULL_CONT_RAISE_BASIS1 || type == FULL_CONT_RAISE_BASIS2)
  {
    // basis1 should be one degree higher than basis2 and C^p-1 continuous
    result[0].reset(ASMmxBase::adjustBasis(*svol,{SplineUtils::AdjustOp::Raise,
                                                  SplineUtils::AdjustOp::Raise,
                                                  SplineUtils::AdjustOp::Raise}));
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
    result[0].reset(ASMmxBase::adjustBasis(*svol,{SplineUtils::AdjustOp::Original,
                                                  SplineUtils::AdjustOp::Lower,
                                                  SplineUtils::AdjustOp::Lower}));
    result[1].reset(ASMmxBase::adjustBasis(*svol,{SplineUtils::AdjustOp::Lower,
                                                  SplineUtils::AdjustOp::Original,
                                                  SplineUtils::AdjustOp::Lower}));
    result[2].reset(ASMmxBase::adjustBasis(*svol,{SplineUtils::AdjustOp::Lower,
                                                  SplineUtils::AdjustOp::Lower,
                                                  SplineUtils::AdjustOp::Original}));
    result[3].reset(ASMmxBase::adjustBasis(*svol,{SplineUtils::AdjustOp::Lower,
                                                  SplineUtils::AdjustOp::Lower,
                                                  SplineUtils::AdjustOp::Lower}));
    itgBasis = 4;
  } else if (type == SUBGRID) {
    // basis1 should be one degree higher than basis2 and C^p-1 continuous
    result[1].reset(new Go::SplineVolume(*svol));
    result[0].reset(ASMmxBase::adjustBasis(*svol,{SplineUtils::AdjustOp::Raise,
                                                  SplineUtils::AdjustOp::Raise,
                                                  SplineUtils::AdjustOp::Raise}));
    for (size_t dir = 0; dir < 3; ++dir) {
      RealArray extraKnots;
      RealArray::const_iterator uit = result[0]->basis(dir).begin();
      double uprev = *(uit++);
      while (uit != result[0]->basis(dir).end())
      {
        double ucurr = *(uit++);
        if (ucurr > uprev)
          extraKnots.push_back(ucurr*0.5 + uprev*0.5);
        uprev = ucurr;
      }
      result[0]->insertKnot(dir, extraKnots);
    }
    itgBasis = 1;
  }

  if (type == FULL_CONT_RAISE_BASIS2 || type == REDUCED_CONT_RAISE_BASIS2)
    std::swap(result[0], result[1]), itgBasis = 1;

  return result;
}


Go::SplineSurface* ASMmxBase::adjustBasis (const Go::SplineSurface& surf,
                                           const std::array<SplineUtils::AdjustOp,2>& ops)
{
  // Create a surface with adjusted basis from surf
  // while keeping lines of reduced continuity
  std::array<Go::BsplineBasis,2> basis{
    SplineUtils::adjustBasis(surf.basis(0),ops[0]),
    SplineUtils::adjustBasis(surf.basis(1),ops[1])
  };

  // Compute parameter values of the Greville points
  std::array<RealArray,2> ug;
  for (size_t idx = 0; idx < 2; ++idx) {
    ug[idx].resize(basis[idx].numCoefs());
    for (size_t i = 0; i < ug[idx].size(); i++)
      ug[idx][i] = basis[idx].grevilleParameter(i);
  }

  int ndim = surf.dimension();
  if (surf.rational())
  {
    std::vector<double> rCoefs(surf.rcoefs_begin(), surf.rcoefs_end());

    // We normally would set coefs as (x*w, y*w, w).
    // However, GoTools uses this representation internally already.

    // Instance a Bspline surface in ndim+1
    Go::SplineSurface surf2(surf.basis(0), surf.basis(1),
                            rCoefs.begin(), ndim+1, false);

    // Interpolate the Bspline surface onto new basis
    RealArray XYZ((ndim+1)*ug[0].size()*ug[1].size());
    surf2.gridEvaluator(XYZ,ug[0],ug[1]);
    std::unique_ptr<Go::SplineSurface> surf3(
      Go::SurfaceInterpolator::regularInterpolation(basis[0],basis[1],
                                                    ug[0],ug[1],
                                                    XYZ,ndim+1,false,XYZ));

    // New rational coefs are (x/w', y/w', w')
    // Apparently, GoTools will rescale coeffs on surface creation.
    return new Go::SplineSurface(surf3->basis(0), surf3->basis(1),
                                 surf3->coefs_begin(), ndim, true);
  }

  // Evaluate the spline surface at all points
  RealArray XYZ(ndim*ug[0].size()*ug[1].size());
  surf.gridEvaluator(XYZ,ug[0],ug[1]);

  // Project the coordinates onto the new basis (the 2nd XYZ is dummy here)
  return Go::SurfaceInterpolator::regularInterpolation(basis[0],basis[1],
                                                      ug[0],ug[1],XYZ,ndim,false,XYZ);
}


Go::SplineVolume* ASMmxBase::adjustBasis (const Go::SplineVolume& svol,
                                          const std::array<SplineUtils::AdjustOp,3>& ops)
{
  // Create a volume with adjusted orders from svol
  // while keeping lines of reduced continuity
  std::array<Go::BsplineBasis,3> basis{
    SplineUtils::adjustBasis(svol.basis(0),ops[0]),
    SplineUtils::adjustBasis(svol.basis(1),ops[1]),
    SplineUtils::adjustBasis(svol.basis(2),ops[2])
  };

  // Compute parameter values of the Greville points
  std::array<RealArray,3> ug;
  for (size_t idx = 0; idx < 3; ++idx) {
    ug[idx].resize(basis[idx].numCoefs());
    for (size_t i = 0; i < ug[idx].size(); i++)
      ug[idx][i] = basis[idx].grevilleParameter(i);
  }

  int ndim = svol.dimension();
  if (svol.rational())
  {
    std::vector<double> rCoefs(svol.rcoefs_begin(), svol.rcoefs_end());

    // We normally would set coefs as (x*w, y*w, w).
    // However, GoTools use this representation internally already.

    // Instance a Bspline surface in ndim+1
    Go::SplineVolume vol2(svol.basis(0), svol.basis(1), svol.basis(2),
                          rCoefs.begin(), ndim+1, false);

    // Interpolate the Bspline surface onto new basis
    RealArray XYZ((ndim+1)*ug[0].size()*ug[1].size()*ug[2].size());
    vol2.gridEvaluator(XYZ,ug[0],ug[1],ug[2]);
    std::unique_ptr<Go::SplineVolume> svol3(
      Go::VolumeInterpolator::regularInterpolation(basis[0],basis[1],basis[2],
                                                   ug[0],ug[1],ug[2],XYZ,
                                                   ndim+1,false,XYZ));

    // New rational coefs are (x/w', y/w', w')
    // Apparently, GoTools will rescale coeffs on surface creation.
    return new Go::SplineVolume(svol3->basis(0), svol3->basis(1), svol3->basis(2),
                                svol3->coefs_begin(), ndim, true);
  }

  RealArray XYZ(ndim*ug[0].size()*ug[1].size()*ug[2].size());
  // Evaluate the spline surface at all points
  svol.gridEvaluator(ug[0],ug[1],ug[2],XYZ);
  // Project the coordinates onto the new basis (the 2nd XYZ is dummy here)
  return Go::VolumeInterpolator::regularInterpolation(basis[0],basis[1],basis[2],
                                                      ug[0],ug[1],ug[2],
                                                      XYZ,ndim,false,XYZ);
}
