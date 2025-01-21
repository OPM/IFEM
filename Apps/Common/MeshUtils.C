//==============================================================================
//!
//! \file MeshUtils.C
//!
//! \date Feb 16 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various mesh quality indicators
//!
//==============================================================================

#include "MeshUtils.h"
#include "ASMbase.h"
#include "SIMbase.h"
#include "Vec3.h"
#include "Vec3Oper.h"


//! \brief A function calculating a quantity for a single element
typedef double(*CellFunction)(const ASMbase& patch, int iel);


//! \brief Function running over mesh cells calling the cellfunction for each element
static bool compute(std::vector<double>& result, const SIMbase& model,
                    const Vector& displacement, CellFunction func)
{
  result.resize(model.getNoElms(),0.0);

  for (int idx = 1; idx <= model.getNoPatches(); idx++) {
    ASMbase* pch = model.getPatch(idx,true);
    if (!pch) continue;

    Vector locvec;
    if (!displacement.empty()) {
      model.extractPatchSolution(displacement,locvec,pch);
      pch->updateCoords(locvec);
    }

    size_t nel = pch->getNoElms(true);
    for (size_t e = 1; e <= nel; e++) {
      int iel = pch->getElmID(e);
      if (iel > 0)
        result[iel-1] = func(*pch,e);
    }

    if (!displacement.empty()) {
      locvec *= -1.0;
      pch->updateCoords(locvec);
    }
  }

  return true;
}

//! \brief Compute the aspect ratio of a cell
//! \details The aspect ratio is defined as the longest edge divided
//!          by the shortest edge
static double aspectRatio(const ASMbase& patch, int iel)
{
  Matrix X;       // Control point coordinates in element
  patch.getElementCoordinates(X,iel);

  if (patch.getNoSpaceDim() == 2) {
    Vec3 ll(X(1,1), X(2,1), 0.0);
    Vec3 lr(X(1,3), X(2,3), 0.0);
    Vec3 tl(X(1,2), X(2,2), 0.0);
    Vec3 tr(X(1,4), X(2,4), 0.0);
    std::vector<double> e(4);
    e[0] = (tl-ll).length();
    e[1] = (tr-lr).length();
    e[2] = (lr-ll).length();
    e[3] = (tr-tl).length();
    return *std::max_element(e.begin(),e.end())/
      *std::min_element(e.begin(),e.end());
  }

  return 1.0;
}

//! \brief Compute the skewness of a cell
//! \details The skewness measures the deviance from a regular
//!          cell in terms of angles
static double skewness(const ASMbase& patch, int iel)
{
  Matrix X;       // Control point coordinates in element
  patch.getElementCoordinates(X,iel);

  if (patch.getNoSpaceDim() == 2) {
    Vec3 ll(X(1,1), X(2,1), 0.0);
    Vec3 lr(X(1,3), X(2,3), 0.0);
    Vec3 tl(X(1,2), X(2,2), 0.0);
    Vec3 tr(X(1,4), X(2,4), 0.0);
    Vec3 v1 = tl-ll;
    Vec3 v2 = tr-lr;
    Vec3 v3 = lr-ll;
    Vec3 v4 = tr-tl;
    std::vector<double> t(4);
    t[0] = acos(v1*v3/(v1.length()*v3.length()))/M_PI*180.0;
    t[1] = acos(v3*v2/(v3.length()*v2.length()))/M_PI*180.0;
    t[2] = acos(v2*v4/(v2.length()*v4.length()))/M_PI*180.0;
    t[3] = acos(v4*v1/(v1.length()*v4.length()))/M_PI*180.0;
    std::vector<double>::const_iterator min = std::min_element(t.begin(), t.end());
    std::vector<double>::const_iterator max = std::max_element(t.begin(), t.end());
    return std::max((*max-90.0)/90.0, (90.0-*min)/90.0);
  }

  return 0.0;
}


namespace MeshUtils
{
  bool computeAspectRatios(std::vector<double>& elmAspects,
                           const SIMbase& model, const Vector& displacement)
  {
    return compute(elmAspects, model, displacement, aspectRatio);
  }

  bool computeMeshSkewness(std::vector<double>& elmSkewness,
                           const SIMbase& model, const Vector& displacement)
  {
    return compute(elmSkewness, model, displacement, skewness);
  }
}
