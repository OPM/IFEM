// $Id$
//==============================================================================
//!
//! \file ImmersedBoundaries.C
//!
//! \date Dec 18 2012
//!
//! \author Dominik Schillinger / ICES and Knut Morten Okstad / SINTEF
//!
//! \brief Utilities for immersed boundary calculations.
//!
//==============================================================================

#include "ImmersedBoundaries.h"
#include "GaussQuadrature.h"
#include <iostream>
#include <cmath>


/*!
  \brief A struct representing an integration cell.
  \details This struct represents the components of each integration cell that
  are stored during refinement of one element.
  Note that index 1 (as e.g. in MidGlob1) denotes a quantity in x or xi
  direction, index 2 in y or eta direction.
*/

struct cell
{
  int    depth;         //!< Depth of the current integration cell
  double MidGlob1;      //!< x-coordinate of the cell midpoint
  double MidGlob2;      //!< y-coordinate of the cell midpoint
  double MidLoc1;       //!< xi-coordinate of the cell midpoint
  double MidLoc2;       //!< eta-coordinate of the cell midpoint
  double CellVerts1[4]; //!< x-coordinates of the cell vertices
  double CellVerts2[4]; //!< y-coordinates of the cell vertices
  //! \brief Default constructor.
  cell() : depth(0), MidGlob1(0.0), MidGlob2(0.0), MidLoc1(0.0), MidLoc2(0.0) {}
};


// ---------------------------------------------------------------
// 4. FUNCTION getQuadraturePoints
// ---------------------------------------------------------------

/*!
  This is the function that you will call:
  Following our discussion, you need to specify the following:
  Global x,y-coordinate pairs of the 4 vertices of the element.
  Maximum depth up to which you want to refine.
  Order of the Gauss integration.
  3 arrays that contain xi,eta-coordinates and weights of the Gauss points.
*/

bool Immersed::getQuadraturePoints (const Geometry& geo,
			            double x1, double y1,
                                    double x2, double y2,
                                    double x3, double y3,
                                    double x4, double y4,
                                    int max_depth, int nGauss,
                                    RealArray& GP1,
                                    RealArray& GP2,
                                    RealArray& GPw)
{
  GP1.clear();
  GP2.clear();
  GPw.clear();

  // Get Gauss point coordinates and weights in 1D, in the domain [-1,1]
  const double* LocGPxi = GaussQuadrature::getCoord(nGauss);
  const double* LocGPw = GaussQuadrature::getWeight(nGauss);
  if (!LocGPxi || !LocGPw) return false;

  std::vector<cell> CellSet(1); // Vector that will contain the cells

  // Fill initial cell of depth 0
  CellSet[0].CellVerts1[0] = x1;
  CellSet[0].CellVerts2[0] = y1;
  CellSet[0].CellVerts1[1] = x2;
  CellSet[0].CellVerts2[1] = y2;
  CellSet[0].CellVerts1[2] = x3;
  CellSet[0].CellVerts2[2] = y3;
  CellSet[0].CellVerts1[3] = x4;
  CellSet[0].CellVerts2[3] = y4;

  CellSet[0].MidGlob1 = (CellSet[0].CellVerts1[0]+CellSet[0].CellVerts1[1]) / 2.0;
  CellSet[0].MidGlob2 = (CellSet[0].CellVerts2[0]+CellSet[0].CellVerts2[3]) / 2.0;

  // Find length of element edges in physical space
  double hx = CellSet[0].CellVerts1[1]-CellSet[0].CellVerts1[0];
  double hy = CellSet[0].CellVerts2[3]-CellSet[0].CellVerts2[0];

  // Loop over levels
  for (int i_depth = 1; i_depth <= max_depth; i_depth++) {

    // Define epsilon of current depth
    // epsilon defines the off-set from each vertex during the inside-outside test
    // The off-set prevents that the cell will be refined when the boundary touches only one vertex.
    double eps1 = 0.001*hx/double(i_depth);
    double eps2 = 0.001*hy/double(i_depth);

    // Loop over all cells of the current depth
    for (int i_cell = CellSet.size()-1; i_cell >= 0; i_cell--) {

      // Check depth
      // If the cell is not part of the current depth -> continue
      if (CellSet[i_cell].depth != i_depth-1) continue;

      // Inside-outside test
      // Check, if vertices are all inside or all outside -> otherwise: Refine!
      int counter = 0;
      double alpha_start = geo.Alpha(CellSet[i_cell].CellVerts1[0]+eps1,CellSet[i_cell].CellVerts2[0]+eps2);
      if (alpha_start == geo.Alpha(CellSet[i_cell].CellVerts1[1]-eps1,CellSet[i_cell].CellVerts2[1]+eps2)) counter++;
      if (alpha_start == geo.Alpha(CellSet[i_cell].CellVerts1[2]-eps1,CellSet[i_cell].CellVerts2[2]-eps2)) counter++;
      if (alpha_start == geo.Alpha(CellSet[i_cell].CellVerts1[3]+eps1,CellSet[i_cell].CellVerts2[3]-eps2)) counter++;
      if (counter == 3) continue;

      // If all tests are passed, cell needs to be refined
      // Refine in the usual order (as suggested by Trond)
      double cellScale = pow(2.0,double(i_depth+1));

      // First cell (push_back to vector of cells)
      // ---------------------------------------------------------------------------
      CellSet.push_back(CellSet[i_cell]);
      CellSet.back().depth = i_depth;
      CellSet.back().CellVerts1[1] = CellSet[i_cell].MidGlob1;
      CellSet.back().CellVerts1[2] = CellSet[i_cell].MidGlob1;
      CellSet.back().CellVerts2[2] = CellSet[i_cell].MidGlob2;
      CellSet.back().CellVerts2[3] = CellSet[i_cell].MidGlob2;

      CellSet.back().MidGlob1 = CellSet[i_cell].MidGlob1 - hx/cellScale;
      CellSet.back().MidGlob2 = CellSet[i_cell].MidGlob2 - hy/cellScale;

      CellSet.back().MidLoc1  = CellSet[i_cell].MidLoc1 - 2.0/cellScale;
      CellSet.back().MidLoc2  = CellSet[i_cell].MidLoc2 - 2.0/cellScale;

      // Second cell (push_back to vector of cells)
      // -------------------------------------------------------------------
      CellSet.push_back(CellSet[i_cell]);
      CellSet.back().depth = i_depth;
      CellSet.back().CellVerts1[0] = CellSet[i_cell].MidGlob1;
      CellSet.back().CellVerts2[2] = CellSet[i_cell].MidGlob2;
      CellSet.back().CellVerts1[3] = CellSet[i_cell].MidGlob1;
      CellSet.back().CellVerts2[3] = CellSet[i_cell].MidGlob2;

      CellSet.back().MidGlob1 = CellSet[i_cell].MidGlob1 + hx/cellScale;
      CellSet.back().MidGlob2 = CellSet[i_cell].MidGlob2 - hy/cellScale;

      CellSet.back().MidLoc1  = CellSet[i_cell].MidLoc1 + 2.0/cellScale;
      CellSet.back().MidLoc2  = CellSet[i_cell].MidLoc2 - 2.0/cellScale;

      // Third cell (push_back to vector of cells)
      // --------------------------------------------------------------------------------
      CellSet.push_back(CellSet[i_cell]);
      CellSet.back().depth = i_depth;
      CellSet.back().CellVerts1[0] = CellSet[i_cell].MidGlob1;
      CellSet.back().CellVerts2[0] = CellSet[i_cell].MidGlob2;
      CellSet.back().CellVerts2[1] = CellSet[i_cell].MidGlob2;
      CellSet.back().CellVerts1[3] = CellSet[i_cell].MidGlob1;

      CellSet.back().MidGlob1 = CellSet[i_cell].MidGlob1 + hx/cellScale;
      CellSet.back().MidGlob2 = CellSet[i_cell].MidGlob2 + hy/cellScale;

      CellSet.back().MidLoc1  = CellSet[i_cell].MidLoc1 + 2.0/cellScale;
      CellSet.back().MidLoc2  = CellSet[i_cell].MidLoc2 + 2.0/cellScale;

      // Fourth cell (change data in current cell -> automatically deletes old cell)
      // --------------------------------------------------------------------------------------
      CellSet[i_cell].depth = i_depth;
      CellSet[i_cell].CellVerts2[0] = CellSet[i_cell].MidGlob2;
      CellSet[i_cell].CellVerts1[1] = CellSet[i_cell].MidGlob1;
      CellSet[i_cell].CellVerts2[1] = CellSet[i_cell].MidGlob2;
      CellSet[i_cell].CellVerts1[2] = CellSet[i_cell].MidGlob1;

      CellSet[i_cell].MidGlob1 = CellSet[i_cell].MidGlob1 - hx/cellScale;
      CellSet[i_cell].MidGlob2 = CellSet[i_cell].MidGlob2 + hy/cellScale;

      CellSet[i_cell].MidLoc1  = CellSet[i_cell].MidLoc1 - 2.0/cellScale;
      CellSet[i_cell].MidLoc2  = CellSet[i_cell].MidLoc2 + 2.0/cellScale;
    }
  }

  // Compute Gauss points
  // Loop over all cells
  for (size_t i_cell = 0; i_cell < CellSet.size(); i_cell++) {
    double cellScale = pow(2.0,double(CellSet[i_cell].depth+1));
    for (int jGP = 0; jGP < nGauss; jGP++)
      for (int iGP = 0; iGP < nGauss; iGP++) {

	// Compute global x,y-coordinates of current Gauss point
	double GlobGP1 = CellSet[i_cell].MidGlob1 + LocGPxi[iGP]*hx/cellScale;
	double GlobGP2 = CellSet[i_cell].MidGlob2 + LocGPxi[jGP]*hy/cellScale;

	// Do inside-outside test, if Gauss point is inside -> append to list
	if (geo.Alpha(GlobGP1,GlobGP2) > 0.0) {
	  GP1.push_back(CellSet[i_cell].MidLoc1 + LocGPxi[iGP]*2.0/cellScale);
	  GP2.push_back(CellSet[i_cell].MidLoc2 + LocGPxi[jGP]*2.0/cellScale);
	  GPw.push_back(LocGPw[iGP]*LocGPw[jGP]*4.0/(cellScale*cellScale));
	}
      }
  }

  return true;
}


// ---------------------------------------------------------------
// 5. FUNCTION getQuadraturePoints (3D version)
// ---------------------------------------------------------------

/*!
  This is the function that you will call:
  Following our discussion, you need to specify the following:
  Global x,y,z-coordinate pairs of the 8 vertices of the element.
  Maximum depth up to which you want to refine.
  Order of the Gauss integration.
  4 arrays that contain xi,eta,zeta-coordinates and weights of the Gauss points.
*/

bool Immersed::getQuadraturePoints (const Geometry& geo,
			            double x1, double y1, double z1,
				    double x2, double y2, double z2,
				    double x3, double y3, double z3,
				    double x4, double y4, double z4,
				    double x5, double y5, double z5,
				    double x6, double y6, double z6,
				    double x7, double y7, double z7,
				    double x8, double y8, double z8,
				    int max_depth, int nGauss,
				    RealArray& GP1,
				    RealArray& GP2,
				    RealArray& GP3,
				    RealArray& GPw)
{
  return false; // TODO 3D version...
}


// Wrapper for processing multiple elements.

bool Immersed::getQuadraturePoints (const Geometry& geometry,
			            const Real3DMat& elmCorner,
				    int max_depth, int p,
				    Real3DMat& quadPoints)
{
  bool ok = true;
  quadPoints.resize(elmCorner.size());
  for (size_t e = 0; e < elmCorner.size(); e++)
  {
    int nsd = 0;
    RealArray GP[4];
    const Real2DMat& X = elmCorner[e];
    switch (X.size()) {
    case 4: // 2D element
      nsd = 2;
      ok = getQuadraturePoints(geometry,
			       X[0][0],X[0][1],
			       X[1][0],X[1][1],
			       X[3][0],X[3][1],
			       X[2][0],X[2][1],max_depth,p,
			       GP[1],GP[2],GP[0]);
      break;
    case 8: // 3D element
      nsd = 3;
      ok = getQuadraturePoints(geometry,
			       X[0][0],X[0][1],X[0][2],
			       X[1][0],X[1][1],X[1][2],
			       X[3][0],X[3][1],X[3][2],
			       X[2][0],X[2][1],X[2][2],
			       X[4][0],X[4][1],X[4][2],
			       X[5][0],X[5][1],X[5][2],
			       X[7][0],X[7][1],X[7][2],
			       X[6][0],X[6][1],X[6][2],max_depth,p,
			       GP[1],GP[2],GP[3],GP[0]);
      break;
    default:
      ok = false;
      GP[0].clear();
      std::cerr <<" *** Immersed::getQuadraturePoints: Invalid element ("
		<< X.size() <<" corners)."<< std::endl;
    }

    // Store the quadrature points
    RealArray xg(nsd+1);
    quadPoints[e].resize(GP[0].size());
    for (size_t i = 0; i < GP[0].size(); i++)
    {
      for (int d = 0; d < nsd; d++)
	xg[d] = GP[d+1][i];
      xg[nsd] = GP[0][i];
      quadPoints[e][i] = xg;
    }
  }

  return ok;
}
