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
#include "ElementBlock.h"
#include <iostream>
#include <array>
#include <cmath>


int  Immersed::stabilization = Immersed::NO_STAB;
bool Immersed::plotCells = false;


/*!
  \brief A struct representing a point in 2D space.
*/

struct vertex
{
  double x; //!< Global x-coordinate
  double y; //!< Global y-coordinate

  //! \brief Default constructor.
  vertex(double X = 0.0, double Y = 0.0) : x(X), y(Y) {}
};

//! \brief Addition of two points.
vertex operator+ (const vertex& a, const vertex& b)
{
  return vertex(a.x+b.x,a.y+b.y);
}

//! \brief Scaling a point coordinate.
vertex operator* (double a, const vertex& b)
{
  return vertex(a*b.x,a*b.y);
}


/*!
  \brief A struct representing an integration cell.
  \details This struct represents the components of each integration cell that
  are stored during refinement of one element.
*/

struct cell
{
  int    depth;        //!< Depth of the current integration cell
  double xi;           //!< xi-coordinate of the cell midpoint
  double eta;          //!< eta-coordinate of the cell midpoint
  vertex CellVerts[4]; //!< Global coordinates of the cell vertices

  //! \brief Default constructor.
  explicit cell(int level = 0) : depth(level), xi(0.0), eta(0.0) {}
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
                                    RealArray& GPw,
                                    ElementBlock* grid)
{
  GP1.clear();
  GP2.clear();
  GPw.clear();

  // Get Gauss point coordinates and weights in 1D, in the domain [-1,1]
  const double* LocGPxi = GaussQuadrature::getCoord(nGauss);
  const double* LocGPw = GaussQuadrature::getWeight(nGauss);
  if (!LocGPxi || !LocGPw) return false;

  // Fill initial cell of depth 0
  cell cell0;
  cell0.CellVerts[0] = vertex(x1,y1);
  cell0.CellVerts[1] = vertex(x2,y2);
  cell0.CellVerts[2] = vertex(x3,y3);
  cell0.CellVerts[3] = vertex(x4,y4);

  std::vector<cell> CellSet; // Vector that will contain the cells
  CellSet.push_back(cell0);

  // Find length of element edges in physical space
  double hx1 = hypot(x2-x1,y2-y1);
  double hx2 = hypot(x3-x4,y3-y4);
  double hy1 = hypot(x4-x1,y4-y1);
  double hy2 = hypot(x3-x2,y3-y2);

  // Loop over levels
  for (int i_depth = 1; i_depth <= max_depth; i_depth++) {
    double cellScale = pow(2.0,double(i_depth+1));

    // Define the tolerance epsilon of current depth.
    // It defines the offset from each vertex during the inside-outside test.
    // The offset prevents that the cell will be refined when the boundary
    // touches only one vertex.
    double epsx = 0.002*(hx1+hx2)/cellScale;
    double epsy = 0.002*(hy1+hy2)/cellScale;

    // Loop over all cells of the current depth
    for (int i_cell = CellSet.size()-1; i_cell >= 0; i_cell--) {
      cell& curCell = CellSet[i_cell];

      // Check depth
      // If the cell is not part of the current depth -> continue
      if (curCell.depth != i_depth-1) continue;

      // Inside-outside test
      // Check, if vertices are all inside or all outside -> otherwise: Refine!
      int counter = 0;
      double alpha_start = geo.Alpha(curCell.CellVerts[0].x+epsx,
                                     curCell.CellVerts[0].y+epsy);
      if (alpha_start == geo.Alpha(curCell.CellVerts[1].x-epsx,
                                   curCell.CellVerts[1].y+epsy)) counter++;
      if (alpha_start == geo.Alpha(curCell.CellVerts[2].x-epsx,
                                   curCell.CellVerts[2].y-epsy)) counter++;
      if (alpha_start == geo.Alpha(curCell.CellVerts[3].x+epsx,
                                   curCell.CellVerts[3].y-epsy)) counter++;
      if (counter == 3) continue;

      // If all tests are passed, cell needs to be refined
      // Refine in the usual order (as suggested by Trond)

      // First cell
      // ----------
      cell cell1(i_depth);
      cell1.CellVerts[0] = curCell.CellVerts[0];
      cell1.CellVerts[1] = 0.50*(curCell.CellVerts[0]+curCell.CellVerts[1]);
      cell1.CellVerts[2] = 0.25*(curCell.CellVerts[0]+curCell.CellVerts[1]+
                                 curCell.CellVerts[2]+curCell.CellVerts[3]);
      cell1.CellVerts[3] = 0.50*(curCell.CellVerts[0]+curCell.CellVerts[3]);

      cell1.xi  = curCell.xi  - 2.0/cellScale;
      cell1.eta = curCell.eta - 2.0/cellScale;

      // Second cell
      // -----------
      cell cell2(i_depth);
      cell2.CellVerts[0] = cell1.CellVerts[1];
      cell2.CellVerts[1] = curCell.CellVerts[1];
      cell2.CellVerts[2] = 0.5*(curCell.CellVerts[1]+curCell.CellVerts[2]);
      cell2.CellVerts[3] = cell1.CellVerts[2];

      cell2.xi  = curCell.xi  + 2.0/cellScale;
      cell2.eta = curCell.eta - 2.0/cellScale;

      // Third cell
      // ----------
      cell cell3(i_depth);
      cell3.CellVerts[0] = cell2.CellVerts[3];
      cell3.CellVerts[1] = cell2.CellVerts[2];
      cell3.CellVerts[2] = curCell.CellVerts[2];
      cell3.CellVerts[3] = 0.5*(curCell.CellVerts[2]+curCell.CellVerts[3]);

      cell3.xi  = curCell.xi  + 2.0/cellScale;
      cell3.eta = curCell.eta + 2.0/cellScale;

      // Fourth cell (change data in current cell automatically deletes the old)
      // -----------------------------------------------------------------------
      curCell.depth = i_depth;
      curCell.CellVerts[0] = cell1.CellVerts[3];
      curCell.CellVerts[1] = cell3.CellVerts[0];
      curCell.CellVerts[2] = cell3.CellVerts[3];

      curCell.xi  = curCell.xi  - 2.0/cellScale;
      curCell.eta = curCell.eta + 2.0/cellScale;

      // Push-back the new cells.
      // Must do that here to not invalidate the curCell reference.
      CellSet.push_back(cell1);
      CellSet.push_back(cell2);
      CellSet.push_back(cell3);

      if (grid)
      {
        // Add grid lines for plotting
        grid->addLine(cell1.CellVerts[1].x,cell1.CellVerts[1].y,0.001,
                      cell3.CellVerts[3].x,cell3.CellVerts[3].y,0.001);
        grid->addLine(cell1.CellVerts[3].x,cell1.CellVerts[3].y,0.001,
                      cell2.CellVerts[2].x,cell2.CellVerts[2].y,0.001);
      }
    }
  }

  // Compute Gauss point coordinates
  for (size_t i_cell = 0; i_cell < CellSet.size(); i_cell++) {
    const cell& c = CellSet[i_cell];
    double cellScale = pow(2.0,double(c.depth+1));
    for (int jGP = 0; jGP < nGauss; jGP++)
      for (int iGP = 0; iGP < nGauss; iGP++) {

        // Compute global x,y-coordinates of current Gauss point
        // using the bilinear interpolation functions, N1...N4
        double xi  = LocGPxi[iGP];
        double eta = LocGPxi[jGP];
        double N1  = 0.25*(1.0-xi)*(1.0-eta);
        double N2  = 0.25*(xi+1.0)*(1.0-eta);
        double N3  = 0.25*(xi+1.0)*(eta+1.0);
        double N4  = 0.25*(1.0-xi)*(eta+1.0);
        vertex Xg  = N1*c.CellVerts[0] + N2*c.CellVerts[1]
                   + N3*c.CellVerts[2] + N4*c.CellVerts[3];

        // Do inside-outside test, if Gauss point is inside -> append to list
        if (geo.Alpha(Xg.x,Xg.y) > 0.0) {
          GP1.push_back(CellSet[i_cell].xi  + LocGPxi[iGP]*2.0/cellScale);
          GP2.push_back(CellSet[i_cell].eta + LocGPxi[jGP]*2.0/cellScale);
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
                                    Real3DMat& quadPoints,
                                    ElementBlock* grid)
{
  bool ok = true;
  quadPoints.resize(elmCorner.size());
  for (size_t e = 0; e < elmCorner.size(); e++)
  {
    int nsd = 0;
    std::array<RealArray,4> GP;
    const Real2DMat& X = elmCorner[e];
    switch (X.size()) {
    case 4: // 2D element
      nsd = 2;
      ok = getQuadraturePoints(geometry,
                               X[0][0],X[0][1],
                               X[1][0],X[1][1],
                               X[3][0],X[3][1],
                               X[2][0],X[2][1],max_depth,p,
                               GP[1],GP[2],GP[0],grid);
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
