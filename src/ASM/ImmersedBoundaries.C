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
#include <iostream>
#include <cmath>


// ---------------------------------------
// 1. CLASS cell:
// ---------------------------------------

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
};


// ---------------------------------------
// 2. FUNCTION FillLocalGP:
// ---------------------------------------

/*!
  \brief 1D Gauss coordinates and weights in parametric space [-1,1].
  \details This function provides the Gauss point information from order 1 to 9.
  You can replace it with your own function, if you want.
  I provided a bit more orders so that you can play with different strategies
  to integrate cut elements (either using higher p or increasing the depth).
*/

static bool FillLocalGP (int GO, RealArray& GPxi, RealArray& GPw)
{
  GPxi.reserve(GO); GPw.reserve(GO);

  switch (GO) {
  case 1:
    GPxi.push_back( 0.00000000000000e-01);  GPw.push_back(2.000000000000000e+00);
    break;
  case 2:
    GPxi.push_back(-5.773502691896258e-01); GPw.push_back(1.000000000000000e+00);
    GPxi.push_back( 5.773502691896258e-01); GPw.push_back(1.000000000000000e+00);
    break;
  case 3:
    GPxi.push_back(-7.745966692414834e-01); GPw.push_back(5.555555555555556e-01);
    GPxi.push_back( 0.000000000000000e-01); GPw.push_back(8.888888888888889e-01);
    GPxi.push_back( 7.745966692414834e-01); GPw.push_back(5.555555555555556e-01);
    break;
  case 4:
    GPxi.push_back(-8.611363115940526e-01); GPw.push_back(3.478548451374539e-01);
    GPxi.push_back(-3.399810435848563e-01); GPw.push_back(6.521451548625461e-01);
    GPxi.push_back( 3.399810435848563e-01); GPw.push_back(6.521451548625461e-01);
    GPxi.push_back( 8.611363115940526e-01); GPw.push_back(3.478548451374539e-01);
    break;
  case 5:
    GPxi.push_back(-9.061798459386640e-01); GPw.push_back(2.369268850561891e-01);
    GPxi.push_back(-5.384693101056831e-01); GPw.push_back(4.786286704993665e-01);
    GPxi.push_back( 0.000000000000000e-01); GPw.push_back(5.688888888888889e-01);
    GPxi.push_back( 5.384693101056831e-01); GPw.push_back(4.786286704993665e-01);
    GPxi.push_back( 9.061798459386640e-01); GPw.push_back(2.369268850561891e-01);
    break;
  case 6:
    GPxi.push_back(-9.324695142031520e-01); GPw.push_back(1.713244923791703e-01);
    GPxi.push_back(-6.612093864662645e-01); GPw.push_back(3.607615730481386e-01);
    GPxi.push_back(-2.386191860831969e-01); GPw.push_back(4.679139345726910e-01);
    GPxi.push_back( 2.386191860831969e-01); GPw.push_back(4.679139345726910e-01);
    GPxi.push_back( 6.612093864662645e-01); GPw.push_back(3.607615730481386e-01);
    GPxi.push_back( 9.324695142031520e-01); GPw.push_back(1.713244923791703e-01);
    break;
  case 7:
    GPxi.push_back(-9.491079123427585e-01); GPw.push_back(1.294849661688697e-01);
    GPxi.push_back(-7.415311855993944e-01); GPw.push_back(2.797053914892767e-01);
    GPxi.push_back(-4.058451513773972e-01); GPw.push_back(3.818300505051189e-01);
    GPxi.push_back( 0.000000000000000e-01); GPw.push_back(4.179591836734694e-01);
    GPxi.push_back( 4.058451513773972e-01); GPw.push_back(3.818300505051189e-01);
    GPxi.push_back( 7.415311855993944e-01); GPw.push_back(2.797053914892767e-01);
    GPxi.push_back( 9.491079123427585e-01); GPw.push_back(1.294849661688697e-01);
    break;
  case 8:
    GPxi.push_back(-9.602898564975362e-01); GPw.push_back(1.012285362903763e-01);
    GPxi.push_back(-7.966664774136267e-01); GPw.push_back(2.223810344533745e-01);
    GPxi.push_back(-5.255324099163290e-01); GPw.push_back(3.137066458778873e-01);
    GPxi.push_back(-1.834346424956498e-01); GPw.push_back(3.626837833783620e-01);
    GPxi.push_back( 1.834346424956498e-01); GPw.push_back(3.626837833783620e-01);
    GPxi.push_back( 5.255324099163290e-01); GPw.push_back(3.137066458778873e-01);
    GPxi.push_back( 7.966664774136267e-01); GPw.push_back(2.223810344533745e-01);
    GPxi.push_back( 9.602898564975362e-01); GPw.push_back(1.012285362903763e-01);
    break;
  case 9:
    GPxi.push_back(-9.681602395076261e-01); GPw.push_back(8.127438836157441e-02);
    GPxi.push_back(-8.360311073266358e-01); GPw.push_back(1.806481606948574e-01);
    GPxi.push_back(-6.133714327005904e-01); GPw.push_back(2.606106964029355e-01);
    GPxi.push_back(-3.242534234038089e-01); GPw.push_back(3.123470770400028e-01);
    GPxi.push_back( 0.000000000000000e-01); GPw.push_back(3.302393550012598e-01);
    GPxi.push_back( 3.242534234038089e-01); GPw.push_back(3.123470770400028e-01);
    GPxi.push_back( 6.133714327005904e-01); GPw.push_back(2.606106964029355e-01);
    GPxi.push_back( 8.360311073266358e-01); GPw.push_back(1.806481606948574e-01);
    GPxi.push_back( 9.681602395076261e-01); GPw.push_back(8.127438836157441e-02);
    break;
  default:
    std::cerr <<" *** Immersed::FillLocalGP: Gauss order "<< GO
              <<" exceeds maximum of 9!"<< std::endl;
    return false;
  }

  return true;
}

// ---------------------------------------
// 2.a FUNCTION FillLocalGLP:
// ---------------------------------------

/*!
  \brief 1D Gauss-Lobatto coordinates and weights in parametric space [-1,1].
  \details This function provides the Gauss-Lobatto point information
  from order 2 to 9.
*/

static bool FillLocalGLP (int GO, RealArray& GPxi, RealArray& GPw)
{
  GPxi.reserve(GO); GPw.reserve(GO);

  switch (GO) {
  case 2:
    GPxi.push_back(-1.000000000000000000); GPw.push_back(1.000000000000000e+00);
    GPxi.push_back( 1.000000000000000000); GPw.push_back(1.000000000000000e+00);
    break;
  case 3:
    GPxi.push_back(-1.000000000000000000); GPw.push_back(0.333333333333333310);
    GPxi.push_back( 0.000000000000000000); GPw.push_back(1.333333333333333300);
    GPxi.push_back( 1.000000000000000000); GPw.push_back(0.333333333333333310);
    break;
  case 4:
    GPxi.push_back(-1.000000000000000000); GPw.push_back(0.166666666666666660);
    GPxi.push_back(-0.447213595499957930); GPw.push_back(0.833333333333333370);
    GPxi.push_back( 0.447213595499957930); GPw.push_back(0.833333333333333370);
    GPxi.push_back( 1.000000000000000000); GPw.push_back(0.166666666666666660);
    break;
  case 5:
    GPxi.push_back(-1.000000000000000000); GPw.push_back(0.100000000000000010);
    GPxi.push_back(-0.654653670707977090); GPw.push_back(0.544444444444444400);
    GPxi.push_back( 0.000000000000000000); GPw.push_back(0.711111111111111140);
    GPxi.push_back( 0.654653670707977090); GPw.push_back(0.544444444444444400);
    GPxi.push_back( 1.000000000000000000); GPw.push_back(0.100000000000000010);
    break;
  case 6:
    GPxi.push_back(-1.000000000000000000); GPw.push_back(0.066666666666666666);
    GPxi.push_back(-0.765055323929464740); GPw.push_back(0.378474956297846940);
    GPxi.push_back(-0.285231516480645100); GPw.push_back(0.554858377035486460);
    GPxi.push_back( 0.285231516480645100); GPw.push_back(0.554858377035486460);
    GPxi.push_back( 0.765055323929464740); GPw.push_back(0.378474956297846940);
    GPxi.push_back( 1.000000000000000000); GPw.push_back(0.066666666666666666);
    break;
  case 7:
    GPxi.push_back(-1.000000000000000000); GPw.push_back(0.047619047619047616);
    GPxi.push_back(-0.830223896278566960); GPw.push_back(0.276826047361565910);
    GPxi.push_back(-0.468848793470714230); GPw.push_back(0.431745381209862610);
    GPxi.push_back( 0.000000000000000000); GPw.push_back(0.487619047619047620);
    GPxi.push_back( 0.468848793470714230); GPw.push_back(0.431745381209862610);
    GPxi.push_back( 0.830223896278566960); GPw.push_back(0.276826047361565910);
    GPxi.push_back( 1.000000000000000000); GPw.push_back(0.047619047619047616);
    break;
  case 8:
    GPxi.push_back(-1.000000000000000000); GPw.push_back(0.035714285714285712);
    GPxi.push_back(-0.871740148509606570); GPw.push_back(0.210704227143506150);
    GPxi.push_back(-0.591700181433142290); GPw.push_back(0.341122692483504410);
    GPxi.push_back(-0.209299217902478850); GPw.push_back(0.412458794658703720);
    GPxi.push_back( 0.209299217902478850); GPw.push_back(0.412458794658703720);
    GPxi.push_back( 0.591700181433142290); GPw.push_back(0.341122692483504410);
    GPxi.push_back( 0.871740148509606570); GPw.push_back(0.210704227143506150);
    GPxi.push_back( 1.000000000000000000); GPw.push_back(0.035714285714285712);
    break;
  case 9:
    GPxi.push_back(-1.000000000000000000); GPw.push_back(0.027777777777777776);
    GPxi.push_back(-0.899757995411460180); GPw.push_back(0.165495361560805580);
    GPxi.push_back(-0.677186279510737730); GPw.push_back(0.274538712500161650);
    GPxi.push_back(-0.363117463826178160); GPw.push_back(0.346428510973046170);
    GPxi.push_back( 0.000000000000000000); GPw.push_back(0.371519274376417240);
    GPxi.push_back( 0.363117463826178160); GPw.push_back(0.346428510973046170);
    GPxi.push_back( 0.677186279510737730); GPw.push_back(0.274538712500161650);
    GPxi.push_back( 0.899757995411460180); GPw.push_back(0.165495361560805580);
    GPxi.push_back( 1.000000000000000000); GPw.push_back(0.027777777777777776);
    break;
  default:
    std::cerr <<" *** Immersed::FillLocalGLP: Gauss-Lobatto order "<< GO
              <<" outside legal range [2,9]!"<< std::endl;
    return false;
  }

  return true;
}


// ----------------------------------------------------------------
// 3. Inside-outside test for the perforated plate benchmark
// ----------------------------------------------------------------

/*!
  (See the example Fig. 38 in the recent immersed boundary paper coming out of
  Tom's group).

  The following global information is needed:
  (a) The center of the hole; (b) the radius of the hole denoted by variable R.
  The function FindAlpha receives the global coordinates X and Y of the point
  under consideration.
  (may be a vertex of an integration element during the set up of the adaptive
  integration structure,
  or an intergration point during the integration of the stiffness matrix).
*/

static double FindAlpha (double X, double Y)
{
  // Set center and radius of the hole here
  double CenterX = 0.0;
  double CenterY = 0.0;
  double R = 1.0;

  // Alpha is used as an indicator here
  // alpha = 0.0 if the point is lying outside the physical domain
  // alpha = 0.0 if the point is lying directly on the boundary (the case Trond mentioned)
  // alpha = 1.0 if the point is lying inside the physical domain

  // Determine radius of the point to be checked
  double r = sqrt( (X-CenterX) * (X-CenterX) + (Y-CenterY) * (Y-CenterY) );

  // Determine if point is located within the hole or not
  double alpha = r > R ? 1.0 : 0.0;

  // Alpha is the penalization parameter in the sense of equation (29) in the paper
  // It can be used either to indicate where a point is located
  // Or can be used as a penalization value for the local stiffness contribution at the current integration point.
  // In the latter case, it needs to have data type double.
  return alpha;
}


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

bool Immersed::getQuadraturePoints (double x1, double y1,
                                    double x2, double y2,
                                    double x3, double y3,
                                    double x4, double y4, int max_depth, int p,
                                    RealArray& GP1,
                                    RealArray& GP2,
                                    RealArray& GPw)
{
  // Initialize some data
  std::vector<cell> CellSet(1); // Vector that will contain the cells
  RealArray LocGPxi, LocGPw; // Gauss coordinates and weights in 1D in [-1,1]
  //if (!FillLocalGLP(p+1,LocGPxi,LocGPw)) // Fill Gauss-Lobatto point information
  if (!FillLocalGP(p,LocGPxi,LocGPw)) // Fill Gauss point information
    return false;

  // Fill initial cell of depth 0
  CellSet[0].CellVerts1[0] = x1;
  CellSet[0].CellVerts2[0] = y1;
  CellSet[0].CellVerts1[1] = x2;
  CellSet[0].CellVerts2[1] = y2;
  CellSet[0].CellVerts1[2] = x3;
  CellSet[0].CellVerts2[2] = y3;
  CellSet[0].CellVerts1[3] = x4;
  CellSet[0].CellVerts2[3] = y4;

  CellSet[0].depth = 0;
  CellSet[0].MidGlob1 = (CellSet[0].CellVerts1[0]+CellSet[0].CellVerts1[1]) / 2.0;
  CellSet[0].MidGlob2 = (CellSet[0].CellVerts2[0]+CellSet[0].CellVerts2[3]) / 2.0;
  CellSet[0].MidLoc1  = 0.0;
  CellSet[0].MidLoc2  = 0.0;

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
      double alpha_start = FindAlpha(CellSet[i_cell].CellVerts1[0]+eps1,CellSet[i_cell].CellVerts2[0]+eps2);
      if (alpha_start == FindAlpha(CellSet[i_cell].CellVerts1[1]-eps1,CellSet[i_cell].CellVerts2[1]+eps2)) counter++;
      if (alpha_start == FindAlpha(CellSet[i_cell].CellVerts1[2]-eps1,CellSet[i_cell].CellVerts2[2]-eps2)) counter++;
      if (alpha_start == FindAlpha(CellSet[i_cell].CellVerts1[3]+eps1,CellSet[i_cell].CellVerts2[3]-eps2)) counter++;
      if (counter == 3) continue;

      // If all tests are passed, cell needs to be refined
      // Refine in the usual order (as suggested by Trond)
      // First cell (push_back to vector of cells)
      // ---------------------------------------------------------------------------
      CellSet.push_back(CellSet[i_cell]);
      CellSet.back().depth = i_depth;
      CellSet.back().CellVerts1[1] = CellSet[i_cell].MidGlob1;
      CellSet.back().CellVerts1[2] = CellSet[i_cell].MidGlob1;
      CellSet.back().CellVerts2[2] = CellSet[i_cell].MidGlob2;
      CellSet.back().CellVerts2[3] = CellSet[i_cell].MidGlob2;
      CellSet.back().MidGlob1 = CellSet[i_cell].MidGlob1 - hx/(pow(2.0,double(i_depth+1)));
      CellSet.back().MidGlob2 = CellSet[i_cell].MidGlob2 - hy/(pow(2.0,double(i_depth+1)));

      CellSet.back().MidLoc1  = CellSet[i_cell].MidLoc1 - 2.0/(pow(2.0,double(i_depth+1)));
      CellSet.back().MidLoc2  = CellSet[i_cell].MidLoc2 - 2.0/(pow(2.0,double(i_depth+1)));

      // Second cell (push_back to vector of cells)
      // -------------------------------------------------------------------
      CellSet.push_back(CellSet[i_cell]);
      CellSet.back().depth = i_depth;
      CellSet.back().CellVerts1[0] = CellSet[i_cell].MidGlob1;
      CellSet.back().CellVerts2[2] = CellSet[i_cell].MidGlob2;
      CellSet.back().CellVerts1[3] = CellSet[i_cell].MidGlob1;
      CellSet.back().CellVerts2[3] = CellSet[i_cell].MidGlob2;

      CellSet.back().MidGlob1 = CellSet[i_cell].MidGlob1 + hx/(pow(2.0,double(i_depth+1)));
      CellSet.back().MidGlob2 = CellSet[i_cell].MidGlob2 - hy/(pow(2.0,double(i_depth+1)));

      CellSet.back().MidLoc1  = CellSet[i_cell].MidLoc1 + 2.0/(pow(2.0,double(i_depth+1)));
      CellSet.back().MidLoc2  = CellSet[i_cell].MidLoc2 - 2.0/(pow(2.0,double(i_depth+1)));

      // Third cell (push_back to vector of cells)
      // --------------------------------------------------------------------------------
      CellSet.push_back(CellSet[i_cell]);
      CellSet.back().depth = i_depth;
      CellSet.back().CellVerts1[0] = CellSet[i_cell].MidGlob1;
      CellSet.back().CellVerts2[0] = CellSet[i_cell].MidGlob2;
      CellSet.back().CellVerts2[1] = CellSet[i_cell].MidGlob2;
      CellSet.back().CellVerts1[3] = CellSet[i_cell].MidGlob1;

      CellSet.back().MidGlob1 = CellSet[i_cell].MidGlob1 + hx/(pow(2.0,double(i_depth+1)));
      CellSet.back().MidGlob2 = CellSet[i_cell].MidGlob2 + hy/(pow(2.0,double(i_depth+1)));

      CellSet.back().MidLoc1  = CellSet[i_cell].MidLoc1 + 2.0/(pow(2.0,double(i_depth+1)));
      CellSet.back().MidLoc2  = CellSet[i_cell].MidLoc2 + 2.0/(pow(2.0,double(i_depth+1)));

      // Fourth cell (change data in current cell -> automatically deletes old cell)
      // --------------------------------------------------------------------------------------
      CellSet[i_cell].depth = i_depth;
      CellSet[i_cell].CellVerts2[0] = CellSet[i_cell].MidGlob2;
      CellSet[i_cell].CellVerts1[1] = CellSet[i_cell].MidGlob1;
      CellSet[i_cell].CellVerts2[1] = CellSet[i_cell].MidGlob2;
      CellSet[i_cell].CellVerts1[2] = CellSet[i_cell].MidGlob1;

      CellSet[i_cell].MidGlob1 = CellSet[i_cell].MidGlob1 - hx/(pow(2.0,double(i_depth+1)));
      CellSet[i_cell].MidGlob2 = CellSet[i_cell].MidGlob2 + hy/(pow(2.0,double(i_depth+1)));

      CellSet[i_cell].MidLoc1  = CellSet[i_cell].MidLoc1 - 2.0/(pow(2.0,double(i_depth+1)));
      CellSet[i_cell].MidLoc2  = CellSet[i_cell].MidLoc2 + 2.0/(pow(2.0,double(i_depth+1)));
    }
  }

  // Compute Gauss points
  // Loop over all cells
  for (size_t i_cell = 0; i_cell < CellSet.size(); i_cell++)
    for (size_t jGP = 0; jGP < LocGPxi.size(); jGP++)
      for (size_t iGP = 0; iGP < LocGPxi.size(); iGP++) {

	// Compute global x,y-coordinates of current Gauss point
	double GlobGP1 = CellSet[i_cell].MidGlob1 + LocGPxi[iGP]*(hx/pow(2.0,double(CellSet[i_cell].depth+1)));
	double GlobGP2 = CellSet[i_cell].MidGlob2 + LocGPxi[jGP]*(hy/pow(2.0,double(CellSet[i_cell].depth+1)));

	// Do inside-outside test, if Gauss point is inside -> append to list
	if (FindAlpha(GlobGP1,GlobGP2) > 0.0) {
	  GP1.push_back(CellSet[i_cell].MidLoc1 + LocGPxi[iGP]*(2.0/(pow(2.0,double(CellSet[i_cell].depth+1)))));
	  GP2.push_back(CellSet[i_cell].MidLoc2 + LocGPxi[jGP]*(2.0/(pow(2.0,double(CellSet[i_cell].depth+1)))));
	  GPw.push_back(LocGPw[iGP]*LocGPw[jGP]*1.0/pow(2.0,double(CellSet[i_cell].depth))*
			1.0/pow(2.0,double(CellSet[i_cell].depth)));
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

bool Immersed::getQuadraturePoints (double x1, double y1, double z1,
				    double x2, double y2, double z2,
				    double x3, double y3, double z3,
				    double x4, double y4, double z4,
				    double x5, double y5, double z5,
				    double x6, double y6, double z6,
				    double x7, double y7, double z7,
				    double x8, double y8, double z8,
				    int max_depth, int p,
				    RealArray& GP1,
				    RealArray& GP2,
				    RealArray& GP3,
				    RealArray& GPw)
{
  return false; // TODO 3D version...
}


// Wrapper for processing multiple elements.

bool Immersed::getQuadraturePoints (const Real3DMat& elmCorner,
				    int max_depth, int p,
				    Real3DMat& quadPoints)
{
  int nsd = 0;
  bool ok = true;
  quadPoints.resize(elmCorner.size());
  for (size_t e = 0; e < elmCorner.size() && ok; e++)
  {
    RealArray GP[4];
    const Real2DMat& X = elmCorner[e];
    switch (X.size()) {
    case 4: // 2D element
      nsd = 2;
      ok = getQuadraturePoints(X[0][0],X[0][1],
			       X[1][0],X[1][1],
			       X[2][0],X[2][1],
			       X[3][0],X[3][1],max_depth,p,
			       GP[1],GP[2],GP[0]);
      break;
    case 8: // 3D element
      nsd = 3;
      ok = getQuadraturePoints(X[0][0],X[0][1],X[0][2],
			       X[1][0],X[1][1],X[1][2],
			       X[2][0],X[2][1],X[2][2],
			       X[3][0],X[3][1],X[3][2],
			       X[4][0],X[4][1],X[4][2],
			       X[5][0],X[5][1],X[5][2],
			       X[6][0],X[6][1],X[6][2],
			       X[7][0],X[7][1],X[7][2],max_depth,p,
			       GP[1],GP[2],GP[3],GP[0]);
      break;
    default:
      std::cerr <<" *** Immersed::getQuadraturePoints: Invalid element ("
		<< X.size() <<" corners)."<< std::endl;
      nsd = 0;
      ok = false;
    }

    // Store the quadrature points
    RealArray xg(nsd+1);
    for (size_t i = 0; i < GP[0].size(); i++)
    {
      for (int d = 0; d < nsd; d++)
	xg[d] = GP[d+1][i];
      xg[nsd] = GP[0][i];
      quadPoints[e].push_back(xg);
    }
  }

  return ok;
}
