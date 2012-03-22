#include "GoTools/geometry/SurfaceInterpolator.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"

#include "ASMs2D.h"
#include "ASMstruct.h"
#include "SplineInterpolator.h"

#include <vector>


Go::SplineSurface*
  ASMs2D::leastsquare_approximation(const Go::BsplineBasis& basis_u,
                                    const Go::BsplineBasis& basis_v,
                                    std::vector<double>& par_u,
                                    std::vector<double>& par_v,
                                    std::vector<double>& wpar_u,
                                    std::vector<double>& wpar_v,
                                    std::vector<double>& points,
                                    int dimension,
                                    bool rational,
                                    std::vector<double>& weights) const
{
  // Check input
  ASSERT(par_u.size()*par_v.size() == points.size()/dimension);
  ASSERT((int)wpar_u.size() == (int)par_u.size());
  ASSERT((int)wpar_v.size() == (int)par_v.size());
  std::vector<double> points2;
  int perknot=dimension;
  if (rational)
  {
    // Include weight information
    // First get weights in interpolation points
    shared_ptr<Go::SplineSurface> denom =
      shared_ptr<Go::SplineSurface>(new Go::SplineSurface(basis_u, basis_v,
                                                          weights.begin(), 1,
                                                          false));
    std::vector<double> wgtval;
    denom->gridEvaluator(wgtval, par_u, par_v);
    size_t nmb_pnt = par_u.size()*par_v.size();
    points2.reserve(nmb_pnt*(dimension+1));
    for (size_t kr=0; kr<nmb_pnt; ++kr)
    {
      for (int kh=0; kh<dimension; kh++)
        points2.push_back(points[kr*dimension+kh]*wgtval[kr]);
      points2.push_back(wgtval[kr]);
    }
    perknot++;
  }
  else
    points2 = points;


  // Interpolate curves in the first parameter direction
  size_t ki;
  std::vector<double> cv_coefs;
  std::vector<double> tg_pnt;
  for (ki=0; ki<par_v.size(); ++ki)
  {
    // Interpolate
    std::vector<double> coefs;
    std::vector<double> pnts;
    pnts.insert(pnts.end(), points2.begin()+ki*perknot*par_u.size(),
                points2.begin()+(ki+1)*perknot*par_u.size());
    SplineInterpolator::leastsquare_approximation(par_u,  wpar_u, pnts,
                                                  tg_pnt, basis_u, coefs);
    cv_coefs.insert(cv_coefs.end(), coefs.begin(), coefs.end());
  }

  // Interpolate the curves to make a surface
  std::vector<double> sf_coefs;
  SplineInterpolator::leastsquare_approximation(par_v, wpar_v, cv_coefs,
                                                tg_pnt, basis_v, sf_coefs);

  // Make surface
  Go::SplineSurface* surf = new Go::SplineSurface(basis_u, basis_v,
                                                  sf_coefs.begin(), dimension,
                                                  rational);

  return surf;
}


struct InterpolateLoop {
  // for bound : bp*p + bm*m + b*1
  double bp, bm, b;
  // start offset for parameters; u1j*j+u1c*count+u1p*p+u1m*m+u1
  int u1j, u1c, u1p, u1m, u1;
  // end offset for parameters; u2j*j+u2c*count+u2p*p+u2m*m+u2
  int u2j, u2c, u2p, u2m, u2;
  // start offset for points; (p1j*j+p1c*count+p1p*p+p1m*m+p1)*perknot
  int p1j, p1c, p1p, p1m, p1;
  // end offset for points; (p2j*j+p2c*count+p2p*p+p2m*m+p2)*perknot
  int p2j, p2c, p2p, p2m, p2;
};


// helper tables with coefficients, see InterpolateLoop struct doxy
static InterpolateLoop equal[] =
  {{0.5f, -0.5f,  1.f,
       2,     2,    0,  0,  0,
       2,     2,    2,  0, -1,
       2,     2,    0,  0,  0,
       2,     2,    2,  0, -1},
   {0.0f,   2.f, -3.f,
       1,     2,    1, -1,  1,
       1,     2,    3, -1,  0,
       1,     2,    1, -1,  1,
       1,     2,    3, -1,  0},
   {0.5f, -0.5f,  1.f,
       2,     2,    1,  1, -2,
       2,     2,    3,  1, -3,
       2,     2,    1,  1, -2,
       2,     2,    3,  1, -3}};

static InterpolateLoop differs[] =
  {{0.5f, -0.5f,  0.5f,
       2,     2,     0,  0,  0,
       2,     2,     2,  0, -1,
       2,     2,     0,  0,  0,
       2,     2,     2,  0, -1},
   {0.0f,  1.0f, -1.0f,
       1,     2,     1, -1,  0,
       1,     2,     3, -1, -1,
       1,     2,     1, -1,  0,
       1,     2,     3, -1, -1},
   {0.0f,  1.0f, -1.0f,
       1,     2,     1,  0,  0,
       1,     2,     3,  0, -1,
       1,     2,     1,  0,  0,
       1,     2,     3,  0, -1},
   {0.5f, -0.5f,  0.5f,
       2,     2,     1,  1, -1,
       2,     2,     3,  1, -2,
       2,     2,     1,  1, -1,
       2,     2,     3,  1, -2}};

static InterpolateLoop no_multiple_skeleton = 
  {0.f, 0.f, 0.f,
     2,   0, 0, 0,  0,
     2,   0, 2, 0, -1,
     2,   0, 0, 0,  0,
     2,   0, 2, 0, -1};


static void quasiInterpolateLoop(const InterpolateLoop& loop,
                                 int& countj,
                                 const std::vector<double>& par_u,
                                 const std::vector<double>& pnts,
                                 std::vector<double>& cv_coefs,
                                 int count, int p, int m, int perknot,
                                 const Go::BsplineBasis& basis)
{
  std::vector<double> tg_pnt;
  int bnd = loop.bp*p+loop.bm*m+loop.b;
  int sbofs_base = loop.u1c*count+loop.u1p*p+loop.u1m*m+loop.u1;
  int ebofs_base = loop.u2c*count+loop.u2p*p+loop.u2m*m+loop.u2;
  int dbofs_base = loop.p1c*count+loop.p1p*p+loop.p1m*m+loop.p1;
  int deofs_base = loop.p2c*count+loop.p2p*p+loop.p2m*m+loop.p2;
  for (int j=0;j<bnd;++j) {
    std::vector<double> pnts_parts;
    std::vector<double> par_u_parts;

    int sofs = loop.u1j*j+sbofs_base;
    int eofs = loop.u2j*j+ebofs_base;
    par_u_parts.insert(par_u_parts.end(),
                       par_u.begin()+sofs,
                       par_u.begin()+eofs);

    sofs = loop.p1j*j+dbofs_base;
    eofs = loop.p2j*j+deofs_base;
    pnts_parts.insert(pnts_parts.end(),
                      pnts.begin()+sofs*perknot,
                      pnts.begin()+eofs*perknot);
    std::vector<double> coefs;
    SplineInterpolator::quasiinterpolate(par_u_parts, pnts_parts,
                                         tg_pnt, basis,
                                         countj++, coefs);
    cv_coefs.insert(cv_coefs.end(), coefs.begin(), coefs.end());
  }
}


static void quasiInterpolate1D(int count_multiple_knots,
                               const std::vector<double>& par,
                               const std::vector<double>& pnts,
                               const Go::BsplineBasis& basis,
                               int p, int perknot,
                               double multi_value,
                               std::vector<double>& coefs)
{
  int end = par.size() + 2*count_multiple_knots - basis.numCoefs()+1;
  // case evaluation depending on multiplicity of knots
  if (count_multiple_knots == 0) {
    InterpolateLoop no_multiple(no_multiple_skeleton);
    no_multiple.b = end;
    int countj=0, count=0, m=0;
    quasiInterpolateLoop(no_multiple, countj, par, pnts,
                         coefs, count, p, m, perknot, basis);
  }
  else if (p > 1) {
    int count = 0;
    int dcount = 0;
    int ti = 0;
    int terminate = end;
    int m = count_multiple_knots+1;

    while (ti < terminate) {
      std::vector<double> par_parts;
      par_parts.insert(par_parts.end(),par.begin()+ti*2,
                       par.begin()+ti*2+(2*p-1));

      if (par_parts.back() == multi_value) {
        terminate = 2*ti+1;

        InterpolateLoop* loop;
        int loops;
        if ((!(p&1) && !(m&1)) || ( (p&1) && (m&1))) {
          loop = equal;
          loops = 3;
        } else {
          loop = differs;
          loops = 4;
        }

        int countj=count;
        for (int i=0;i<loops;++i)
          quasiInterpolateLoop(loop[i], countj, par, pnts,
                               coefs, count, p, m, perknot, basis);
        dcount = p+(m-2);
      } else {
        InterpolateLoop loop = {0};
        loop.b = 1;
        int pos=ti;
        if (dcount > 0) {
          loop.u1 = (dcount+ti-(m-1))*2;
          loop.u2 = (dcount+ti-(m-1))*2+(2*p-1);
          loop.p1 = (dcount+ti-(m-1))*2;
          loop.p2 = ((dcount+ti-(m-1))*2+(2*p-1));
          pos += dcount;
        } else {
          loop.u1 = (ti+dcount)*2;
          loop.u2 = (ti+dcount)*2+(2*p-1);
          loop.p1 = (ti+dcount)*2;
          loop.p2 = ((ti+dcount)*2+(2*p-1));
          count++;
        }
        quasiInterpolateLoop(loop, pos, par, pnts,
                             coefs, count, p, m, perknot, basis);
      }
      ti++;
    }
  }// end if p>1
}


Go::SplineSurface*
  ASMs2D::quasiInterpolation(const Go::BsplineBasis& basis_u,
                             const Go::BsplineBasis& basis_v,
                             std::vector<double>& par_u,
                             std::vector<double>& par_v,
                             std::vector<double>& points,
                             int dimension,
                             bool rational,
                             std::vector<double>& weights) const
{
  // Check knot multiplicity for case evaluation
  std::vector<double> knots_simple_u;
  basis_u.knotsSimple(knots_simple_u);
  std::vector<double> knots_simple_v;
  basis_v.knotsSimple(knots_simple_v);
  int count_multipl_knots_u;
  int count_multipl_knots_v;
  count_multipl_knots_u = basis_u.numCoefs()-basis_u.order()+2 - knots_simple_u.size();
  count_multipl_knots_v = basis_v.numCoefs()-basis_v.order()+2 - knots_simple_v.size();

  // Check input
  ASSERT(par_u.size()*par_v.size() == points.size()/dimension);
  ASSERT(2*(basis_u.numCoefs()-basis_u.order()+1)+1-2*count_multipl_knots_u == (int)par_u.size());
  ASSERT(2*(basis_v.numCoefs()-basis_v.order()+1)+1-2*count_multipl_knots_v == (int)par_v.size());
  ASSERT(2*basis_u.order()-3 <= (int)par_u.size());
  ASSERT(2*basis_v.order()-3 <= (int)par_v.size());
  if (count_multipl_knots_u > 0)
    ASSERT( (par_u.size()+1)*0.5 - (2*(basis_u.order()-1)-1) >= 0);
  if (count_multipl_knots_v > 0)
    ASSERT( (par_v.size()+1)*0.5 - (2*(basis_v.order()-1)-1) >= 0);

  std::vector<double> points2;
  int perknot=dimension;
  if (rational)
  {
    // Include weight information
    // First get weights in interpolation points
    shared_ptr<Go::SplineSurface> denom =
      shared_ptr<Go::SplineSurface>(new Go::SplineSurface(basis_u, basis_v,
                                                          weights.begin(), 1,
                                                          false));
    std::vector<double> wgtval;
    denom->gridEvaluator(wgtval, par_u, par_v);
    size_t nmb_pnt = par_u.size()*par_v.size();
    points2.reserve(nmb_pnt*(dimension+1));
    for (size_t kr=0; kr<nmb_pnt; ++kr)
    {
      for (int kh=0; kh<dimension; kh++)
        points2.push_back(points[kr*dimension+kh]*wgtval[kr]);
      points2.push_back(wgtval[kr]);
    }
    perknot++;
  }
  else
    points2 = points;

  // interpolate in the first direction
  size_t ki;
  int p,q;
  p = basis_u.order()-1;
  q = basis_v.order()-1;

  std::vector<double> cv_coefs;
  std::vector<double> tg_pnt;
  std::vector<double> multi_u(knots_simple_u.size());
  std::vector<double> multi_idx_u(knots_simple_u.size());
  double multi_value_u=-1.0;
  multi_idx_u[0]=0;
  multi_u[0]=1;
  multi_idx_u[knots_simple_u.size()-1]=0;
  multi_u[knots_simple_u.size()-1]=1;

  // create multiplicity index for case evaluation
  for (int i = 1;i<(int)knots_simple_u.size()-1;i++)
  {
    int knot_intern_multipl;
    knot_intern_multipl = basis_u.knotMultiplicity(knots_simple_u[i]);
    multi_u[i]=knot_intern_multipl;
    if (knot_intern_multipl > 1) {
      multi_idx_u[i]=1;
      multi_value_u = knots_simple_u[i];
    } else
      multi_idx_u[i]= 0;
  }

  //----------------------------------------------------


  for (ki=0; ki<par_v.size(); ++ki)
  {
    // Interpolate over local patch of size 2*degree_p-1 (moving window priciple)
    std::vector<double> pnts;
    pnts.insert(pnts.end(), points2.begin()+ki*perknot*par_u.size(),
                points2.begin()+(ki+1)*perknot*par_u.size());
    quasiInterpolate1D(count_multipl_knots_u, par_u, pnts,
                       basis_u, p, perknot, multi_value_u, cv_coefs);
  }

  //-------------------------------------------------------------------------
  // Interpolate the curves to make a surface
  // same principle as in u-direction
  //----------------------------------------------------
  //
  std::vector<double> multi_v(knots_simple_v.size());
  std::vector<double> multi_idx_v(knots_simple_v.size());
  double multi_value_v=-1.0;
  multi_idx_v[0]=0;
  multi_v[0]=1;
  multi_idx_v[knots_simple_v.size()-1]=0;
  multi_v[knots_simple_v.size()-1]=1;
  for (int i = 1;i<(int)knots_simple_v.size()-1;i++)
  {
    int knot_intern_multipl;
    knot_intern_multipl = basis_v.knotMultiplicity(knots_simple_v[i]);
    multi_v[i]=knot_intern_multipl;
    if (knot_intern_multipl > 1) {  
      multi_idx_v[i]=1;
      multi_value_v = knots_simple_v[i];
    } else
      multi_idx_v[i]= 0;
  }
  //----------------------------------------------------

  std::vector<double> sf_coefs;
  int ucount = (par_u.size() + 2*count_multipl_knots_u - basis_u.numCoefs()+1)*(2*p-1)*perknot;
  quasiInterpolate1D(count_multipl_knots_v, par_v, cv_coefs,
                     basis_v, q, ucount, multi_value_v, sf_coefs);

  //eval sf_coefs to generate sf_coefs_solution
  //boundary conditions
  int n = basis_u.numCoefs();
  int m = basis_v.numCoefs();
  int spxl = 2*p-1;
  int spyl = 2*q-1;
  int qpxl = 2*(n-p)+1;
  int qpyl = 2*(m-q)+1;
  int invxl = qpxl-n+1;
  int invyl = qpyl-m+1;
  int gmxl = spxl*invxl;
  int gmyl = spyl*invyl;

  std::vector<std::vector<double> > gm;
  std::vector<double> tmp;
  gm.resize(gmyl);
  for(int i = 0; i < gmyl ; i++)
    gm[i].resize(gmxl*perknot);

  std::vector<double> sf_coefs_solution;
  sf_coefs_solution.resize(n*m*perknot);

  if (invyl == 1 && invxl == 1)
  {
    sf_coefs_solution = sf_coefs;
  }
  else
  {
    for (int b = 0; b < invyl; b++)
    {
      for (int a = 0; a < invxl; a++)
      {
        int starti = b*(2*q-1);
        int startj = a*(2*p-1);
        int endi, endj;
        if (invyl == 1) {
          if (a == 0) {
            endi = starti+2*q-1;
            endj = (startj+p)*perknot;
          } else if (a == invxl-1) {
            endi = starti+2*q-1;
            endj = (startj+2*p-1)*perknot;
            startj = (startj+p-1)*perknot;
          }
          else {
            endi = (starti+2*q-1);
            endj = (startj+p)*perknot;
            startj = (startj+p-1)*perknot;
          }
        }
        else if (invxl == 1) {
          if (b == 0) {
            endi = starti+q;
            endj = (startj+2*p-1)*perknot;
          } else if (b == invyl-1) {
            endi = starti+2*q-1;
            starti = starti+q-1;
            endj = (startj+2*p-1)*perknot;
          } else {
            endi = starti+q;
            starti = starti+q-1;
            endj = (startj+2*p-1)*perknot;
          }
        }
        else if (b==0 && a == 0) { //top left corner
          endi = starti+q;
          endj = (startj+p)*perknot;
        } else if (b == 0 && a == invxl-1) { //top right corner
          endi = starti+q;
          endj = (startj+2*p-1)*perknot;
          startj = (startj+p-1)*perknot;
        } else if (b == invyl-1 && a == 0) { //bottom left corner
          endi = (starti+2*q-1);
          starti = (starti+q-1);
          endj = (startj+p)*perknot;
        } else if (b == invyl-1 && a == invxl-1) { //bottom right corner
          endi = (starti+2*q-1);
          starti = (starti+q-1);
          endj = (startj+2*p-1)*perknot;
          startj = (startj+p-1)*perknot;
        } else if (b==0) { // top edge
          endi = starti+q;
          endj = (startj+p)*perknot;
          startj = (startj+p-1)*perknot;
        } else if (a == invxl-1) { // right edge
          endi = starti+q;
          starti = starti+q-1;
          endj = (startj+2*p-1)*perknot;
          startj = (startj+p-1)*perknot;
        } else if (a == 0) { // left edge
          endi = starti+q;
          starti = starti+q-1;
          endj = (startj+p)*perknot;
        } else if (b == invyl-1) {// bottom edge
          endi = starti+2*q-1;
          starti = starti+q-1;
          endj = (startj+p)*perknot;
          startj = (startj+p-1)*perknot;
        } else {
          // interior elements
          endi = starti+q;
          starti = starti+q-1;
          endj = (startj+p)*perknot;
          startj = (startj+p-1)*perknot;
        }
        for (int i = starti;i<endi;i++)
          for (int j = startj;j<endj;j++)
            gm[i][j]= 1;
      }
    }


    int count = -1;
    for (int i = 0;i<gmyl;i++)
      for (int j = 0;j<gmxl*perknot;j++)
        if (gm[i][j] == 1) {
          count = count+1;
          sf_coefs_solution[count] = sf_coefs[i*(gmxl*perknot)+j];
        }
  }

  ///////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////


  // Make surface
  Go::SplineSurface* surf = new Go::SplineSurface(basis_u, basis_v,
                                                  sf_coefs_solution.begin(),
                                                  dimension, rational);

  return surf;
}


Go::SplineSurface*
  ASMs2D::VariationDiminishingSplineApproximation(const Go::BsplineBasis& basis_u,
                                                  const Go::BsplineBasis& basis_v,
                                                  std::vector<double>& par_u,
                                                  std::vector<double>& par_v,
                                                  std::vector<double>& points,
                                                  int dimension,
                                                  bool rational,
                                                  std::vector<double>& weights) const
{
  // Check input
  ASSERT(par_u.size()*par_v.size() == points.size()/dimension);
  ASSERT(basis_u.numCoefs() == (int)par_u.size());
  ASSERT(basis_v.numCoefs() == (int)par_v.size());

  size_t sizepoints;
  if (rational)
  {
    sizepoints = par_u.size()*par_v.size();
  }
  else
    sizepoints = points.size();

  std::vector<double> local_coefs;
  size_t countpoints=0;
  if (rational)
  {
    for (size_t i = 0; i<sizepoints;i++)
    {
      for (int j = 0; j<dimension;j++){
        local_coefs.push_back(points[countpoints]*weights[i]);
        countpoints++;}
        local_coefs.push_back(weights[i]);
    }
  }
  else
    local_coefs = points;


  // Make surface
  Go::SplineSurface* surf = new Go::SplineSurface(basis_u, basis_v,
                                                  local_coefs.begin(),
                                                  dimension, rational);

  return surf;
}
