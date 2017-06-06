#include "GoTools/trivariate/SplineVolume.h"
#include "SplineInterpolator.h"
#include "MatVec.h"

/*!
  \brief Global projection method (Least-Square Approximation).
  \param[in] basis_u Basis values in the first parameter direction
  \param[in] basis_v Basis values in the second parameter direction
  \param[in] basis_w Basis values in the third parameter direction
  \param[in] par_u Gauss values in the first parameter direction
  \param[in] par_v Gauss values in the second parameter direction
  \param[in] par_w Gauss values in the third parameter direction
  \param[in] wpar_u Gauss weights in the first parameter direction
  \param[in] wpar_v Gauss weights in the second parameter direction
  \param[in] wpar_w Gauss weights in the third parameter direction
  \param[in] points Secondary solution field evaluated at Gauss points
  \param[in] dimension Dimension of the secondary solution field
  \param[in] rational Value marks NURBS geometry
  \param[in] weights NURBS weights for the projective control points
  \return Spline volume object representing the projected field
*/

static Go::SplineVolume*
leastsquare_approximation(const Go::BsplineBasis& basis_u,
                          const Go::BsplineBasis& basis_v,
                          const Go::BsplineBasis& basis_w,
                          const RealArray& par_u,
                          const RealArray& par_v,
                          const RealArray& par_w,
                          const RealArray& wpar_u,
                          const RealArray& wpar_v,
                          const RealArray& wpar_w,
                          const RealArray& points,
                          int dimension, bool rational,
                          const RealArray& weights)
{
  // Check input
  ASSERT(par_u.size()*par_v.size()*par_w.size() == points.size()/dimension);
  ASSERT((int)wpar_u.size() == (int)par_u.size());
  ASSERT((int)wpar_v.size() == (int)par_v.size());
  ASSERT((int)wpar_w.size() == (int)par_w.size());

  std::vector<double> points2;
  int perknot=dimension;
  if (rational)
  {
    Go::SplineVolume denom(basis_u, basis_v, basis_w, weights.begin(), 1,false);
    std::vector<double> wgtval;
    denom.gridEvaluator(par_u, par_v, par_w, wgtval);
    size_t nmb_pnt = par_u.size()*par_v.size()*par_w.size();
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

  // Interpolate surfaces in the second parameter direction and
  // curves in the first parameter direction
  size_t ki, kj;
  std::vector<double> sf_coefs;
  std::vector<double> tg_pnt;
  for (kj=0; kj<par_w.size(); ++kj)
  {
    std::vector<double> cv_coefs;
    for (ki=0; ki<par_v.size(); ++ki)
    {
      // Interpolate
      std::vector<double> coefs;
      std::vector<double> pnts;
      pnts.insert(pnts.end(),
          points2.begin()+(kj*par_v.size()+ki)*perknot*par_u.size(),
          points2.begin()+(kj*par_v.size()+ki+1)*perknot*par_u.size());

      SplineInterpolator::leastsquare_approximation(par_u, wpar_u, pnts, tg_pnt,
                                                    basis_u, coefs);
      cv_coefs.insert(cv_coefs.end(), coefs.begin(), coefs.end());
    }

    // Interpolate the curves to make a surface
    std::vector<double> coefs2;

    SplineInterpolator::leastsquare_approximation(par_v, wpar_v, cv_coefs, tg_pnt,
                                                  basis_v, coefs2);
    sf_coefs.insert(sf_coefs.end(), coefs2.begin(), coefs2.end());
  }

  // Interpolate surfaces to create volume
  std::vector<double> vol_coefs;
  SplineInterpolator::leastsquare_approximation(par_w, wpar_w, sf_coefs, tg_pnt,
                                                basis_w, vol_coefs);

  // Make volume
  return new Go::SplineVolume(basis_u, basis_v, basis_w, vol_coefs.begin(),
                              dimension, rational);
}


/*!
  \brief Local projection method (Quasi-Interpolation).
  \param[in] basis_u Basis values in the first parameter direction
  \param[in] basis_v Basis values in the second parameter direction
  \param[in] basis_w Basis values in the third parameter direction
  \param[in] par_u Grevielle sites in the first parameter direction
  \param[in] par_v Grevielle sites in the second parameter direction
  \param[in] par_w Grevielle sites in the third parameter direction
  \param[in] points Secondary solution field evaluated at Greville points
  \param[in] dimension Dimension of the secondary solution field
  \param[in] rational Value marks NURBS geometry
  \param[in] weights NURBS weights for the projective control points
  \return Spline volume object representing the projected field
*/

static Go::SplineVolume*
quasiInterpolation(const Go::BsplineBasis& basis_u,
                   const Go::BsplineBasis& basis_v,
                   const Go::BsplineBasis& basis_w,
                   const RealArray& par_u,
                   const RealArray& par_v,
                   const RealArray& par_w,
                   const RealArray& points,
                   int dimension, bool rational,
                   const RealArray& weights)
{
  std::vector< double > knots_simple_u;
  basis_u.knotsSimple(knots_simple_u);
  std::vector< double > knots_simple_v;
  basis_v.knotsSimple(knots_simple_v);

  int  count_multipl_knots_u;
  int  count_multipl_knots_v;

  count_multipl_knots_u = basis_u.numCoefs()-basis_u.order()+2 - knots_simple_u.size();
  count_multipl_knots_v = basis_v.numCoefs()-basis_v.order()+2 - knots_simple_v.size();

  ASSERT(2*basis_u.order()-3 <= (int)par_u.size());
  ASSERT(2*basis_v.order()-3 <= (int)par_v.size());

  // Check input
  ASSERT(par_u.size()*par_v.size()*par_w.size() == points.size()/dimension);
  ASSERT(2*(basis_u.numCoefs()-basis_u.order()+1)+1-2*count_multipl_knots_u == (int)par_u.size());
  ASSERT(2*(basis_v.numCoefs()-basis_v.order()+1)+1-2*count_multipl_knots_v == (int)par_v.size());

  if (count_multipl_knots_u > 0)
    ASSERT( (par_u.size()+1)*0.5 - (2*(basis_u.order()-1)-1) >= 0);
  if (count_multipl_knots_v > 0)
    ASSERT( (par_v.size()+1)*0.5 - (2*(basis_v.order()-1)-1) >= 0);

  std::vector<double> points2;
  int perknot=dimension;
  if (rational)
  {
    Go::SplineVolume denom(basis_u, basis_v, basis_w, weights.begin(), 1,false);
    std::vector<double> wgtval;
    denom.gridEvaluator(par_u, par_v, par_w, wgtval);
    size_t nmb_pnt = par_u.size()*par_v.size()*par_w.size();
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


  // Interpolate surfaces in the second parameter direction and
  // curves in the first parameter direction
  size_t ki, kj;
  int ui;
  std::vector<double> volinput_coefs;
  std::vector<double> tg_pnt;
  int p,q;
  p = basis_u.order()-1;
  q = basis_v.order()-1;

  for (kj=0; kj<par_w.size(); ++kj)
  {
    std::vector<double> cv_coefs;
    std::vector<double> multi_u(knots_simple_u.size());
    std::vector<double> multi_idx_u(knots_simple_u.size());
    double multi_value_u=-1.0;
    multi_idx_u[0]=0;
    multi_u[0]=1;
    multi_idx_u[knots_simple_u.size()-1]=0;
    multi_u[knots_simple_u.size()-1]=1;
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

    for (ki=0; ki<par_v.size(); ++ki)
    {
      // Interpolate
      std::vector<double> pnts;

      pnts.insert(pnts.end(), points2.begin()+(kj*par_v.size()+ki)*perknot*par_u.size(),points2.begin()+(kj*par_v.size()+ki+1)*perknot*par_u.size());
      int ui_end =  par_u.size()+2*count_multipl_knots_u-basis_u.numCoefs()+1;

      if (count_multipl_knots_u == 0) {
        for (ui = 0; ui< ui_end;ui++) {
          std::vector<double> pnts_parts;
          std::vector<double> par_u_parts;
          par_u_parts.insert(par_u_parts.end(),par_u.begin()+ui*2,par_u.begin()+ui*2+(2*p-1));
          pnts_parts.insert(pnts_parts.end(), pnts.begin()+ui*2*perknot, pnts.begin()+(ui*2+(2*p-1))*perknot);
          std::vector<double> coefs;

          SplineInterpolator::quasiinterpolate(par_u_parts, pnts_parts, tg_pnt,
                                               basis_u, ui, coefs);
          cv_coefs.insert(cv_coefs.end(), coefs.begin(), coefs.end());
        }
      }
      else
      {
        if (p>1) {
          int count = 0;
          int dcount = 0;
          int count_intern = 0;
          int ti = 0;
          int terminate=ui_end;
          int m = count_multipl_knots_u+1;

          while(ti<terminate)
          {
            std::vector<double> par_u_parts;
            par_u_parts.insert(par_u_parts.end(),par_u.begin()+ti*2,par_u.begin()+ti*2+(2*p-1));

            if ( par_u_parts[par_u_parts.size()-1]==multi_value_u)
            {
              terminate = 2*ti+1;

              if ( (!(p&1) && !(m&1)) || ( (p&1) && (m&1) ) )
              {
                int countj=0;
                for (int j= 0; j < (p-m+2)*0.5;j++){
                  std::vector<double> pnts_parts;
                  std::vector<double> par_u_parts;
                  par_u_parts.insert(par_u_parts.end(),par_u.begin()+(j+count)*2,par_u.begin()+(j+count)*2+(2*p-1));//dp(j+count,:)
                  pnts_parts.insert(pnts_parts.end(), pnts.begin()+(j+count)*2*perknot, pnts.begin()+((j+count)*2+(2*p-1))*perknot);
                  std::vector<double> coefs;
                  SplineInterpolator::quasiinterpolate(par_u_parts, pnts_parts,
                                                       tg_pnt, basis_u, 
                                                       j+count+countj, coefs);
                  cv_coefs.insert(cv_coefs.end(), coefs.begin(), coefs.end());
                  count_intern = j+count+countj;
                }
                countj = count_intern+1;
                for (int j= 0; j < 2*m-3;j++){
                  std::vector<double> pnts_parts;
                  std::vector<double> par_u_parts;
                  par_u_parts.insert(par_u_parts.end(),par_u.begin()+p-m+1+j+2*count,par_u.begin()+3*p-m-1+j+2*count+1);//dist_vec(p-m+1+j+2*count:3*p-m-1+j+2*count)
                  pnts_parts.insert(pnts_parts.end(), pnts.begin()+(p-m+1+j+2*count)*perknot, pnts.begin()+(3*p-m-1+j+2*count+1)*perknot);
                  std::vector<double> coefs;

                  SplineInterpolator::quasiinterpolate(par_u_parts, pnts_parts, tg_pnt, 
                                                       basis_u, j+countj, coefs);
                  cv_coefs.insert(cv_coefs.end(), coefs.begin(), coefs.end());
                  count_intern = j+countj;
                }
                countj = count_intern+1;
                for (int j = 0; j<(p-m+2)*0.5;j++){
                  std::vector<double> pnts_parts;
                  std::vector<double> par_u_parts;
                  par_u_parts.insert(par_u_parts.end(),par_u.begin()+(j+count+((p-m+2)*0.5)+m-2)*2,par_u.begin()+(j+count+((p-m+2)*0.5)+m-2)*2+(2*p-1));//dp(j+count+((p-m+2)*0.5)+m-2,:)
                  pnts_parts.insert(pnts_parts.end(), pnts.begin()+(j+count+((p-m+2)*0.5)+m-2)*2*perknot, pnts.begin()+((j+count+((p-m+2)*0.5)+m-2)*2+(2*p-1))*perknot);
                  std::vector<double> coefs;

                  SplineInterpolator::quasiinterpolate(par_u_parts, pnts_parts, tg_pnt, 
                                                       basis_u, j+countj, coefs);
                  cv_coefs.insert(cv_coefs.end(), coefs.begin(), coefs.end());
                  count_intern = j+countj;
                }
              }
              else
              {
                int countj=0;
                for (int j=0; j<(p-m+1)*0.5;j++){
                  std::vector<double> pnts_parts;
                  std::vector<double> par_u_parts;
                  par_u_parts.insert(par_u_parts.end(),par_u.begin()+(j+count)*2,par_u.begin()+(j+count)*2+(2*p-1));//dp(j+count,:)
                  pnts_parts.insert(pnts_parts.end(), pnts.begin()+(j+count)*2*perknot, pnts.begin()+((j+count)*2+(2*p-1))*perknot);
                  std::vector<double> coefs;

                  SplineInterpolator::quasiinterpolate(par_u_parts, pnts_parts, 
                                                       tg_pnt, basis_u,
                                                       j+count+countj, coefs);
                  cv_coefs.insert(cv_coefs.end(), coefs.begin(), coefs.end());
                  count_intern  = j+count+countj;
                }
                countj = count_intern+1;
                for (int j = 0;j<m-1;j++){
                  std::vector<double> pnts_parts;
                  std::vector<double> par_u_parts;
                  par_u_parts.insert(par_u_parts.end(),par_u.begin()+p-m+j+2*count,par_u.begin()+3*p-m-2+2*count+j+1);//dist_vec(p-m+j+2*count:3*p-m-2+2*count+j)
                  pnts_parts.insert(pnts_parts.end(), pnts.begin()+(p-m+j+2*count)*perknot, pnts.begin()+(3*p-m-2+2*count+j+1)*perknot);
                  std::vector<double> coefs;

                  SplineInterpolator::quasiinterpolate(par_u_parts, pnts_parts, 
                                                       tg_pnt, basis_u,
                                                       j+countj, coefs);
                  cv_coefs.insert(cv_coefs.end(), coefs.begin(), coefs.end());
                  count_intern = j+countj;
                }
                countj = count_intern+1;
                for (int j=0;j<m-1;j++){
                  std::vector<double> pnts_parts;
                  std::vector<double> par_u_parts;
                  par_u_parts.insert(par_u_parts.end(),par_u.begin()+p+j+2*count,par_u.begin()+3*p-2+j+2*count+1);//dist_vec(p+j+2*count:3*p-2+j+2*count)
                  pnts_parts.insert(pnts_parts.end(), pnts.begin()+(p+j+2*count)*perknot, pnts.begin()+(3*p-2+j+2*count+1)*perknot);
                  std::vector<double> coefs;

                  SplineInterpolator::quasiinterpolate(par_u_parts, pnts_parts,
                                                       tg_pnt, basis_u,
                                                       j+countj, coefs);
                  cv_coefs.insert(cv_coefs.end(), coefs.begin(), coefs.end());
                  count_intern = j+countj;
                }
                countj = count_intern+1;
                for (int j=0;j<(p-m+1)*0.5;j++){
                  std::vector<double> pnts_parts;
                  std::vector<double> par_u_parts;
                  par_u_parts.insert(par_u_parts.end(),par_u.begin()+(j+count+((p-m+1)*0.5)+m-1)*2,par_u.begin()+(j+count+((p-m+1)*0.5)+m-1)*2+(2*p-1));//dp(j+count+((p-m+1)*0.5)+m-1,:)
                  pnts_parts.insert(pnts_parts.end(), pnts.begin()+(j+count+((p-m+1)*0.5)+m-1)*2*perknot, pnts.begin()+((j+count+((p-m+1)*0.5)+m-1)*2+(2*p-1))*perknot);
                  std::vector<double> coefs;

                  SplineInterpolator::quasiinterpolate(par_u_parts, pnts_parts,
                                                       tg_pnt, basis_u,
                                                       j+countj, coefs);
                  cv_coefs.insert(cv_coefs.end(), coefs.begin(), coefs.end());
                }			    
              }//else
              dcount = p+(m-2);
            }//end if multivalue 
            else
            {
              if (dcount > 0)
              {
                std::vector<double> pnts_parts;
                std::vector<double> par_u_parts;
                par_u_parts.insert(par_u_parts.end(),par_u.begin()+(dcount+ti-(m-1))*2,par_u.begin()+(dcount+ti-(m-1))*2+(2*p-1));
                pnts_parts.insert(pnts_parts.end(), pnts.begin()+(dcount+ti-(m-1))*2*perknot, pnts.begin()+((dcount+ti-(m-1))*2+(2*p-1))*perknot);

                std::vector<double> coefs;

                SplineInterpolator::quasiinterpolate(par_u_parts, pnts_parts, 
                                                      tg_pnt, basis_u,
                                                      ti+dcount, coefs);
                cv_coefs.insert(cv_coefs.end(), coefs.begin(), coefs.end());
              }
              else
              {
                std::vector<double> pnts_parts;
                std::vector<double> par_u_parts;
                par_u_parts.insert(par_u_parts.end(),par_u.begin()+(ti+dcount)*2,par_u.begin()+(ti+dcount)*2+(2*p-1));
                pnts_parts.insert(pnts_parts.end(), pnts.begin()+(ti+dcount)*2*perknot, pnts.begin()+((ti+dcount)*2+(2*p-1))*perknot);
                std::vector<double> coefs;

                SplineInterpolator::quasiinterpolate(par_u_parts, pnts_parts,
                                                     tg_pnt, basis_u, ti, coefs);
                cv_coefs.insert(cv_coefs.end(), coefs.begin(), coefs.end());

                count = count+1;
              }
            }
            ti++;
          }// end while
        }// end if p>1
      }//else

    }//end ki


    //-------------------------------------------------------------------------

    // Interpolate the curves to make a surface
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

    std::vector<double> sf_coefs;
    int ucount = (par_u.size() + 2*count_multipl_knots_u - basis_u.numCoefs()+1)*(2*p-1)*perknot;
    int vi_end =  par_v.size() + 2*count_multipl_knots_v - basis_v.numCoefs()+1 ;

    if (count_multipl_knots_v == 0)
    {
      for (int vi = 0; vi< vi_end;vi++)	
      {
        std::vector<double> coefs_parts;
        std::vector<double> par_v_parts;
        par_v_parts.insert(par_v_parts.end(),par_v.begin()+vi*2,par_v.begin()+vi*2+(2*q-1));
        coefs_parts.insert(coefs_parts.end(), cv_coefs.begin()+vi*2*ucount, cv_coefs.begin()+(vi*2+(2*q-1))*ucount);
        std::vector<double> sf_coefs_parts;

        SplineInterpolator::quasiinterpolate(par_v_parts, coefs_parts, tg_pnt,
                                             basis_v, vi, sf_coefs_parts);
        sf_coefs.insert(sf_coefs.end(), sf_coefs_parts.begin(), sf_coefs_parts.end());
      }
    }
    else
    {
      if (q>1)
      {
        int count = 0;
        int dcount = 0;
        int count_intern = 0;
        int ti = 0;
        int terminate=vi_end;
        int m = count_multipl_knots_v+1;

        while(ti<terminate)
        {
          std::vector<double> par_v_parts;
          par_v_parts.insert(par_v_parts.end(),par_v.begin()+ti*2,par_v.begin()+ti*2+(2*q-1));

          if ( par_v_parts[par_v_parts.size()-1]==multi_value_v)
          {
            terminate = 2*ti+1;

            if ( (!(q&1) && !(m&1)) || ( (q&1) && (m&1) ) )
            {	
              int countj=0;
              for (int j= 0; j < (q-m+2)*0.5;j++){
                std::vector<double> coefs_parts;
                std::vector<double> par_v_parts;
                par_v_parts.insert(par_v_parts.end(),par_v.begin()+(j+count)*2,par_v.begin()+(j+count)*2+(2*q-1));//dp(j+count,:)
                coefs_parts.insert(coefs_parts.end(), cv_coefs.begin()+(j+count)*2*ucount, cv_coefs.begin()+((j+count)*2+(2*q-1))*ucount);
                std::vector<double> sf_coefs_parts;

                SplineInterpolator::quasiinterpolate(par_v_parts, coefs_parts, 
                                                     tg_pnt, basis_v,
                                                     j+count+countj, sf_coefs_parts);
                sf_coefs.insert(sf_coefs.end(), sf_coefs_parts.begin(), sf_coefs_parts.end());
                count_intern = j+count+countj;
              }
              countj = count_intern+1;
              for (int j= 0; j < 2*m-3;j++){
                std::vector<double> coefs_parts;
                std::vector<double> par_v_parts;
                par_v_parts.insert(par_v_parts.end(),par_v.begin()+q-m+1+j+2*count,par_v.begin()+3*q-m-1+j+2*count+1);//dist_vec(q-m+1+j+2*count:3*q-m-1+j+2*count)
                coefs_parts.insert(coefs_parts.end(), cv_coefs.begin()+(q-m+1+j+2*count)*ucount, cv_coefs.begin()+(3*q-m-1+j+2*count+1)*ucount);
                std::vector<double> sf_coefs_parts;

                SplineInterpolator::quasiinterpolate(par_v_parts, coefs_parts, 
                                                     tg_pnt, basis_v,
                                                     j+countj, sf_coefs_parts);
                sf_coefs.insert(sf_coefs.end(), sf_coefs_parts.begin(), sf_coefs_parts.end());
                count_intern = j+countj;
              }
              countj = count_intern+1;
              for (int j = 0; j<(q-m+2)*0.5;j++){
                std::vector<double> coefs_parts;
                std::vector<double> par_v_parts;
                par_v_parts.insert(par_v_parts.end(),par_v.begin()+(j+count+((q-m+2)*0.5)+m-2)*2,par_v.begin()+(j+count+((q-m+2)*0.5)+m-2)*2+(2*q-1));//dp(j+count+((q-m+2)*0.5)+m-2,:)
                coefs_parts.insert(coefs_parts.end(), cv_coefs.begin()+(j+count+((q-m+2)*0.5)+m-2)*2*ucount, cv_coefs.begin()+((j+count+((q-m+2)*0.5)+m-2)*2+(2*q-1))*ucount);
                std::vector<double> sf_coefs_parts;

                SplineInterpolator::quasiinterpolate(par_v_parts, coefs_parts,
                                                     tg_pnt, basis_v,
                                                     j+countj, sf_coefs_parts);
                sf_coefs.insert(sf_coefs.end(), sf_coefs_parts.begin(), sf_coefs_parts.end());
                count_intern = j+countj;
              }
            }
            else
            {
              int countj=0;
              for (int j=0; j<(q-m+1)*0.5;j++){
                std::vector<double> coefs_parts;
                std::vector<double> par_v_parts;
                par_v_parts.insert(par_v_parts.end(),par_v.begin()+(j+count)*2,par_v.begin()+(j+count)*2+(2*q-1));//dp(j+count,:)
                coefs_parts.insert(coefs_parts.end(), cv_coefs.begin()+(j+count)*2*ucount, cv_coefs.begin()+((j+count)*2+(2*q-1))*ucount);

                std::vector<double> sf_coefs_parts;

                SplineInterpolator::quasiinterpolate(par_v_parts, coefs_parts,
                                                     tg_pnt, basis_v,
                                                     j+count+countj, sf_coefs_parts);
                sf_coefs.insert(sf_coefs.end(), sf_coefs_parts.begin(), sf_coefs_parts.end());
                count_intern = j+count+countj;
              }
              countj = count_intern+1;
              for (int j = 0;j<m-1;j++){
                std::vector<double> coefs_parts;
                std::vector<double> par_v_parts;
                par_v_parts.insert(par_v_parts.end(),par_v.begin()+q-m+j+2*count,par_v.begin()+3*q-m-2+2*count+j+1);//dist_vec(q-m+j+2*count:3*q-m-2+2*count+j)
                coefs_parts.insert(coefs_parts.end(), cv_coefs.begin()+(q-m+j+2*count)*ucount, cv_coefs.begin()+(3*q-m-2+2*count+j+1)*ucount);

                std::vector<double> sf_coefs_parts;

                SplineInterpolator::quasiinterpolate(par_v_parts, coefs_parts,
                                                     tg_pnt, basis_v,
                                                     j+countj, sf_coefs_parts);
                sf_coefs.insert(sf_coefs.end(), sf_coefs_parts.begin(), sf_coefs_parts.end());
                count_intern = j+countj;
              }
              countj = count_intern+1;
              for (int j=0;j<m-1;j++){
                std::vector<double> coefs_parts;
                std::vector<double> par_v_parts;
                par_v_parts.insert(par_v_parts.end(),par_v.begin()+q+j+2*count,par_v.begin()+3*q-2+j+2*count+1);//dist_vec(q+j+2*count:3*q-2+j+2*count)
                coefs_parts.insert(coefs_parts.end(), cv_coefs.begin()+(q+j+2*count)*ucount, cv_coefs.begin()+(3*q-2+j+2*count+1)*ucount);

                std::vector<double> sf_coefs_parts;

                SplineInterpolator::quasiinterpolate(par_v_parts, coefs_parts,
                                                     tg_pnt, basis_v,
                                                     j+countj, sf_coefs_parts);
                sf_coefs.insert(sf_coefs.end(), sf_coefs_parts.begin(), sf_coefs_parts.end());
                count_intern = j+countj;
              }
              countj = count_intern+1;
              for (int j=0;j<(q-m+1)*0.5;j++){
                std::vector<double> coefs_parts;
                std::vector<double> par_v_parts;
                par_v_parts.insert(par_v_parts.end(),par_v.begin()+(j+count+((q-m+1)*0.5)+m-1)*2,par_v.begin()+(j+count+((q-m+1)*0.5)+m-1)*2+(2*q-1));//dp(j+count+((q-m+1)*0.5)+m-1,:)
                coefs_parts.insert(coefs_parts.end(), cv_coefs.begin()+(j+count+((q-m+1)*0.5)+m-1)*2*ucount, cv_coefs.begin()+((j+count+((q-m+1)*0.5)+m-1)*2+(2*q-1))*ucount);

                std::vector<double> sf_coefs_parts;

                SplineInterpolator::quasiinterpolate(par_v_parts, coefs_parts,
                                                     tg_pnt, basis_v,
                                                     j+countj, sf_coefs_parts);
                sf_coefs.insert(sf_coefs.end(), sf_coefs_parts.begin(), sf_coefs_parts.end());
              }		    
            }//else

            dcount = q+(m-2);
          }//end if multivalue 
          else
          {
            if (dcount > 0)
            {
              std::vector<double> coefs_parts;
              std::vector<double> par_v_parts;

              par_v_parts.insert(par_v_parts.end(),par_v.begin()+(dcount+ti-(m-1))*2,par_v.begin()+(dcount+ti-(m-1))*2+(2*q-1));
              coefs_parts.insert(coefs_parts.end(), cv_coefs.begin()+(dcount+ti-(m-1))*2*ucount, cv_coefs.begin()+((dcount+ti-(m-1))*2+(2*q-1))*ucount);
              std::vector<double> sf_coefs_parts;

              SplineInterpolator::quasiinterpolate(par_v_parts, coefs_parts, 
                                                   tg_pnt, basis_v,
                                                   ti+dcount, sf_coefs_parts);
              sf_coefs.insert(sf_coefs.end(), sf_coefs_parts.begin(), sf_coefs_parts.end());
            }
            else
            {
              std::vector<double> coefs_parts;
              std::vector<double> par_v_parts;
              par_v_parts.insert(par_v_parts.end(),par_v.begin()+(ti+dcount)*2,par_v.begin()+(ti+dcount)*2+(2*q-1));
              coefs_parts.insert(coefs_parts.end(), cv_coefs.begin()+(ti+dcount)*2*ucount, cv_coefs.begin()+((ti+dcount)*2+(2*q-1))*ucount);

              std::vector<double> sf_coefs_parts;

              SplineInterpolator::quasiinterpolate(par_v_parts, coefs_parts,
                                                   tg_pnt, basis_v,
                                                   ti, sf_coefs_parts);
              sf_coefs.insert(sf_coefs.end(), sf_coefs_parts.begin(), sf_coefs_parts.end());
              count = count+1;
            }
          }
          ti++;
        }// end while
      }// end if q>1
    }//else


    ///////////////////////////////////////////////////////////////////

    //eval sf_coefs to generate sf_coefs_solution
    //boundary conditions
    int n, m, spxl, spyl, qpxl, qpyl, invxl, invyl, gmxl, gmyl;
    n = basis_u.numCoefs();
    m = basis_v.numCoefs();
    spxl = 2*p-1;
    spyl = 2*q-1;
    qpxl = 2*(n-p)+1;
    qpyl = 2*(m-q)+1;
    invxl = qpxl-n+1;
    invyl = qpyl-m+1;
    gmxl = spxl*invxl;
    gmyl = spyl*invyl;

    std::vector<std::vector<double> > gm;
    std::vector<double> tmp;
    for(int i = 0; i < gmyl ; i++)
    { 
      for(int  j=0 ; j < gmxl*perknot  ;j++)
      {
        tmp.push_back( 0.0 );
      }
      gm.push_back(tmp);
    }

    std::vector<double> sf_coefs_solution;
    for(int i = 0; i < n*m*perknot ; i++)
      sf_coefs_solution.push_back( 0.0 );

    int starti, startj;
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
          starti = b*(2*q-1);
          startj = a*(2*p-1);
          {if (invyl == 1)
            {if (a == 0)
              for (int i = starti;i<(starti+2*q-1);i++)
                for (int j = startj;j<(startj+p)*perknot;j++)
                  gm[i][j] = 1;
              else if (a == invxl-1)
                for (int i = starti;i<(starti+2*q-1);i++)
                  for (int j = (startj+p-1)*perknot; j< (startj+2*p-1)*perknot;j++)
                    gm[i][j] = 1;
              for (int i = starti;i<(starti+2*q-1);i++)
                for (int j = (startj+p-1)*perknot;j<(startj+p)*perknot;j++)
                  gm[i][j] = 1;}
              else if (invxl == 1) // 
              {if (b == 0)
                for (int i = starti; i< (starti+q);i++)
                  for (int j = startj;j<(startj+2*p-1)*perknot;j++)
                    gm[i][j] = 1;
                else if (b == invyl-1)
                  for (int i = (starti+q-1);i<(starti+2*q-1);i++)
                    for (int j = startj;j<(startj+2*p-1)*perknot;j++)
                      gm[i][j] = 1;
                for (int i = (starti+q-1);i<(starti+q);i++)
                  for (int j = startj;j<(startj+2*p-1)*perknot;j++)
                    gm[i][j] = 1;}
                else if (b==0 && a == 0) //top left corner
                  for (int i = starti;i<(starti+q);i++)
                    for (int j = startj;j<(startj+p)*perknot;j++)
                      gm[i][j]= 1;
                else if (b == 0 && a == invxl-1) //top right corner
                  for (int i = starti;i<(starti+q);i++)
                    for (int j = (startj+p-1)*perknot;j<(startj+2*p-1)*perknot;j++)
                      gm[i][j]= 1;
                else if (b == invyl-1 && a == 0) //bottom left corner
                  for (int i = (starti+q-1);i<(starti+2*q-1);i++)
                    for (int j = startj;j<(startj+p)*perknot;j++)
                      gm[i][j]= 1;
                else if (b == invyl-1 && a == invxl-1) //bottom right corner
                  for (int i = (starti+q-1);i<(starti+2*q-1);i++)
                    for (int j = (startj+p-1)*perknot;j<(startj+2*p-1)*perknot;j++)
                      gm[i][j]= 1;
                else if (b==0) // top edge
                  for (int i = starti;i<(starti+q);i++)
                    for (int j = (startj+p-1)*perknot;j<(startj+p)*perknot;j++)
                      gm[i][j]= 1;
                else if (a == invxl-1) // right edge
                  for (int i = (starti+q-1);i<(starti+q);i++)
                    for (int j = (startj+p-1)*perknot;j<(startj+2*p-1)*perknot;j++)
                      gm[i][j]= 1;
                else if (a == 0) // left edge
                  for (int i = (starti+q-1);i<(starti+q);i++)
                    for (int j = startj;j<(startj+p)*perknot;j++)
                      gm[i][j]= 1;
                else if (b == invyl-1) // bottom edge
                  for (int i = (starti+q-1);i<(starti+2*q-1);i++)
                    for (int j = (startj+p-1)*perknot;j<(startj+p)*perknot;j++)
                      gm[i][j]= 1;}// interior elements
                for (int i = (starti+q-1);i<(starti+q);i++)
                  for (int j = (startj+p-1)*perknot;j<(startj+p)*perknot;j++)
                    gm[i][j]= 1;

        }
      }

      int count = -1;
      for (int i = 0;i<gmyl;i++)
        for (int j = 0;j<gmxl*perknot;j++)
          if (gm[i][j] == 1)
          {count = count+1;
            sf_coefs_solution[count] = sf_coefs[i*(gmxl*perknot)+j];}
    }

    volinput_coefs.insert(volinput_coefs.end(), sf_coefs_solution.begin(), sf_coefs_solution.end());

  }

  // Interpolate surfaces to create volume
  std::vector<double> vol_coefs;
  SplineInterpolator::interpolate(par_w, volinput_coefs, tg_pnt, basis_w, vol_coefs);

  // Make volume
  return new Go::SplineVolume(basis_u, basis_v, basis_w, vol_coefs.begin(),
                              dimension, rational);
}


/*!
  \brief Local projection method (Variation Diminishing Spline Approximation).
  \param[in] basis_u Basis values in the first parameter direction
  \param[in] basis_v Basis values in the second parameter direction
  \param[in] basis_w Basis values in the third parameter direction
  \param[in] par_u Grevielle sites in the first parameter direction
  \param[in] par_v Grevielle sites in the second parameter direction
  \param[in] par_w Grevielle sites in the third parameter direction
  \param[in] points Secondary solution field evaluated at Greville points
  \param[in] dimension Dimension of the secondary solution field
  \param[in] rational Value marks NURBS geometry
  \param[in] weights NURBS weights for the projective control points
  \return Spline volume object representing the projected field

  \note VariationDiminishingSplineApproximation only for function values! 
*/

static Go::SplineVolume*
VariationDiminishingSplineApproximation(const Go::BsplineBasis& basis_u,
					const Go::BsplineBasis& basis_v,
					const Go::BsplineBasis& basis_w,
					const RealArray& par_u,
					const RealArray& par_v,
					const RealArray& par_w,
					const RealArray& points,
					int dimension, bool rational,
					const RealArray& weights)
{
  // Check input
  ASSERT(par_u.size()*par_v.size()*par_w.size() == points.size()/dimension);
  ASSERT(basis_u.numCoefs() == (int)par_u.size());
  ASSERT(basis_v.numCoefs() == (int)par_v.size());
  ASSERT(basis_w.numCoefs() == (int)par_w.size());

  std::vector<double> local_coefs;
  if (rational)
  {
    ASSERT(weights.size() == points.size()/dimension);
    size_t sizepoints = par_u.size()*par_v.size()*par_w.size();
    size_t countpoints = 0;
    for (size_t i = 0; i < sizepoints; i++)
    {
      for (int j = 0; j < dimension; j++, countpoints++)
        local_coefs.push_back(points[countpoints]*weights[i]);
      local_coefs.push_back(weights[i]);
    }
  }
  else
    local_coefs = points;

  // Make volume
  return new Go::SplineVolume(basis_u, basis_v, basis_w, local_coefs.begin(),
                              dimension, rational);
}
