//==============================================================================
//!
//! \file SplineInterpolator.C
//!
//! \date Jan 2012
//!
//! \author Anette Stahl
//!
//! \brief Implementation of interpolation/projection schemes for B-splines
//!
//==============================================================================


#include <vector>
#include "GoTools/geometry/BsplineBasis.h"
#include "ASMstruct.h"
#include "DenseMatrix.h"
#include "SplineInterpolator.h"


//Global interpolation
void SplineInterpolator::interpolate(const std::vector<double>& params,
                                     const std::vector<double>& points,
                                     const std::vector<double>& tangent_points,
                                     const Go::BsplineBasis& basis,
                                     std::vector<double>& coefs)
{
  int num_points = (int)params.size();
  int dimension = (int)points.size() / num_points;
  int num_coefs =(int)basis.numCoefs();
  int order = basis.order();
  DEBUG_ERROR_IF(num_coefs < order,
      "Insufficient number of points.");

  coefs.resize(dimension*num_coefs);

  int i, j;

  DenseMatrix A(num_coefs, num_coefs);
  Matrix b(num_coefs, dimension);

  // setting up interpolation matrix A
  int ti = 0; // index to first unused element of tangent_points
  std::vector<double> tmp(2*order);
  for (i = 0; i < num_points; ++i) {
    bool der = false;
    double par = params[i];
    int ki = basis.knotIntervalFuzzy(par); // knot-interval of param.
    basis.computeBasisValues(params[i], &tmp[0], 1);

    for (j = 0; j < order; ++j)
      if ((ki-order+1+j>=0) && (ki-order+1+j<=num_coefs)) {
        A(i+ti+1,ki-order+1+j+1) = tmp[2*j];
        if (der)
          A(i+ti+1,ki-order+1+j+1) = tmp[2*j+1];
      }
    if (der)
      ++ti;
  }

  // generating right-hand side
  for (i = 0; i < num_coefs; ++i) {
    for (int k=0;k<dimension;++k)
      b(i+1,k+1) = points[i*dimension+k];
  }

  // Now we are ready to solve Ac = b.  b will be overwritten by solution
  A.solve(b);
  // copy results
  int k=0;
  for (int i=0;i<num_coefs;++i)
    for (int j=0;j<dimension;++j)
      coefs[k++] = b(i+1,j+1);
}


//Local Approximation - Quasi-Interpolation
void SplineInterpolator::quasiinterpolate(const std::vector<double>& params,
                                          const std::vector<double>& points,
                                          const std::vector<double>& tangent_points,
                                          const Go::BsplineBasis& basis,
                                          int index, std::vector<double>& coefs)
{
  int num_points = (int)params.size();
  int dimension = (int)points.size() / num_points;
  int num_coefs =(int)params.size();// !!!basis_.numCoefs();
  int order = basis.order();

  DEBUG_ERROR_IF(num_coefs < order,
      "Insufficient number of points.");

  coefs.resize(dimension*num_coefs);
  int i, j;
  DenseMatrix A(num_coefs, num_coefs);
  Matrix b(num_coefs, dimension);

  // setting up interpolation matrix A
  std::vector<double> tmp(2*order);
  for (i = 0; i < num_points; ++i) {
    bool der = false;
    double par = params[i];
    int ki = basis.knotIntervalFuzzy(par);

    basis.computeBasisValues(params[i], &tmp[0], 1);
    for (j = 0; j < order; ++j)
    {
      if((ki-order+1+j-index>=0) && (ki-order+1+j-index<num_coefs)) 
      {
        A(i+1,ki-order+1+j-index+1) = tmp[2*j];
        if (der)
          A(i+1+1,ki-order+1+j+1) = tmp[2*j+1];
      }
    }
  }

  // generating right-hand side
  for (i = 0; i < num_coefs; ++i) {
    for (int k=0;k<dimension;++k)
      b(i+1,k+1) = points[i*dimension+k];
  }

  A.solve(b);
  int k=0;
  for (int i=0;i<num_coefs;++i)
    for (int j=0;j<dimension;++j)
      coefs[k++] = b(i+1,j+1);
}
  
  
//Global Approxiamtion - Least-Square Method
void SplineInterpolator::leastsquare_approximation(const std::vector<double>& params,
                                                   const std::vector<double>& paramsweights,
                                                   const std::vector<double>& points,
                                                   const std::vector<double>& tangent_points,
                                                   const Go::BsplineBasis& basis,
                                                   std::vector<double>& coefs)
{
  int num_points = (int)params.size();
  int dimension = (int)points.size() / num_points;
  int num_coefs =(int)basis.numCoefs();
  int order = basis.order();

  DEBUG_ERROR_IF(num_coefs < order,
      "Insufficient number of points.");

  coefs.resize(dimension*num_coefs);
  //std::cout << " dimension" << dimension << std::endl;
  int i, j;

  Matrix A(num_points, num_coefs);
  Matrix b(num_points, dimension);

  // setting up interpolation matrix A
  int ti = 0; // index to first unused element of tangent_points
  std::vector<double> tmp(2*order);
  for (i = 0; i < num_points; ++i) {

    double par = params[i];
    int ki = basis.knotIntervalFuzzy(par); // knot-interval of param.
    basis.computeBasisValues(params[i], &tmp[0], 1);
    for (j = 0; j < order; ++j)
      if ((ki+1+j>=order) && (ki-order+1+j<num_coefs)) {
        A(i+ti+1,ki-order+1+j+1) = tmp[2*j];
      }
  }

  // generating right-hand side
  ti = 0;
  for (i = 0; i < num_points; ++i) {
    for (int k=0;k<dimension;++k)
      b(i+1,k+1) = points[i*dimension+k];
  }

  DenseMatrix Amass(num_points,num_points);
  Matrix bw(num_coefs,dimension);
  Matrix AwT(num_coefs,num_points);
  Matrix Aws(num_points,num_coefs);
  Matrix AwsT(num_coefs,num_points);

  // create Mass Matrix and weighted A Matrix
  for (i = 1; i <= num_points; ++i)
  {      
    for (j = 1; j <= num_coefs; ++j)
    {
      AwT(j,i) = paramsweights[i-1]*A(i,j);
      Aws(i,j) = sqrt(paramsweights[i-1])*A(i,j);
      AwsT(j,i) = Aws(i,j);
    }
  }

  Amass.getMat() = AwsT*Aws;

  bw = AwT*b;

  // Now we are ready to solve Ac = b.  b will be overwritten by solution
  Amass.solve(bw);
  int k=0;
  for (int i=0;i<num_coefs;++i)
    for (int j=0;j<dimension;++j)
      coefs[k++] = bw(i+1,j+1);
}
