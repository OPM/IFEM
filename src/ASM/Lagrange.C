// $Id$
//==============================================================================
//!
//! \file Lagrange.C
//!
//! \date Feb 10 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Evaluation of Lagrange basis functions.
//!
//==============================================================================

#include "Lagrange.h"


bool Lagrange::evalPol (int polnum, double xi, double& retval) const
{
  // Degree of polynomials = number of polynomials
  size_t degree = points.size();

  // Check if the polynomial number is feasible
  if (polnum < 0 || (size_t)polnum >= degree)
  {
    std::cerr <<" *** Lagrange::evalPol: Poynomial number out of range: "
	      << polnum <<", must be in the interval [0,"<< degree
	      <<">"<< std::endl;
    return false;
  }

  // Evaluate value of polynomial number polnum in point x
  retval = 1.0;
  for (size_t j = 0; j < degree; j++)
    if (j != (size_t)polnum)
      retval *= (xi - points[j]) / (points[polnum] - points[j]);

  return true;
}


bool Lagrange::evalDer (int polnum, double xi, double& retval) const
{
  // Check if evaluation point is inside the range of interpolation points
  if (xi < points.front() || xi > points.back())
  {
    std::cerr <<" *** Lagrange::evalDer: Evaluation point out of range: "
	      << xi <<", must be in the interval ["<< points.front()
	      <<","<< points.back() <<"]"<< std::endl;
    return false;
  }

  // Degree of polynomials = number of polynomials
  size_t degree = points.size();

  // Check if the polynomial number is feasible
  if (polnum < 0 || (size_t)polnum >= degree)
  {
    std::cerr <<" *** Lagrange::evalDer: Poynomial number out of range: "
	      << polnum <<", must be in the interval [0,"<< degree
	      <<">"<< std::endl;
    return false;
  }

  // Evaluate derivative of polynomial number polnum in point x
  retval = 0.0;
  for (size_t k = 0; k < degree; k++)
    if (k != (size_t)polnum)
    {
      double prod = 1.0;
      for (size_t j = 0; j < degree; j++)
	if (j != (size_t)polnum && j != k)
	  prod *= (xi - points[j]) / (points[polnum] - points[j]);
      retval += prod / (points[polnum] - points[k]);
    }

  return true;
}


bool Lagrange::computeBasis (RealArray& val,
			     int p1, double x1,
			     int p2, double x2,
			     int p3, double x3)
{
  // Define the Lagrangian interpolation points
  RealArray points1(p1), points2(p2), points3(p3);
  for (int i = 0; i < p1; i++) points1[i] = -1.0 + (i+i)/double(p1-1);
  for (int j = 0; j < p2; j++) points2[j] = -1.0 + (j+j)/double(p2-1);
  for (int k = 0; k < p3; k++) points3[k] = -1.0 + (k+k)/double(p3-1);

  return Lagrange::computeBasis(val,0,points1,x1,points2,x2,points3,x3);
}


bool Lagrange::computeBasis (RealArray& val,
			     Matrix& derval,
			     int p1, double x1,
			     int p2, double x2,
			     int p3, double x3)
{
  // Define the Lagrangian interpolation points
  RealArray points1(p1), points2(p2), points3(p3);
  for (int i = 0; i < p1; i++) points1[i] = -1.0 + (i+i)/double(p1-1);
  for (int j = 0; j < p2; j++) points2[j] = -1.0 + (j+j)/double(p2-1);
  for (int k = 0; k < p3; k++) points3[k] = -1.0 + (k+k)/double(p3-1);

  return Lagrange::computeBasis(val,&derval,points1,x1,points2,x2,points3,x3);
}


bool Lagrange::computeBasis (RealArray& val,
			     Matrix* derval,
			     const RealArray& points1, double x1,
			     const RealArray& points2, double x2,
			     const RealArray& points3, double x3)
{
  int p1 = points1.size();
  int p2 = points2.size();
  int p3 = points3.size();

  // Number of nodes per element
  int n1  = p1 > 0 ? p1 : 1;
  int n2  = p2 > 0 ? p2 : 1;
  int n3  = p3 > 0 ? p3 : 1;
  int nen = n1*n2*n3;

  // Set up the 1D Lagrangian polynomials
  Lagrange L1(points1);
  Lagrange L2(points2);
  Lagrange L3(points3);

  // Vectors of values for the one dimensional polynomials
  RealArray tempval1(n1,1.0), tempval2(n2,1.0), tempval3(n3,1.0);

  int i, j, k;

  // Evaluating values in each direction
  for (i = 0; i < p1; i++)
    if (!L1.evalPol(i,x1,tempval1[i]))
      return false;

  for (j = 0; j < p2; j++)
    if (!L2.evalPol(j,x2,tempval2[j]))
      return false;

  for (k = 0; k < p3; k++)
    if (!L3.evalPol(k,x3,tempval3[k]))
      return false;

  // Evaluate values of the 3-dimensional basis functions in the point x
  val.resize(nen);

  size_t count = 0;
  for (k = 0; k < n3; k++)
    for (j = 0; j < n2; j++)
      for (i = 0; i < n1; i++, count++)
	val[count] = tempval1[i]*tempval2[j]*tempval3[k];

  if (!derval) return true;

  Matrix& deriv = *derval;

  // Vectors of derivative values for the one dimensional polynomials
  RealArray tempder1(p1), tempder2(p2), tempder3(p3);

  // Evaluating values in each direction
  for (i = 0; i < p1; i++)
    if (!L1.evalDer(i,x1,tempder1[i]))
      return false;

  for (j = 0; j < p2; j++)
    if (!L2.evalDer(j,x2,tempder2[j]))
      return false;

  for (k = 0; k < p3; k++)
    if (!L3.evalDer(k,x3,tempder3[k]))
      return false;

  // Evaluate derivatives of the 3-dimensional basis functions in the point x
  derval->resize(nen, p3 > 0 ? 3 : (p2 > 0 ? 2 : 1));

  count = 1;
  for (k = 0; k < n3; k++)
    for (j = 0; j < n2; j++)
      for (i = 0; i < n1; i++, count++)
      {
	if (p1 > 0) deriv(count,1) = tempder1[i]*tempval2[j]*tempval3[k];
	if (p2 > 0) deriv(count,2) = tempval1[i]*tempder2[j]*tempval3[k];
	if (p3 > 0) deriv(count,3) = tempval1[i]*tempval2[j]*tempder3[k];
      }

  return true;
}
