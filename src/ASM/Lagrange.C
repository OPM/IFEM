// $Id$
//==============================================================================
//!
//! \file Lagrange.C
//!
//! \date Feb 10 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Evaluation of %Lagrange basis functions.
//!
//==============================================================================

#include "Lagrange.h"


/*!
  \brief Class representing 1D %Lagrange polynomials.
*/

class Lagrange1D
{
  const RealArray& points; //!< Natural coordinates of the interpolation points

public:
  //! \brief Constructor initializing the reference to natural coordinates.
  //! \param[in] p Natural interpolation point coordinates in range [-1,1]
  explicit Lagrange1D(const RealArray& p) : points(p) {}

  //! \brief Evaluates a 1D %Lagrange polynomial or its derivatives
  //! \param[in] polnum Which polynomial of the basis to evaluate
  //! \param[in] xi Natural coordinate of the evaluation point
  //! \param[out] retval The computed polynomial value
  //! \param[in] deriv If \e true, evaluate the first derivatives
  bool eval(int polnum, double xi, double& retval, bool deriv = false) const
  {
    // Degree of polynomials = number of polynomials
    int degree = points.size();

    // Check if the polynomial number is feasible
    if (polnum < 0 || polnum >= degree)
    {
      std::cerr <<" *** Lagrange1D::eval: Poynomial number out of range: "
                << polnum <<", must be in the interval [0,"<< degree
                <<">"<< std::endl;
      return false;
    }

    // Check if evaluation point is inside the range of interpolation points
    if (deriv && (xi < points.front() || xi > points.back()))
    {
      std::cerr <<" *** Lagrange1D::eval: Evaluation point out of range: "
                << xi <<", must be in the interval ["<< points.front()
                <<","<< points.back() <<"]"<< std::endl;
      return false;
    }

    // Evaluate value of polynomial number polnum in point x
    retval = deriv ? 0.0 : 1.0;
    for (int j = 0; j < degree; j++)
      if (j != polnum)
      {
        if (deriv)
        {
          double prod = 1.0;
          for (int i = 0; i < degree; i++)
            if (i != polnum && i != j)
              prod *= (xi - points[i]) / (points[polnum] - points[i]);
          retval += prod / (points[polnum] - points[j]);
        }
        else
          retval *= (xi - points[j]) / (points[polnum] - points[j]);
      }

    return true;
  }
};


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

  return Lagrange::computeBasis(val,nullptr,points1,x1,points2,x2,points3,x3);
}


bool Lagrange::computeBasis (RealArray& val, Matrix& derval,
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


bool Lagrange::computeBasis (RealArray& val, Matrix* derval,
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

  // Set up the 1D Lagrangian polynomials
  Lagrange1D L1(points1);
  Lagrange1D L2(points2);
  Lagrange1D L3(points3);

  // Vectors of values for the one dimensional polynomials
  RealArray tempval1(n1,1.0), tempval2(n2,1.0), tempval3(n3,1.0);

  int i, j, k;

  // Evaluating values in each direction
  for (i = 0; i < p1; i++)
    if (!L1.eval(i,x1,tempval1[i]))
      return false;

  for (j = 0; j < p2; j++)
    if (!L2.eval(j,x2,tempval2[j]))
      return false;

  for (k = 0; k < p3; k++)
    if (!L3.eval(k,x3,tempval3[k]))
      return false;

  // Evaluate values of the 3-dimensional basis functions in the point x
  size_t nen = n1*n2*n3, count = 0;
  if (val.size() < nen)
    val.resize(nen);
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
    if (!L1.eval(i,x1,tempder1[i],true))
      return false;

  for (j = 0; j < p2; j++)
    if (!L2.eval(j,x2,tempder2[j],true))
      return false;

  for (k = 0; k < p3; k++)
    if (!L3.eval(k,x3,tempder3[k],true))
      return false;

  // Evaluate derivatives of the 3-dimensional basis functions in the point x
  size_t ndim = p3 > 0 ? 3 : (p2 > 0 ? 2 : 1);
  if (deriv.rows() < nen || deriv.cols() < ndim)
    deriv.resize(nen,ndim);

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
