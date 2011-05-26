// $Id$
//==============================================================================
//!
//! \file Legendre.C
//!
//! \date Mar 19 2009
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Various utility methods for Spectral elements.
//!
//==============================================================================

#include "Legendre.h"
#include "DenseMatrix.h"


bool Legendre::GL (RealArray& weights, RealArray& points, int n)
{
  points.resize(n);
  weights.resize(n);
  if (n < 2) return false;

  DenseMatrix A(n,n);
  A(1,2) = real(1);

  int i;
  if (n > 2)
  {
    for (i = 2; i < n; i++)
    {
      A(i,i-1) = real(i-1)/real(2*i-1);
      A(i,i+1) = real(i)  /real(2*i-1);
    }
    A(n,n-1) = real(n-1)/real(2*n-1);
  }

  Vector eig_complex, evalpoints;
  if (!A.solveEigNon(evalpoints,eig_complex)) return false;
  std::sort(evalpoints.begin(),evalpoints.end());

  points = evalpoints;

  real L;
  for (i = 0; i < n; i++)
    if (!LegendreDerEval(n,evalpoints[i],L))
      return false;
    else
      weights[i] = real(2)/((real(1)-evalpoints[i]*evalpoints[i])*L*L);

  return true;
}



bool Legendre::GLL (Vector& weights, Vector& points, int n)
{
  weights.resize(n);
  points.resize(n);

  points(1) = -real(1);
  points(n) =  real(1);

  int i;
  real L, Ld;

  if (n < 3)
  {
    for (i = 1; i <= n; i++)
      if (!LegendreEval(n-1,points(i),L))
	return false;
      else
	weights(i) = real(2)/((n-1)*n*L*L);

    return true;
  }

  RealArray xw, xg;
  if (!GL(xw,xg,n-1)) return false;

  const real tol = 1.0e-8;

  for (i = 2; i < n; i++)
  {
    points(i) = (xg[i-2] + xg[i-1])/real(2);
    real ptemp, res = 1;
    while (res > tol)
    {
      ptemp = points(i);
      if (!LegendreEval(n-1,ptemp,L)) return false;
      if (!LegendreDerEval(n-1,ptemp,Ld)) return false;
      points(i) = ptemp + ((real(1)-ptemp*ptemp)*Ld)/((n-1)*n*L);
      res = fabs(points(i) - ptemp);
    }
  }

  for (i = 1; i <= n; i++)
    if (!LegendreEval(n-1,points(i),L))
      return false;
    else
      weights(i) = real(2)/((n-1)*n*L*L);

  return true;
}


bool Legendre::LegendreEval (int n, real x, real& retval)
{
  if (n < 1)
  {
    std::cerr <<" *** Legendre::LegendreEval: Polynomial index "<< n
	      <<" out of range: n < 1"<< std::endl;
    return false;
  }

  if (x < -real(1) || x > real(1))
  {
    std::cerr <<" *** Legendre::LegendreEval: Evaluation point "<< x
	      <<" out of range: [-1,1]"<< std::endl;
    return false;
  }

  Vector val(n+1);
  val(1) = real(1);
  val(2) = x;
  for (int i = 2; i <= n; i++)
    val(i+1) = real(2*i-1)/real(i)*x*val(i) - real(i-1)/real(i)*val(i-1);

  retval = val(n+1);
  return true;
}


bool Legendre::LegendreDerEval (int n, real x, real& retval)
{
  if (n < 1)
  {
    std::cerr <<" *** Legendre::LegendreDerEval: Polynomial index "<< n
	      <<" out of range: n < 1"<< std::endl;
    return false;
  }

  if (x < -real(1) || x > real(1))
  {
    std::cerr <<" *** Legendre::LegendreDerEval: Evaluation point "<< x
	      <<" out of range: [-1,1]"<< std::endl;
    return false;
  }

  if (x == real(1) || x == -real(1))
  {
    retval = pow(x,n-1)*real(n)*real(n+1)/real(2);
    return true;
  }

  Vector val(n+1);
  val(1) = real(1);
  val(2) = x;
  for (int i = 2; i <= n; i++)
    val(i+1) = real(2*i-1)/real(i)*x*val(i) - real(i-1)/real(i)*val(i-1);

  retval = real(n)/(real(1)-x*x)*val(n) - real(n)*x/(real(1)-x*x)*val(n+1);
  return true;
}


bool Legendre::basisDerivatives (int n, Matrix& D)
{
  D.resize(n,n);

  Vector w,p;
  if (!GLL(w,p,n)) return false;

  real l1,l2;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      if (i == j)
	D(i,j) = real(0);
      else
      {
	if (!LegendreEval(n-1,p(i),l1)) return false;
	if (!LegendreEval(n-1,p(j),l2)) return false;
	D(i,j) = l1/(l2*(p(i)-p(j)));
      }

  D(1,1) = -real(n)*real(n-1)/real(4);
  D(n,n) = -D(1,1);
  return true;
}
