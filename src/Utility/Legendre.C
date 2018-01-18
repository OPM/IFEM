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
  A(1,2) = Real(1);

  int i;
  if (n > 2)
  {
    for (i = 2; i < n; i++)
    {
      A(i,i-1) = Real(i-1)/Real(2*i-1);
      A(i,i+1) = Real(i)  /Real(2*i-1);
    }
    A(n,n-1) = Real(n-1)/Real(2*n-1);
  }

  RealArray eig_complex, evalpoints;
  if (!A.solveEigNon(evalpoints,eig_complex)) return false;
  std::sort(evalpoints.begin(),evalpoints.end());

  points = evalpoints;

  Real L;
  for (i = 0; i < n; i++)
    if (!Legendre::derEval(n,evalpoints[i],L))
      return false;
    else
      weights[i] = Real(2)/((Real(1)-evalpoints[i]*evalpoints[i])*L*L);

  return true;
}



bool Legendre::GLL (RealArray& weights, RealArray& points, int n)
{
  weights.resize(n);
  points.resize(n);

  points.front() = -Real(1);
  points.back()  =  Real(1);

  int i;
  Real L, Ld;

  if (n < 3)
  {
    for (i = 0; i < n; i++)
      if (!Legendre::eval(n-1,points[i],L))
        return false;
      else
        weights[i] = Real(2)/((n-1)*n*L*L);

    return true;
  }

  RealArray xw, xg;
  if (!GL(xw,xg,n-1)) return false;

  const Real tol = 1.0e-8;

  for (i = 1; i+1 < n; i++)
  {
    points[i] = (xg[i-1] + xg[i])/Real(2);
    Real ptemp, res = 1;
    while (res > tol)
    {
      ptemp = points[i];
      if (!Legendre::eval(n-1,ptemp,L)) return false;
      if (!Legendre::derEval(n-1,ptemp,Ld)) return false;
      points[i] = ptemp + ((Real(1)-ptemp*ptemp)*Ld)/((n-1)*n*L);
      res = fabs(points[i] - ptemp);
    }
  }

  for (i = 0; i < n; i++)
    if (!Legendre::eval(n-1,points[i],L))
      return false;
    else
      weights[i] = Real(2)/((n-1)*n*L*L);

  return true;
}


bool Legendre::eval (int n, Real x, Real& retval)
{
  if (n < 1)
  {
    std::cerr <<" *** Legendre::eval: Polynomial index "<< n
              <<" out of range: n < 1"<< std::endl;
    return false;
  }

  if (x < -Real(1) || x > Real(1))
  {
    std::cerr <<" *** Legendre::eval: Evaluation point "<< x
              <<" out of range: [-1,1]"<< std::endl;
    return false;
  }

  RealArray val(n+1);
  val[0] = Real(1);
  val[1] = x;
  for (int i = 2; i <= n; i++)
    val[i] = Real(2*i-1)/Real(i)*x*val[i-1] - Real(i-1)/Real(i)*val[i-2];

  retval = val.back();
  return true;
}


bool Legendre::derEval (int n, Real x, Real& retval)
{
  if (n < 1)
  {
    std::cerr <<" *** Legendre::derEval: Polynomial index "<< n
              <<" out of range: n < 1"<< std::endl;
    return false;
  }

  if (x < -Real(1) || x > Real(1))
  {
    std::cerr <<" *** Legendre::derEval: Evaluation point "<< x
              <<" out of range: [-1,1]"<< std::endl;
    return false;
  }

  if (x == Real(1) || x == -Real(1))
  {
    retval = pow(x,Real(n-1))*Real(n)*Real(n+1)/Real(2);
    return true;
  }

  RealArray val(n+1);
  val[0] = Real(1);
  val[1] = x;
  for (int i = 2; i <= n; i++)
    val[i] = Real(2*i-1)/Real(i)*x*val[i-1] - Real(i-1)/Real(i)*val[i-2];

  retval = Real(n)/(Real(1)-x*x)*val[n-1] - Real(n)*x/(Real(1)-x*x)*val[n];
  return true;
}


bool Legendre::basisDerivatives (int n, Matrix& D)
{
  D.resize(n,n);

  Vector w,p;
  if (!GLL(w,p,n)) return false;

  Real l1,l2;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      if (i == j)
        D(i,j) = Real(0);
      else
      {
        if (!Legendre::eval(n-1,p(i),l1)) return false;
        if (!Legendre::eval(n-1,p(j),l2)) return false;
        D(i,j) = l1/(l2*(p(i)-p(j)));
      }

  D(1,1) = -Real(n)*Real(n-1)/Real(4);
  D(n,n) = -D(1,1);
  return true;
}
