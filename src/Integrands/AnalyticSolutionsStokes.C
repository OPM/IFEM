// $Id$
//==============================================================================
//!
//! \file AnalyticSolutionsStokes.C
//!
//! \date Des 10 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Analytic solutions for some Stokes and Navier-Stokes problems.
//!
//==============================================================================

#include "AnalyticSolutionsStokes.h"
#include "Vec3.h"

// Define pi
#ifndef M_PI
#define PI 3.141592653589793
#else
#define PI M_PI
#endif


/*!
  \class Poiseuille

  Channel flow.
*/

Poiseuille::Poiseuille(double P, double diam, double len, double visc)
  : Pin(P), D(diam), L(len), mu(visc)
{
  scalSol    = new Pressure(*this);
  vecSol     = new Velocity(*this);
  scalSecSol = new PressureGrad(*this);
  vecSecSol  = new VelocityGrad(*this);
}


// Pressure solution
double Poiseuille::Pressure::evaluate(const Vec3& x) const
{
  return params.Pin*(1.0-x[0]/params.L);
}


// Velocity solution
Vec3 Poiseuille::Velocity::evaluate(const Vec3& x) const
{
  double dPdx  = -params.Pin/params.L;
  double coeff = -dPdx/(2.0*params.mu);

  Vec3 u;
  u[0] = coeff*x[1]*(params.D-x[1]);

  return u;
}


// Pressure gradient
Vec3 Poiseuille::PressureGrad::evaluate(const Vec3& x) const
{
  Vec3 dPdX;
  dPdX[0] = -params.Pin/params.L;

  return dPdX;
}


// Velocity gradient
Tensor Poiseuille::VelocityGrad::evaluate(const Vec3& x) const
{
  Tensor dUdX(3);
  dUdX(1,2) = params.Pin/(2.0*params.L*params.mu)*(1.0-2.0*x[1]);

  return dUdX;
}


/*!
  \class TestSolution

  Analytical test solution using trigonometrical functions.
*/

TestSolution::TestSolution(double dens, double visc)
  : rho(dens), mu(visc)
{
  scalSol    = new Pressure();
  vecSol     = new Velocity();
  scalSecSol = new PressureGrad();
  vecSecSol  = new VelocityGrad();
}


// Pressure
double TestSolution::Pressure::evaluate(const Vec3& x) const
{
  const Vec4* X = dynamic_cast<const Vec4*>(&x);
  double t = X ? X->t : 0.0;

  return cos(PI*x[0])*sin(PI*x[1])*sin(t);
}


// Velocity
Vec3 TestSolution::Velocity::evaluate(const Vec3& x) const
{
  const Vec4* X = dynamic_cast<const Vec4*>(&x);
  double t = X ? X->t : 0.0;

  Vec3 u;
  u[0] = pow(sin(PI*x[0]),2.0)*sin(2.0*PI*x[1])*sin(t);
  u[1] = -sin(2.0*PI*x[0])*pow(sin(PI*x[1]),2.0)*sin(t);

  return u;
}


// Pressure gradient
Vec3 TestSolution::PressureGrad::evaluate(const Vec3& x) const
{
  const Vec4* X = dynamic_cast<const Vec4*>(&x);
  double t = X ? X->t : 0.0;

  Vec3 dPdX;
  dPdX[0] = -PI*sin(PI*x[0])*sin(PI*x[1])*sin(t);
  dPdX[1] =  PI*cos(PI*x[0])*cos(PI*x[1])*sin(t);

  return dPdX;
}


// Velocity gradient
Tensor TestSolution::VelocityGrad::evaluate(const Vec3& x) const
{
  const Vec4* X = dynamic_cast<const Vec4*>(&x);
  double t = X ? X->t : 0.0;

  Tensor dUdX(3);
  dUdX(1,1) =  2.0*PI*sin(PI*x[0])*cos(PI*x[0])*sin(2.0*PI*x[1])*sin(t);
  dUdX(1,2) =  2.0*PI*pow(sin(PI*x[0]),2.0)*cos(2.0*PI*x[1])*sin(t);
  dUdX(2,1) = -2.0*PI*cos(2.0*PI*x[0])*pow(sin(PI*x[1]),2.0)*sin(t);
  dUdX(2,2) = -2.0*PI*sin(2.0*PI*x[0])*sin(PI*x[1])*cos(PI*x[1])*sin(t);

  return dUdX;
}
