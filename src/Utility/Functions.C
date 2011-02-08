// $Id: Functions.C,v 1.7 2011-01-05 12:45:14 kmo Exp $
//==============================================================================
//!
//! \file Functions.C
//!
//! \date Jun 1 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Specific function implementations.
//!
//==============================================================================

#include "Functions.h"
#include "Vec3.h"
#include <string.h>
#include <stdlib.h>


PressureField::PressureField (real p, int dir) : pdir(dir)
{
  pressure = new ConstFunc(p);
}


real ConstTimeFunc::evaluate (const Vec3& x) const
{
  const Vec4* X = dynamic_cast<const Vec4*>(&x);
  return (*tfunc)(X ? X->t : real(0));
}


real SpaceTimeFunc::evaluate (const Vec3& x) const
{
  const Vec4* X = dynamic_cast<const Vec4*>(&x);
  return (*sfunc)(x) * (*tfunc)(X ? X->t : real(0));
}


real LinearTinitFunc::evaluate (const Vec3& x) const
{
  const Vec4* X = dynamic_cast<const Vec4*>(&x);
  if (X && X->t < Tinit)
    return value*X->t/Tinit;
  else
    return value;
}


real LinearXFunc::evaluate (const Vec3& X) const
{
  return a*X.x + b;
}


real LinearYFunc::evaluate (const Vec3& X) const
{
  return a*X.y + b;
}


real LinearZFunc::evaluate (const Vec3& X) const
{
  return a*X.z + b;
}


real LinearRotZFunc::evaluate (const Vec3& _X) const
{
  // Always return zero if the argument has no time component
  const Vec4* X = dynamic_cast<const Vec4*>(&_X);
  if (!X) return real(0);

  real x = X->x - x0;
  real y = X->y - y0;
  real c = cos(A*X->t);
  real s = sin(A*X->t);
  return rX ? x*c-y*s-x : x*s+y*c-y;
}


real StepXFunc::evaluate (const Vec3& X) const
{
  return X.x < x0 || X.x > x1 ? 0.0 : fv;
}



real StepXYFunc::evaluate (const Vec3& X) const
{
  return X.x < x0 || X.x > x1 || X.y < y0 || X.y > y1 ? 0.0 : fv;
}


const RealFunc* utl::parseRealFunc (char* cline, real A)
{
  // Check for spatial variation
  int linear = 0;
  if (cline)
    if (strcmp(cline,"X") == 0)
      linear = 1;
    else if (strcmp(cline,"Y") == 0)
      linear = 2;
    else if (strcmp(cline,"Z") == 0)
      linear = 3;
    else if (strcmp(cline,"XrotZ") == 0)
      linear = 4;
    else if (strcmp(cline,"YrotZ") == 0)
      linear = 5;
    else if (strcmp(cline,"Tinit") == 0)
      linear = 6;
    else if (strcmp(cline,"StepX") == 0)
      linear = 7;
    else if (strcmp(cline,"StepXY") == 0)
      linear = 8;

  real C = A;
  const RealFunc* f = 0;
  if (linear && (cline = strtok(NULL," ")))
  {
    C = real(1);
    std::cout <<"("<< A <<"*";
    if (linear < 4)
      std::cout << char('W' + linear) <<" + "<< cline <<")";
    else if (linear < 6)
      std::cout << char('W' + linear-3) <<"RotZ("<< cline <<"))";
    switch (linear) {
    case 1:
      f = new LinearXFunc(A,atof(cline));
      break;
    case 2:
      f = new LinearYFunc(A,atof(cline));
      break;
    case 3:
      f = new LinearZFunc(A,atof(cline));
      break;
    case 4:
      f = new LinearRotZFunc(true,A,atof(cline),atof(strtok(NULL," ")));
      break;
    case 5:
      f = new LinearRotZFunc(false,A,atof(cline),atof(strtok(NULL," ")));
      break;
    case 6:
      std::cout <<"RampT("<< cline <<"))";
      f = new LinearTinitFunc(A,atof(cline));
      break;
    case 7:
      std::cout <<"StepX("<< cline <<"))";
      f = new StepXFunc(A,atof(cline),atof(strtok(NULL," ")));
      break;
    case 8:
      std::cout <<"StepXY("<< cline <<"))";
      f = new StepXYFunc(A,atof(cline),atof(strtok(NULL," ")));
      break;
    }
    cline = strtok(NULL," ");
  }
  else
    std::cout << C;

  // Check for time variation
  if (!cline) return f;

  const ScalarFunc* s = 0;
  double freq = atof(cline);
  if (cline = strtok(NULL," "))
  {
    std::cout <<" * sin("<< freq <<"*t + "<< cline <<")";
    s = new SineFunc(C,freq,atof(cline));
  }
  else
  {
    std::cout <<" * "<< freq <<"*t";
    s = new LinearFunc(C*freq);
  }

  if (f)
    return new SpaceTimeFunc(f,s);
  else
    return new ConstTimeFunc(s);
}
