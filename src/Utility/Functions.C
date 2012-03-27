// $Id$
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
#include "Vec3Oper.h"
#include <cstring>
#include <fstream>
#include <algorithm>


PressureField::PressureField (real p, int dir) : pdir(dir)
{
  pressure = new ConstFunc(p);
}


real RampFunc::evaluate (const real& x) const
{
  return x < xmax ? fval*x/xmax : fval;
}


real DiracFunc::evaluate (const real& x) const
{
  return fabs(x-xmax) < 1.0e-4 ? amp : real(0);
}


real StepFunc::evaluate (const real& x) const
{
  return x >= xmax ? amp : real(0);
}


real SineFunc::evaluate (const real& x) const
{
  return scale*sin(freq*x+phase);
}


real ConstTimeFunc::evaluate (const Vec3& X) const
{
  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  return (*tfunc)(Xt ? Xt->t : real(0));
}


real SpaceTimeFunc::evaluate (const Vec3& X) const
{
  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  return (*sfunc)(X) * (*tfunc)(Xt ? Xt->t : real(0));
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


real QuadraticXFunc::evaluate (const Vec3& X) const
{
  real val = (a-b)/real(2);
  return max*(a-X.x)*(X.x-b)/(val*val);
}


real QuadraticYFunc::evaluate (const Vec3& X) const
{
  real val = (a-b)/real(2);
  return max*(a-X.y)*(X.y-b)/(val*val);
}


real QuadraticZFunc::evaluate (const Vec3& X) const
{
  real val = (a-b)/real(2);
  return max*(a-X.z)*(X.z-b)/(val*val);
}


real LinearRotZFunc::evaluate (const Vec3& X) const
{
  // Always return zero if the argument has no time component
  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  if (!Xt) return real(0);

  real x = X.x - x0;
  real y = X.y - y0;
  real c = cos(A*Xt->t);
  real s = sin(A*Xt->t);
  return rX ? x*c-y*s-x : x*s+y*c-y;
}


real StepXFunc::evaluate (const Vec3& X) const
{
  return X.x < x0 || X.x > x1 ? real(0) : fv;
}


real StepXYFunc::evaluate (const Vec3& X) const
{
  return X.x < x0 || X.x > x1 || X.y < y0 || X.y > y1 ? real(0) : fv;
}


Interpolate1D::Interpolate1D (const char* file, int dir_) : dir(dir_)
{
  std::ifstream is(file);
  while (is.good()) {
    real x, v;
    is >> x >> v;
    grid.push_back(x);
    values.push_back(v);
  }
}


real Interpolate1D::evaluate (const Vec3& X) const
{
  real x = X[dir];
  std::vector<real>::const_iterator xb =
    std::find_if(grid.begin(),grid.end()-1,
		 std::bind2nd(std::greater<real>(),x));

  size_t pos = xb-grid.begin();
  real x1 = *(xb-1);
  real x2 = *xb;
  real val1 = values[pos-1];
  real val2 = values[pos];

  double delta = x2 - x1;
  if (fabs(delta) < 1.0e-8)
    return val1;

  return (val1*(x-x1) + val2*(x2-x)) / delta;
}


/*!
  The functions are assumed on the general form
  \f[ f({\bf X},t) = A * g({\bf X}) * h(t) \f]

  The character string \a cline is assumed to contain first the definition
  of the spatial function \a g( \b X ) and then the time function \a h(t).
  Either of the two components may be omitted, for creating a space-function
  constant in time, or a time function constant in space.
*/

const RealFunc* utl::parseRealFunc (char* cline, real A)
{
  // Check for spatial variation
  int linear    = 0;
  int quadratic = 0;
  if (!cline)
    linear = -1;
  else if (strcmp(cline,"X") == 0)
    linear = 1;
  else if (strcmp(cline,"Y") == 0)
    linear = 2;
  else if (strcmp(cline,"Z") == 0)
    linear = 3;
  else if (strcmp(cline,"XrotZ") == 0)
    linear = 4;
  else if (strcmp(cline,"YrotZ") == 0)
    linear = 5;
  else if (strcmp(cline,"StepX") == 0)
    linear = 6;
  else if (strcmp(cline,"StepXY") == 0)
    linear = 7;
  else if (strcmp(cline,"Interpolate1D") == 0)
    linear = 8;
  else if (strcmp(cline,"quadX") == 0)
    quadratic = 1;
  else if (strcmp(cline,"quadY") == 0)
    quadratic = 2;
  else if (strcmp(cline,"quadZ") == 0)
    quadratic = 3;

  real C = A;
  const RealFunc* f = 0;
  if (linear > 0 && (cline = strtok(NULL," ")))
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
      cline = strtok(NULL," ");
      break;
    case 2:
      f = new LinearYFunc(A,atof(cline));
      cline = strtok(NULL," ");
      break;
    case 3:
      f = new LinearZFunc(A,atof(cline));
      cline = strtok(NULL," ");
      break;
    case 4:
      f = new LinearRotZFunc(true,A,atof(cline),atof(strtok(NULL," ")));
      cline = strtok(NULL," ");
      break;
    case 5:
      f = new LinearRotZFunc(false,A,atof(cline),atof(strtok(NULL," ")));
      cline = strtok(NULL," ");
      break;
    case 6:
      {
	double x0 = atof(cline);
	double x1 = atof(strtok(NULL," "));
	std::cout <<"StepX("<< x0 <<","<< x1 <<"))";
	f = new StepXFunc(A,x0,x1);
      }
      cline = strtok(NULL," ");
      break;
    case 7:
      {
	double x0 = atof(cline);
	double y0 = atof(strtok(NULL," "));
	cline = strtok(NULL," ");
	if (cline && cline[0] == 't')
	{
	  double x1 = atof(strtok(NULL," "));
	  double y1 = atof(strtok(NULL," "));
	  std::cout <<"StepXY(["<< x0<<","<<x1 <<"]x["<< y0<<","<<y1 <<"]))";
	  f = new StepXYFunc(A,x1,y1,x0,y0);
	  cline = strtok(NULL," ");
	}
	else
	{
	  std::cout <<"StepXY([-inf,"<< x0 <<"]x[-inf,"<< y0 <<"]))";
	  f = new StepXYFunc(A,x0,y0);
	}
      }
      break;
    case 8:
      {
        int dir = atoi(strtok(NULL, " "));
        std::cout <<"Interpolate1D("<< cline <<","<< (char)('X'+dir) <<")";
        f = new Interpolate1D(cline,dir);
        cline = strtok(NULL," ");
      }
      break;
    }
  }
  else if (quadratic > 0 && (cline = strtok(NULL," ")))
  {
    C = real(1);
    real a = atof(cline);
    real b = atof(strtok(NULL," "));
    real val = (a-b)*(a-b)/real(4);
    char var = 'W' + quadratic;
    std::cout << A/val <<" * ("<< a <<"-"<< var <<")*("<< b <<"-"<< var <<")";
    switch (quadratic) {
    case 1:
      f = new QuadraticXFunc(A,a,b);
      break;
    case 2:
      f = new QuadraticYFunc(A,a,b);
      break;
    case 3:
      f = new QuadraticZFunc(A,a,b);
      break;
    }
    cline = strtok(NULL," ");
  }
  else // constant in space
  {
    std::cout << C;
    if (linear < 0)
      f = new ConstFunc(C);
  }

  // Check for time variation
  if (!cline) return f; // constant in time

  const ScalarFunc* s = 0;
  if (strncmp(cline,"Ramp",4) == 0 || strcmp(cline,"Tinit") == 0)
  {
    real xmax = atof(strtok(NULL," "));
    std::cout <<" * Ramp(t,"<< xmax <<")";
    s = new RampFunc(C,xmax);
  }
  else if (strncmp(cline,"Dirac",5) == 0)
  {
    real xmax = atof(strtok(NULL," "));
    std::cout <<" * Dirac(t,"<< xmax <<")";
    s = new DiracFunc(C,xmax);
  }
  else if (strncmp(cline,"Step",4) == 0)
  {
    real xmax = atof(strtok(NULL," "));
    std::cout <<" * Step(t,"<< xmax <<")";
    s = new StepFunc(C,xmax);
  }
  else if (strcmp(cline,"sin") == 0)
  {
    real freq = atof(strtok(NULL," "));
    if ((cline = strtok(NULL," ")))
    {
      real phase = atof(cline);
      std::cout <<" * sin("<< freq <<"*t + "<< phase <<")";
      s = new SineFunc(C,freq,phase);
    }
    else
    {
      std::cout <<" * sin("<< freq <<"*t)";
      s = new SineFunc(C,freq);
    }
  }
  else // linear in time
  {
    real scale = atof(cline);
    std::cout <<" * "<< scale <<"*t";
    s = new LinearFunc(C*scale);
  }

  if (f)
    return new SpaceTimeFunc(f,s);
  else
    return new ConstTimeFunc(s);
}


RealFunc* utl::parseRealFunc (const std::string& func, const std::string& type)
{
  if (func.empty()) return NULL;

  std::cout <<": ";
  real p = real(0);
  if (type == "constant")
    p = atof(func.c_str());
  else if (type == "linear")
  {
    p = atof(func.c_str());
    std::cout << p <<"*t";
    return new ConstTimeFunc(new LinearFunc(p));
  }
  else if (type.substr(0,10) == "expression")
  {
    std::cout << func;
    return new EvalFunction(func.c_str());
  }
  else
  {
    std::string tmp(func);
    p = atof(strtok(const_cast<char*>(tmp.c_str())," "));
    char* sfun = strtok(NULL," ");
    if (sfun) return const_cast<RealFunc*>(parseRealFunc(sfun,p));
  }

  std::cout << p;
  return new ConstFunc(p);
}


VecFunc* utl::parseVecFunc (const std::string& func, const std::string& type)
{
  if (func.empty()) return NULL;

  if (type == "constant")
  {
    Vec3 v;
    std::string tmp(func);
    char* s = strtok(const_cast<char*>(tmp.c_str())," ");
    for (int i = 0; i < 3 && s; i++, s = strtok(NULL," "))
      v[i] = atof(s);
    std::cout <<": "<< v;
    return new ConstVecFunc(v);
  }
  else if (type.substr(0,10) == "expression")
  {
    std::cout <<": "<< func;
    return new VecFuncExpr(func.c_str());
  }

  return NULL;
}


TractionFunc* utl::parseTracFunc (const std::string& func,
				  const std::string& type, int dir)
{
  if (func.empty()) return NULL;

  std::cout <<": ";
  real p = real(0);
  const RealFunc* f = 0;
  if (type == "constant")
  {
    p = atof(func.c_str());
    std::cout << p;
  }
  else if (type == "linear")
  {
    p = atof(func.c_str());
    f = new ConstTimeFunc(new LinearFunc(p));
    std::cout << p <<"*t";
  }
  else if (type.substr(0,10) == "expression")
  {
    std::cout << func;
    f = new EvalFunction(func.c_str());
  }
  else
  {
    std::string tmp(func);
    p = atof(strtok(const_cast<char*>(tmp.c_str())," "));
    char* sfun = strtok(NULL," ");
    if (sfun)
      f = parseRealFunc(sfun,p);
    else
      std::cout << p;
  }

  return f ? new PressureField(f,dir) : new PressureField(p,dir);
}
