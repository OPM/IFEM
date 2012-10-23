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
#include <sstream>
#include <algorithm>


PressureField::PressureField (Real p, int dir) : pdir(dir)
{
  pressure = new ConstFunc(p);
}


Real RampFunc::evaluate (const Real& x) const
{
  return x < xmax ? fval*x/xmax : fval;
}


Real DiracFunc::evaluate (const Real& x) const
{
  return fabs(x-xmax) < 1.0e-4 ? amp : Real(0);
}


Real StepFunc::evaluate (const Real& x) const
{
  return x >= xmax ? amp : Real(0);
}


Real SineFunc::evaluate (const Real& x) const
{
  return scale*sin(freq*x+phase);
}


Real ConstTimeFunc::evaluate (const Vec3& X) const
{
  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  return (*tfunc)(Xt ? Xt->t : Real(0));
}


Real SpaceTimeFunc::evaluate (const Vec3& X) const
{
  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  return (*sfunc)(X) * (*tfunc)(Xt ? Xt->t : Real(0));
}


Real LinearXFunc::evaluate (const Vec3& X) const
{
  return a*X.x + b;
}


Real LinearYFunc::evaluate (const Vec3& X) const
{
  return a*X.y + b;
}


Real LinearZFunc::evaluate (const Vec3& X) const
{
  return a*X.z + b;
}


Real QuadraticXFunc::evaluate (const Vec3& X) const
{
  Real val = (a-b)/Real(2);
  return max*(a-X.x)*(X.x-b)/(val*val);
}


Real QuadraticYFunc::evaluate (const Vec3& X) const
{
  Real val = (a-b)/Real(2);
  return max*(a-X.y)*(X.y-b)/(val*val);
}


Real QuadraticZFunc::evaluate (const Vec3& X) const
{
  Real val = (a-b)/Real(2);
  return max*(a-X.z)*(X.z-b)/(val*val);
}


Real LinearRotZFunc::evaluate (const Vec3& X) const
{
  // Always return zero if the argument has no time component
  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  if (!Xt) return Real(0);

  Real x = X.x - x0;
  Real y = X.y - y0;
  Real c = cos(A*Xt->t);
  Real s = sin(A*Xt->t);
  return rX ? x*c-y*s-x : x*s+y*c-y;
}


Real StepXFunc::evaluate (const Vec3& X) const
{
  return X.x < x0 || X.x > x1 ? Real(0) : fv;
}


Real StepXYFunc::evaluate (const Vec3& X) const
{
  return X.x < x0 || X.x > x1 || X.y < y0 || X.y > y1 ? Real(0) : fv;
}


Interpolate1D::Interpolate1D (const char* file, int dir_) : dir(dir_)
{
  std::ifstream is(file);
  while (is.good() && !is.eof()) {
    char temp[1024];
    is.getline(temp, 1024);
    if (is.eof())
      continue;
    std::stringstream str(temp);
    Real x, v;
    str >> x >> v;
    grid.push_back(x);
    values.push_back(v);
  }
}


Real Interpolate1D::evaluate (const Vec3& X) const
{
  Real x = X[dir];
  std::vector<Real>::const_iterator xb =
    std::find_if(grid.begin(),grid.end()-1,
		 std::bind2nd(std::greater<Real>(),x));

  size_t pos = xb-grid.begin();
  Real x1 = *(xb-1);
  Real x2 = *xb;
  Real val1 = values[pos-1];
  Real val2 = values[pos];

  double delta = x2 - x1;
  double alpha = (x2-x)/delta;

  return (val1*alpha + val2*(1-alpha));
}


/*!
  The functions are assumed on the general form
  \f[ f({\bf X},t) = A * g({\bf X}) * h(t) \f]

  The character string \a cline is assumed to contain first the definition
  of the spatial function \a g( \b X ) and then the time function \a h(t).
  Either of the two components may be omitted, for creating a space-function
  constant in time, or a time function constant in space.
*/

const RealFunc* utl::parseRealFunc (char* cline, Real A)
{
  // Check for spatial variation
  int linear    = 0;
  int quadratic = 0;
  if (!cline)
    linear = -1;
  else if (strcasecmp(cline,"X") == 0)
    linear = 1;
  else if (strcasecmp(cline,"Y") == 0)
    linear = 2;
  else if (strcasecmp(cline,"Z") == 0)
    linear = 3;
  else if (strcasecmp(cline,"XrotZ") == 0)
    linear = 4;
  else if (strcasecmp(cline,"YrotZ") == 0)
    linear = 5;
  else if (strcasecmp(cline,"StepX") == 0)
    linear = 6;
  else if (strcasecmp(cline,"StepXY") == 0)
    linear = 7;
  else if (strcasecmp(cline,"Interpolate1D") == 0)
    linear = 8;
  else if (strcasecmp(cline,"quadX") == 0)
    quadratic = 1;
  else if (strcasecmp(cline,"quadY") == 0)
    quadratic = 2;
  else if (strcasecmp(cline,"quadZ") == 0)
    quadratic = 3;

  Real C = A;
  const RealFunc* f = 0;
  if (linear > 0 && (cline = strtok(NULL," ")))
  {
    C = Real(1);
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
    C = Real(1);
    Real a = atof(cline);
    Real b = atof(strtok(NULL," "));
    Real val = (a-b)*(a-b)/Real(4);
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

  std::cout <<" * ";
  const ScalarFunc* s = parseTimeFunc(cline,NULL,C);

  if (f)
    return new SpaceTimeFunc(f,s);
  else
    return new ConstTimeFunc(s);
}


const ScalarFunc* utl::parseTimeFunc (const char* type, char* cline, Real C)
{
  if (strncasecmp(type,"expr",4) == 0 && cline != NULL)
  {
    std::cout << cline;
    return new EvalFunc(cline,"t");
  }
  else if (strncasecmp(type,"Ramp",4) == 0 || strcmp(type,"Tinit") == 0)
  {
    Real xmax = atof(strtok(cline," "));
    std::cout <<"Ramp(t,"<< xmax <<")";
    return new RampFunc(C,xmax);
  }
  else if (strncasecmp(type,"Dirac",5) == 0)
  {
    Real xmax = atof(strtok(cline," "));
    std::cout <<"Dirac(t,"<< xmax <<")";
    return new DiracFunc(C,xmax);
  }
  else if (strncasecmp(type,"Step",4) == 0)
  {
    Real xmax = atof(strtok(cline," "));
    std::cout <<"Step(t,"<< xmax <<")";
    return new StepFunc(C,xmax);
  }
  else if (strcasecmp(type,"sin") == 0)
  {
    Real freq = atof(strtok(cline," "));
    if ((cline = strtok(NULL," ")))
    {
      Real phase = atof(cline);
      std::cout <<"sin("<< freq <<"*t + "<< phase <<")";
      return new SineFunc(C,freq,phase);
    }
    else
    {
      std::cout <<"sin("<< freq <<"*t)";
      return new SineFunc(C,freq);
    }
  }
  else // linear in time
  {
    Real scale = atof(type);
    std::cout << scale <<"*t";
    return new LinearFunc(C*scale);
  }
}


ScalarFunc* utl::parseTimeFunc (const char* func, const std::string& type)
{
  char* cstr = NULL;
  const ScalarFunc* sf = NULL;
  if (type == "expression")
  {
    std::cout <<"(expression) ";
    if (func) cstr = strdup(func);
    sf = parseTimeFunc("expression",cstr);
  }
  else if (type == "linear")
    sf = parseTimeFunc(func);
  else
  {
    if (func) cstr = strdup(func);
    sf = parseTimeFunc(type.c_str(),cstr);
  }
  std::cout << std::endl;
  if (cstr) free(cstr);

  return const_cast<ScalarFunc*>(sf);
}


RealFunc* utl::parseRealFunc (const std::string& func, const std::string& type)
{
  if (func.empty()) return NULL;

  std::cout <<": ";
  Real p = Real(0);
  if (type == "expression")
  {
    std::cout << func;
    return new EvalFunction(func.c_str());
  }
  else if (type == "linear")
  {
    p = atof(func.c_str());
    std::cout << p <<"*t";
    return new ConstTimeFunc(new LinearFunc(p));
  }
  else if (type == "constant" || func.find_first_of("\t ") == std::string::npos)
  {
    p = atof(func.c_str());
    std::cout << p;
    return new ConstFunc(p);
  }

  std::string tmp(func);
  p = atof(strtok(const_cast<char*>(tmp.c_str())," "));
  char* funcType = const_cast<char*>(type.c_str());
  return const_cast<RealFunc*>(parseRealFunc(funcType,p));
}


VecFunc* utl::parseVecFunc (const std::string& func, const std::string& type)
{
  if (func.empty()) return NULL;

  if (type == "expression")
  {
    std::cout <<": "<< func;
    return new VecFuncExpr(func.c_str());
  }
  else if (type == "constant")
  {
    Vec3 v;
    std::string tmp(func);
    char* s = strtok(const_cast<char*>(tmp.c_str())," ");
    for (int i = 0; i < 3 && s; i++, s = strtok(NULL," "))
      v[i] = atof(s);
    std::cout <<": "<< v;
    return new ConstVecFunc(v);
  }

  return NULL;
}


TractionFunc* utl::parseTracFunc (const std::string& func,
				  const std::string& type, int dir)
{
  if (func.empty()) return NULL;

  std::cout <<": ";
  Real p = Real(0);
  const RealFunc* f = 0;
  if (type == "expression")
  {
    std::cout << func;
    f = new EvalFunction(func.c_str());
  }
  else if (type == "linear")
  {
    p = atof(func.c_str());
    f = new ConstTimeFunc(new LinearFunc(p));
    std::cout << p <<"*t";
  }
  else if (type == "constant" || func.find_first_of("\t ") == std::string::npos)
  {
    p = atof(func.c_str());
    std::cout << p;
  }
  else
  {
    std::string tmp(func);
    p = atof(strtok(const_cast<char*>(tmp.c_str())," "));
    f = parseRealFunc(const_cast<char*>(type.c_str()),p);
  }

  return f ? new PressureField(f,dir) : new PressureField(p,dir);
}
