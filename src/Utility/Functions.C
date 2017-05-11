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
#include "ExprFunctions.h"
#include "FieldFunctions.h"
#include "Vec3Oper.h"
#include "IFEM.h"
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


Real SpaceTimeFunc::deriv (const Vec3& X, int dir) const
{
  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  return sfunc->deriv(X,dir) * (*tfunc)(Xt ? Xt->t : Real(0));
}


Real SpaceTimeFunc::dderiv (const Vec3& X, int dir1, int dir2) const
{
  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  return sfunc->dderiv(X,dir1,dir2) * (*tfunc)(Xt ? Xt->t : Real(0));
}


Real LinearXFunc::evaluate (const Vec3& X) const
{
  return a*X.x + b;
}


Real LinearXFunc::deriv (const Vec3&, int dir) const
{
  return dir == 1 ? a : Real(0);
}


Real LinearYFunc::evaluate (const Vec3& X) const
{
  return a*X.y + b;
}


Real LinearYFunc::deriv (const Vec3&, int dir) const
{
  return dir == 2 ? a : Real(0);
}


Real LinearZFunc::evaluate (const Vec3& X) const
{
  return a*X.z + b;
}


Real LinearZFunc::deriv (const Vec3&, int dir) const
{
  return dir == 3 ? a : Real(0);
}


Real QuadraticXFunc::evaluate (const Vec3& X) const
{
  Real val = (a-b)/Real(2);
  return max*(a-X.x)*(X.x-b)/(val*val);
}


Real QuadraticXFunc::deriv (const Vec3& X, int dir) const
{
  if (dir != 1) return Real(0);

  Real val = (a-b)/Real(2);
  return max*(a+b-Real(2)*X.x)/(val*val);
}


Real QuadraticXFunc::dderiv (const Vec3&, int dir1, int dir2) const
{
  if (dir1 != 1 || dir2 != 1) return Real(0);

  Real val = (a-b)/Real(2);
  return -max*Real(2)/(val*val);
}


Real QuadraticYFunc::evaluate (const Vec3& X) const
{
  Real val = (a-b)/Real(2);
  return max*(a-X.y)*(X.y-b)/(val*val);
}


Real QuadraticYFunc::deriv (const Vec3& X, int dir) const
{
  if (dir != 2) return Real(0);

  Real val = (a-b)/Real(2);
  return max*(a+b-Real(2)*X.y)/(val*val);
}


Real QuadraticYFunc::dderiv (const Vec3&, int dir1, int dir2) const
{
  if (dir1 != 2 || dir2 != 2) return Real(0);

  Real val = (a-b)/Real(2);
  return -max*Real(2)/(val*val);
}


Real QuadraticZFunc::evaluate (const Vec3& X) const
{
  Real val = (a-b)/Real(2);
  return max*(a-X.z)*(X.z-b)/(val*val);
}


Real QuadraticZFunc::dderiv (const Vec3&, int dir1, int dir2) const
{
  if (dir1 != 3 || dir2 != 3) return Real(0);

  Real val = (a-b)/Real(2);
  return -max*Real(2)/(val*val);
}


Real QuadraticZFunc::deriv (const Vec3& X, int dir) const
{
  if (dir != 3) return Real(0);

  Real val = (a-b)/Real(2);
  return max*(a+b-Real(2)*X.z)/(val*val);
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


Real LinearRotZFunc::deriv (const Vec3& X, int dir) const
{
  if (dir > 2) return Real(0);

  // Always return zero if the argument has no time component
  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  if (!Xt) return Real(0);

  if (dir == 1)
    return rX ?  cos(A*Xt->t) - Real(1) : sin(A*Xt->t);
  else if (dir == 2)
    return rX ? -sin(A*Xt->t) : cos(A*Xt->t) - Real(1);
  else
    return Real(0);
}


Real StepXFunc::evaluate (const Vec3& X) const
{
  return X.x < x0 || X.x > x1 ? Real(0) : fv;
}


Real StepXYFunc::evaluate (const Vec3& X) const
{
  return X.x < x0 || X.x > x1 || X.y < y0 || X.y > y1 ? Real(0) : fv;
}


Interpolate1D::Interpolate1D (const char* file, int dir_, int col, Real ramp) :
  dir(dir_), time(ramp)
{
  std::ifstream is(file);
  if (!is)
  {
    std::cerr <<" *** Interpolate1D: Failed to open file "<< file
              <<", function will be identically zero."<< std::endl;
    return;
  }

  while (is.good() && !is.eof())
  {
    char temp[1024];
    is.getline(temp,1024);
    if (is.eof()) return;
    if (temp[0] == '#') continue;
    std::stringstream str(temp);
    Real x, v;
    str >> x >> v;
    if (col < 2)
      std::swap(x,v);
    else for (int i = 2; i < col; i++)
      str >> v;
    grid.push_back(x);
    values.push_back(v);
  }
}


Real Interpolate1D::evaluate (const Vec3& X) const
{
  if (grid.empty())
    return Real(0);
  else if (grid.size() == 1 || X[dir] <= grid.front())
    return values.front();

  std::vector<Real>::const_iterator xb = std::lower_bound(grid.begin(),
                                                          grid.end(),X[dir]);
  if (xb == grid.end())
    return values.back();

  size_t pos = xb - grid.begin();
  Real x1 = *(xb-1);
  Real x2 = *xb;
  Real v1 = values[pos-1];
  Real v2 = values[pos];
  Real res = v1 + (v2-v1)*(X[dir]-x1)/(x2-x1);
  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  if (Xt && time > Real(0) && Xt->t < time)
    res *= Xt->t/time;

  return res;
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
  else if (strcasecmp(cline,"Field") == 0)
    linear = 9;
  else if (strcasecmp(cline,"quadX") == 0)
    quadratic = 1;
  else if (strcasecmp(cline,"quadY") == 0)
    quadratic = 2;
  else if (strcasecmp(cline,"quadZ") == 0)
    quadratic = 3;

  Real C = A;
  const RealFunc* f = 0;
  if (linear > 0 && (cline = strtok(nullptr," ")))
  {
    C = Real(1);
    IFEM::cout <<"("<< A <<"*";
    if (linear < 4)
      IFEM::cout << char('W' + linear) <<" + "<< cline <<")";
    else if (linear < 6)
      IFEM::cout << char('W' + linear-3) <<"RotZ("<< cline <<"))";
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
      f = new LinearRotZFunc(true,A,atof(cline),atof(strtok(nullptr," ")));
      break;
    case 5:
      f = new LinearRotZFunc(false,A,atof(cline),atof(strtok(nullptr," ")));
      break;
    case 6:
      {
        double x0 = atof(cline);
        double x1 = atof(strtok(nullptr," "));
        IFEM::cout <<"StepX("<< x0 <<","<< x1 <<"))";
        f = new StepXFunc(A,x0,x1);
      }
      break;
    case 7:
      {
        double x0 = atof(cline);
        double y0 = atof(strtok(nullptr," "));
        cline = strtok(nullptr," ");
        if (cline && cline[0] == 't')
        {
          double x1 = atof(strtok(nullptr," "));
          double y1 = atof(strtok(nullptr," "));
          IFEM::cout <<"StepXY(["<< x0 <<","<< x1
                     <<"]x["<< y0 <<","<< y1 <<"]))";
          f = new StepXYFunc(A,x1,y1,x0,y0);
        }
        else
        {
          IFEM::cout <<"StepXY([-inf,"<< x0 <<"]x[-inf,"<< y0 <<"]))";
          f = new StepXYFunc(A,x0,y0);
        }
      }
      break;
    case 8:
      {
        int dir = atoi(strtok(nullptr," ")), col = 2;
        IFEM::cout <<"Interpolate1D("<< cline;
        const char* t = strtok(nullptr," ");
        if (t && t[0] == 'c')
        {
          col = atoi(t+1);
          t = strtok(nullptr," ");
          IFEM::cout <<",column #"<< col;
        }
        IFEM::cout <<","<< (char)('X'+dir);
        if (t)
        {
          double time = atof(t);
          IFEM::cout <<")*Ramp("<< time;
          f = new Interpolate1D(cline,dir,col,time);
        }
        else
          f = new Interpolate1D(cline,dir,col);
        IFEM::cout <<")";
      }
      break;
    case 9:
      {
        std::string basis, field;
        basis = strtok(nullptr, " ");
        field = strtok(nullptr, " ");
        IFEM::cout <<"Field("<< cline <<","<< basis <<","<< field <<")";
        f = new FieldFunction(cline,basis,field);
      }
      break;
    }
    if (cline && (linear != 7 || cline[0] == 't'))
      cline = strtok(nullptr," ");
  }
  else if (quadratic > 0 && (cline = strtok(nullptr," ")))
  {
    C = Real(1);
    Real a = atof(cline);
    Real b = atof(strtok(nullptr," "));
    Real val = (a-b)*(a-b)/Real(4);
    char var = 'W' + quadratic;
    IFEM::cout << A/val <<" * ("<< a <<"-"<< var <<")*("<< b <<"-"<< var <<")";
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
    cline = strtok(nullptr," ");
  }
  else // constant in space
  {
    IFEM::cout << C;
    if (linear < 0)
      f = new ConstFunc(C);
  }

  // Check for time variation
  if (!cline) return f; // constant in time

  IFEM::cout <<" * ";
  const ScalarFunc* s = parseTimeFunction(cline,nullptr,C);

  if (f)
    return new SpaceTimeFunc(f,s);
  else
    return new ConstTimeFunc(s);
}


const ScalarFunc* utl::parseTimeFunction (const char* type, char* cline, Real C)
{
  if (strncasecmp(type,"expr",4) == 0 && cline != nullptr)
  {
    IFEM::cout << cline;
    EvalFunc::numError = 0;
    ScalarFunc* sf = new EvalFunc(cline,"t");
    if (EvalFunc::numError > 0)
    {
      delete sf;
      sf = nullptr;
    }
    return sf;
  }
  else if (strncasecmp(type,"Ramp",4) == 0 || strcmp(type,"Tinit") == 0)
  {
    Real xmax = atof(strtok(cline," "));
    IFEM::cout <<"Ramp(t,"<< xmax <<")";
    return new RampFunc(C,xmax);
  }
  else if (strncasecmp(type,"Dirac",5) == 0)
  {
    Real xmax = atof(strtok(cline," "));
    IFEM::cout <<"Dirac(t,"<< xmax <<")";
    return new DiracFunc(C,xmax);
  }
  else if (strncasecmp(type,"Step",4) == 0)
  {
    Real xmax = atof(strtok(cline," "));
    IFEM::cout <<"Step(t,"<< xmax <<")";
    return new StepFunc(C,xmax);
  }
  else if (strcasecmp(type,"sin") == 0)
  {
    Real freq = atof(strtok(cline," "));
    if ((cline = strtok(nullptr," ")))
    {
      Real phase = atof(cline);
      IFEM::cout <<"sin("<< freq <<"*t + "<< phase <<")";
      return new SineFunc(C,freq,phase);
    }
    else
    {
      IFEM::cout <<"sin("<< freq <<"*t)";
      return new SineFunc(C,freq);
    }
  }
  else // linear in time
  {
    Real scale = atof(type);
    IFEM::cout << scale <<"*t";
    return new LinearFunc(C*scale);
  }
}


ScalarFunc* utl::parseTimeFunc (const char* func, const std::string& type)
{
  char* cstr = nullptr;
  const ScalarFunc* sf = nullptr;
  if (type == "expression")
  {
    IFEM::cout <<"(expression) ";
    if (func) cstr = strdup(func);
    sf = parseTimeFunction("expression",cstr);
  }
  else if (type == "linear")
    sf = parseTimeFunction(func,cstr);
  else
  {
    if (func) cstr = strdup(func);
    sf = parseTimeFunction(type.c_str(),cstr);
  }
  IFEM::cout << std::endl;
  if (cstr) free(cstr);

  return const_cast<ScalarFunc*>(sf);
}


RealFunc* utl::parseRealFunc (const std::string& func, const std::string& type)
{
  if (func.empty()) return nullptr;

  IFEM::cout <<": ";
  Real p = Real(0);
  if (type == "expression")
  {
    IFEM::cout << func;
    EvalFunc::numError = 0;
    RealFunc* rf = new EvalFunction(func.c_str());
    if (EvalFunc::numError > 0)
    {
      delete rf;
      rf = nullptr;
    }
    return rf;
  }
  else if (type == "linear")
  {
    p = atof(func.c_str());
    IFEM::cout << p <<"*t";
    return new ConstTimeFunc(new LinearFunc(p));
  }
  else if (type == "constant" || func.find_first_of("\t ") == std::string::npos)
  {
    p = atof(func.c_str());
    IFEM::cout << p;
    return new ConstFunc(p);
  }

  std::string tmp(func);
  p = atof(strtok(const_cast<char*>(tmp.c_str())," "));
  char* funcType = const_cast<char*>(type.c_str());
  return const_cast<RealFunc*>(parseRealFunc(funcType,p));
}


VecFunc* utl::parseVecFunc (const std::string& func, const std::string& type,
                            const std::string& variables)
{
  if (func.empty()) return nullptr;

  if (type == "expression")
  {
    IFEM::cout <<": "<< func;
    EvalFunc::numError = 0;
    VecFunc* vf = new VecFuncExpr(func,variables);
    if (EvalFunc::numError > 0)
    {
      delete vf;
      vf = nullptr;
    }
    return vf;
  }
  else if (type == "constant")
  {
    Vec3 v;
    std::string tmp(func);
    char* s = strtok(const_cast<char*>(tmp.c_str())," ");
    for (int i = 0; i < 3 && s; i++, s = strtok(nullptr," "))
      v[i] = atof(s);
    IFEM::cout <<": "<< v;
    return new ConstVecFunc(v);
  }

  return nullptr;
}


TractionFunc* utl::parseTracFunc (const std::string& func,
                                  const std::string& type, int dir)
{
  if (func.empty()) return nullptr;

  IFEM::cout <<": ";
  Real p = Real(0);
  const RealFunc* f = 0;
  if (type == "expression")
  {
    IFEM::cout << func;
    EvalFunc::numError = 0;
    f = new EvalFunction(func.c_str());
    if (EvalFunc::numError > 0)
    {
      delete f;
      return nullptr;
    }
  }
  else if (type == "linear")
  {
    p = atof(func.c_str());
    f = new ConstTimeFunc(new LinearFunc(p));
    IFEM::cout << p <<"*t";
  }
  else if (type == "constant" || func.find_first_of("\t ") == std::string::npos)
  {
    p = atof(func.c_str());
    IFEM::cout << p;
  }
  else
  {
    std::string tmp(func);
    p = atof(strtok(const_cast<char*>(tmp.c_str())," "));
    f = parseRealFunc(const_cast<char*>(type.c_str()),p);
  }

  return f ? new PressureField(f,dir) : new PressureField(p,dir);
}
