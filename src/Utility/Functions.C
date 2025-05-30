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
#include "Chebyshev.h"
#include "ExprFunctions.h"
#include "FieldFunctions.h"
#include "TensorFunction.h"
#include "TractionField.h"
#include "StringUtils.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "IFEM.h"
#include "tinyxml2.h"
#include <cstring>
#include <fstream>
#include <sstream>
#include <algorithm>


static const Real zTol = Real(1.0e-12); //!< Zero tolerance on function values

//! \brief Creates a scalar function by parsing a character string.
static const ScalarFunc* parseFunction(const char* type, char* cline,
                                       Real C = Real(1), bool print = true,
                                       bool parseConstant = true);


PressureField::PressureField (Real p, int dir) : pdir(dir), pdfn(nullptr)
{
  pressure = p > zTol || p < -zTol ? new ConstFunc(p) : nullptr;
}


LinearFunc::LinearFunc (const char* file, int c, Real s) : scale(s)
{
  std::ifstream is(file);
  if (!is)
  {
    std::cerr <<"\n *** LinearFunc: Failed to open file "<< file
              <<", function will be f(x) = "<< scale <<"*x"<< std::endl;
    return;
  }

  char temp[1024];
  while (is.good() && is.getline(temp,1024))
  {
    if (temp[0] == '#') continue;
    std::stringstream str(temp);
    Real x, v;
    str >> x >> v;
    if (c < 2)
      std::swap(x,v);
    else for (int i = 2; i < c; i++)
      str >> v;
    if (fvals.empty() || fvals.back().first <= x)
      fvals.push_back({x,v});
    else if (fvals.back().first-x < zTol*fvals.back().first)
    {
      x = 0.5*(fvals.back().first+x);
      fvals.back().first = x;
      fvals.push_back({x,v});
    }
    else
    {
      std::cerr <<"\n *** LinearFunc: x-values aren't monotonically increasing";
      for (size_t i = 0; i < fvals.size(); i++)
        if (i+5 == fvals.size())
          std::cerr <<"\n           :";
        else if (i+5 > fvals.size())
          std::cerr <<"\n     Line "<< i+1 <<": "<< fvals[i].first;
      std::cerr <<"\n     Line "<< fvals.size()+1 <<": "<< x
                <<"\n\n     Only the first "<< fvals.size()
                <<" points will be used."<< std::endl;
      return;
    }
  }
}


LinearFunc::LinearFunc (const std::vector<Real>& x, const std::vector<Real>& y)
{
  scale = Real(1);
  fvals.reserve(x.size());
  for (size_t i = 0; i < x.size() && i < y.size(); i++)
    fvals.push_back({x[i],y[i]});
}


size_t LinearFunc::locate (Real x) const
{
  if (fvals.size() < 2)
    return 0;

  auto&& compArgs = [](const Point& a, Real b) { return a.first < b; };
  return std::lower_bound(fvals.begin(),fvals.end(),x,compArgs) - fvals.begin();
}


bool LinearFunc::isZero () const
{
  for (const Point& v : fvals)
    if (v.second > zTol || v.second < -zTol)
      return false;

  return scale <= zTol && scale >= -zTol;
}


Real LinearFunc::deriv (Real x) const
{
  if (fvals.empty())
    return scale;

  size_t ix = this->locate(x);
  if (ix == 0 || ix >= fvals.size())
    return Real(0);

  // Need to interpolate
  Real x1 = fvals[ix-1].first;
  Real x2 = fvals[ix  ].first;
  Real f1 = fvals[ix-1].second;
  Real f2 = fvals[ix  ].second;
  return scale*(f2-f1)/(x2-x1);
}


Real LinearFunc::evaluate (const Real& x) const
{
  if (fvals.empty())
    return scale*x;

  size_t ix = this->locate(x);
  if (ix == 0)
    return scale*fvals.front().second;
  else if (ix >= fvals.size())
    return scale*fvals.back().second;

  // Need to interpolate
  Real x1 = fvals[ix-1].first;
  Real x2 = fvals[ix  ].first;
  Real f1 = fvals[ix-1].second;
  Real f2 = fvals[ix  ].second;
  return scale*(f1 + (f2-f1)*(x-x1)/(x2-x1));
}


LinVecFunc::LinVecFunc (const char* file, int c)
{
  std::ifstream is(file);
  if (!is)
  {
    std::cerr <<"\n *** LinVecFunc: Failed to open file "<< file
              <<", function will be identically zero."<< std::endl;
    return;
  }
  else if (c < 2)
  {
    std::cerr <<"\n *** LinVecFunc: Column index ("<< c <<") should be > 1"
              <<", function will be identically zero."<< std::endl;
    return;
  }

  char temp[1024];
  while (is.good() && is.getline(temp,1024))
  {
    if (temp[0] == '#') continue;
    std::stringstream str(temp);
    Real x;
    Vec3 v;
    str >> x;
    for (int i = 2; i < c; i++) str >> v.x;
    str >> v;
    if (fvals.empty() || fvals.back().first <= x)
      fvals.push_back({x,v});
    else if (fvals.back().first-x < zTol*fvals.back().first)
    {
      x = 0.5*(fvals.back().first+x);
      fvals.back().first = x;
      fvals.push_back({x,v});
    }
    else
    {
      std::cerr <<"\n *** LinVecFunc: x-values aren't monotonically increasing";
      for (size_t i = 0; i < fvals.size(); i++)
        if (i+5 == fvals.size())
          std::cerr <<"\n           :";
        else if (i+5 > fvals.size())
          std::cerr <<"\n     Line "<< i+1 <<": "<< fvals[i].first;
      std::cerr <<"\n     Line "<< fvals.size()+1 <<": "<< x
                <<"\n\n     Only the first "<< fvals.size()
                <<" points will be used."<< std::endl;
      return;
    }
  }
}


bool LinVecFunc::isZero () const
{
  for (const Point& v : fvals)
    if (!v.second.isZero(zTol))
      return false;

  return true;
}


Vec3 LinVecFunc::evaluate (const Real& x) const
{
  if (fvals.empty())
    return Vec3();

  size_t ix = 0;
  if (fvals.size() > 1)
  {
    auto&& compArgs = [](const Point& a, Real b) { return a.first < b; };
    ix = std::lower_bound(fvals.begin(),fvals.end(),x,compArgs) - fvals.begin();
  }

  if (ix == 0)
    return fvals.front().second;
  else if (ix >= fvals.size())
    return fvals.back().second;

  // Need to interpolate
  Real x1 = fvals[ix-1].first;
  Real x2 = fvals[ix  ].first;
  Vec3 f1 = fvals[ix-1].second;
  Vec3 f2 = fvals[ix  ].second;
  return f1 + (f2-f1)*((x-x1)/(x2-x1));
}


Real RampFunc::evaluate (const Real& x) const
{
  return x < xmax ? fval*x/xmax : fval;
}


Real RampFunc::deriv (Real x) const
{
  return x < xmax ? fval/xmax : Real(0);
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


Real SineFunc::deriv (Real x) const
{
  return freq*scale*cos(freq*x+phase);
}


Real ConstTimeFunc::evaluate (const Vec3& X) const
{
  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  return (*tfunc)(Xt ? Xt->t : Real(0));
}


Real ConstTimeFunc::deriv (const Vec3& X, int dir) const
{
  if (dir < 4) return Real(0);

  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  return tfunc->deriv(Xt ? Xt->t : Real(0));
}


Real SpaceTimeFunc::evaluate (const Vec3& X) const
{
  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  return (*sfunc)(X) * (*tfunc)(Xt ? Xt->t : Real(0));
}


Real SpaceTimeFunc::deriv (const Vec3& X, int dir) const
{
  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  if (dir < 4)
    return sfunc->deriv(X,dir) * (*tfunc)(Xt ? Xt->t : Real(0));
  else
    return (*sfunc)(X) * tfunc->deriv(Xt ? Xt->t : Real(0));
}


Real SpaceTimeFunc::dderiv (const Vec3& X, int dir1, int dir2) const
{
  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  if (dir1 < 4 && dir2 < 4)
    return sfunc->dderiv(X,dir1,dir2) * (*tfunc)(Xt ? Xt->t : Real(0));
  else if (dir1 < 4)
    return sfunc->deriv(X,dir1) * tfunc->deriv(Xt ? Xt->t : Real(0));
  else if (dir2 < 4)
    return sfunc->deriv(X,dir2) * tfunc->deriv(Xt ? Xt->t : Real(0));
  else
    return Real(0);
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
  return X[d-'X'] < x0 || X[d-'X'] > x1 ? Real(0) : fv;
}


Real StepXYFunc::evaluate (const Vec3& X) const
{
  return X.x < x0 || X.x > x1 || X.y < y0 || X.y > y1 ? Real(0) : fv;
}


Real Interpolate1D::evaluate (const Vec3& X) const
{
  Real res = lfunc(X[dir]);

  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  if (Xt && time > Real(0) && Xt->t < time)
    res *= Xt->t/time;

  return res;
}


Real Interpolate1D::deriv (const Vec3& X, int ddir) const
{
  Real res = Real(0);
  if (ddir == dir+1)
    res = lfunc.deriv(X[dir]);
  else if (ddir > 3)
    res = lfunc(X[dir]);
  else
    return res;

  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  if (Xt && time > Real(0) && Xt->t < time)
    res *= (ddir < 4 ? Xt->t : Real(1))/time;
  else if (ddir > 3)
    res = Real(0);

  return res;
}


/*!
  The functions are assumed on the general form
  \f[ f({\bf X},t) = A * g({\bf X}) * h(t) \f]

  The character string \a cline is assumed to contain first the definition
  of the spatial function \a g( \b X ) and then the time function \a h(t).
  Either of the two components may be omitted, for creating a space-function
  constant in time, or a time function constant in space.

  \note This method invokes strtok() with nullptr as first argument,
  assuming the initial call with a non-null string has been performed outside.
*/

const RealFunc* utl::parseRealFunc (char* cline, Real A, bool print)
{
  // Check for spatial variation
  int linear    = 0;
  int quadratic = 0;
  char stepDir  = 'X';
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
  else if (strcasecmp(cline,"StepX") == 0 ||
           strcasecmp(cline,"StepY") == 0 ||
           strcasecmp(cline,"StepZ") == 0)
  {
    linear = 6;
    stepDir = toupper(cline[4]);
  }
  else if (strcasecmp(cline,"StepXY") == 0)
    linear = 7;
  else if (strcasecmp(cline,"Interpolate1D") == 0)
    linear = 8;
  else if (strcasecmp(cline,"Field") == 0)
    linear = 9;
  else if (strcasecmp(cline,"Chebyshev") == 0)
    linear = 10;
  else if (strcasecmp(cline,"quadX") == 0)
    quadratic = 1;
  else if (strcasecmp(cline,"quadY") == 0)
    quadratic = 2;
  else if (strcasecmp(cline,"quadZ") == 0)
    quadratic = 3;

  Real C = A;
  const RealFunc* f = nullptr;
  if (linear+quadratic > 0)
    cline = strtok(nullptr," "); // get first function parameter

  if (linear > 0 && cline)
  {
    C = Real(1);
    if (print) {
      IFEM::cout <<"("<< A <<"*";
      if (linear < 4)
        IFEM::cout << char('W' + linear) <<" + "<< cline <<")";
      else if (linear < 6)
        IFEM::cout << char('W' + linear-3) <<"RotZ("<< cline <<"))";
    }
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
        if (print)
          IFEM::cout <<"Step"<< stepDir <<"("<< x0 <<","<< x1 <<"))";
        f = new StepXFunc(A,x0,x1,stepDir);
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
          if (print)
            IFEM::cout <<"StepXY(["<< x0 <<","<< x1
                       <<"]x["<< y0 <<","<< y1 <<"]))";
          f = new StepXYFunc(A,x1,y1,x0,y0);
        }
        else
        {
          if (print)
            IFEM::cout <<"StepXY([-inf,"<< x0 <<"]x[-inf,"<< y0 <<"]))";
          f = new StepXYFunc(A,x0,y0);
        }
      }
      break;
    case 8:
      {
        int dir = atoi(strtok(nullptr," ")), col = 2;
        if (print)
          IFEM::cout <<"Interpolate1D("<< cline;
        const char* t = strtok(nullptr," ");
        if (t && t[0] == 'c')
        {
          col = atoi(t+1);
          t = strtok(nullptr," ");
          if (print)
            IFEM::cout <<",column #"<< col;
        }
        if (print)
          IFEM::cout <<","<< (char)('X'+dir);
        if (t)
        {
          double time = atof(t);
          if (print)
            IFEM::cout <<")*Ramp("<< time;
          f = new Interpolate1D(cline,dir,col,time);
        }
        else
          f = new Interpolate1D(cline,dir,col);
        if (print)
          IFEM::cout <<")";
      }
      break;
    case 9:
      {
        std::string basis, field;
        basis = strtok(nullptr," ");
        field = strtok(nullptr," ");
        char* lev = strtok(nullptr, " ");
        int level = 0;
        if (lev) {
          level = atoi(lev);
          if (strstr(lev,"f"))
            level |= FieldFuncBase::FIXED_LEVEL;
        }
        if (print)
          IFEM::cout <<"Field("<< cline <<","<< basis <<","<< field <<")";
        f = new FieldFunction(cline,basis,field,level);
      }
      break;
    case 10:
      {
        if (print)
          IFEM::cout <<"Chebyshev("<< cline <<")";
        f = new ChebyshevFunc(cline, true);
      }
      break;
    }
    if (cline && (linear != 7 || cline[0] == 't'))
      cline = strtok(nullptr," ");
  }
  else if (quadratic > 0 && cline)
  {
    C = Real(1);
    Real a = atof(cline);
    Real b = atof(strtok(nullptr," "));
    Real val = (a-b)*(a-b)/Real(4);
    char var = 'W' + quadratic;
    if (print)
      IFEM::cout << A/val <<" * ("<< a <<"-"<< var
                 <<")*("<< b <<"-"<< var <<")";
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
    if (print)
      IFEM::cout << C;
    if (linear < 0)
      f = new ConstFunc(C);
  }

  // Check for time variation
  if (!cline) return f; // constant in time

  if (print)
    IFEM::cout <<" * ";
  const ScalarFunc* s = parseFunction(cline,nullptr,C,print,false);

  if (f)
    return new SpaceTimeFunc(f,s);
  else
    return new ConstTimeFunc(s);
}


static const ScalarFunc* parseFunction (const char* type, char* cline,
                                        Real C, bool print, bool parseConstant)
{
  if (strncasecmp(type,"expr",4) == 0)
  {
    ScalarFunc* sf = nullptr;
    if ((cline = strtok(cline,":")))
    {
      if (print)
        IFEM::cout << cline;
      EvalFuncScalar<Real>::numError = 0;
      sf = new EvalFunc(cline,"t",C);
      if (EvalFunc::numError > 0)
      {
        delete sf;
        sf = nullptr;
      }
    }
    // The derivative can be specified as a second expression after the colon
    if (sf && (cline = strtok(nullptr,":")))
    {
      if (print)
        IFEM::cout <<" (derivative: "<< cline <<")";
      static_cast<EvalFunc*>(sf)->addDerivative(cline,"t");
    }
    return sf;
  }

  // Lambda function for extracting the constant from the strtok buffer.
  auto&& getConstant = [&cline,print,parseConstant](Real c)
  {
    if (parseConstant && (cline = strtok(nullptr," ")))
      if ((c = atof(cline)) != Real(1) && print)
        IFEM::cout << c <<"*";
    return c;
  };

  if (strncasecmp(type,"Ramp",4) == 0 || strcmp(type,"Tinit") == 0)
  {
    Real xmax = (cline = strtok(cline," ")) ? atof(cline) : Real(1);
    Real Cval = getConstant(C);
    if (print)
      IFEM::cout <<"Ramp(t,"<< xmax <<")";
    return new RampFunc(Cval,xmax);
  }
  else if (strncasecmp(type,"Dirac",5) == 0)
  {
    Real xmax = (cline = strtok(cline," ")) ? atof(cline) : Real(1);
    Real Cval = getConstant(C);
    if (print)
      IFEM::cout <<"Dirac(t,"<< xmax <<")";
    return new DiracFunc(Cval,xmax);
  }
  else if (strncasecmp(type,"Step",4) == 0)
  {
    Real xmax = (cline = strtok(cline," ")) ? atof(cline) : Real(1);
    Real Cval = getConstant(C);
    if (print)
      IFEM::cout <<"Step(t,"<< xmax <<")";
    return new StepFunc(Cval,xmax);
  }
  else if (strcasecmp(type,"sin") == 0)
  {
    Real freq = (cline = strtok(cline," ")) ? atof(cline) : Real(1);
    if ((cline = strtok(nullptr," ")))
    {
      Real phase = atof(cline);
      Real Cval  = getConstant(C);
      if (print)
      {
        IFEM::cout <<"sin("<< freq <<"*t";
        if (phase > 1.0e-16)
          IFEM::cout <<" + "<< phase;
        else if (phase < -1.0-16)
          IFEM::cout <<" - "<< -phase;
        IFEM::cout <<")";
      }
      return new SineFunc(Cval,freq,phase);
    }
    else
    {
      if (print && C != Real(1))
        IFEM::cout << C <<"*";
      if (print)
        IFEM::cout <<"sin("<< freq <<"*t)";
      return new SineFunc(C,freq);
    }
  }
  else if (strncasecmp(type,"PiecewiseLin",12) == 0)
  {
    char* fname = strtok(cline," ");
    if (!fname) return nullptr;
    int colum = (cline = strtok(nullptr," ")) ? atoi(cline) : 2;
    Real Cval = getConstant(C);
    if (print)
      IFEM::cout <<"PiecewiseLin(t,"<< fname <<","<< colum <<")";
    return new LinearFunc(fname,colum,Cval);
  }
  else if (strncasecmp(type,"Constant",8) == 0 && cline)
  {
    Real value = atof(cline);
    if (print)
      IFEM::cout << value;
    return new ConstantFunc(C*value);
  }
  else // linear in time
  {
    Real scale;
    if (cline)
      scale = atof(strncasecmp(type,"Lin",3) == 0 ? cline : type);
    else
      scale = Real(1);

    if (print)
      IFEM::cout << scale <<"*t";
    return new LinearFunc(C*scale);
  }
}


IntFunc* utl::parseIntFunc (const std::string& func, const std::string& type)
{
  class Identity : public IntFunc
  {
    Real s; // Scaling factor

  public:
    Identity(Real sf) : s(sf) {}

  protected:
    Real evaluate(const int& x) const override { return s*Real(x > 1 ? x : 0); }
  };

  // This class can be used as an element activation function, where a bunch of
  // elements (nv) is activated in the first parameter direction, assuming the
  // number of elements in the first direction is nu. If nw equals 1, all nv
  // elements in the v-w parameter plane are activated instantly in one step.
  // Otherwise, the nv elements are divided into nw groups, each consisting
  // of nv/nw elements, which then are activated one by one.
  class Parametric : public IntFunc
  {
    int nu; // Number of elements in first parameter direction
    int nv; // Number of elements to activate in v-w plane before advancing in u
    int nw; // Number of element divisions per layer in w-direction
    Real s; // Scaling factor

  public:
    Parametric(int u, int v, int w, Real sf) : nu(u), nv(v), nw(w), s(sf) {}

  protected:
    Real evaluate(const int& x) const override
    {
      int ix = (x-1)%nu;
      int nx = (x-1)/nu;
      if (ix < 1 && nx < nv) return Real(0);

      Real delta = nw > 1 && nw < nv ? Real(1+nx/(nv/nw))/Real(nw) : Real(0);
      return s*(Real(ix + nu*(nx/nv)) + delta);
    }
  };

  if (type == "parametric" && !func.empty())
  {
    char* prms = strdup(func.c_str());
    char* s = nullptr;
    int ne1 = atoi(strtok(prms," "));
    int ne2 = (s = strtok(nullptr," ")) ? atoi(s) : 1;
    int ne3 = (s = strtok(nullptr," ")) ? atoi(s) : 1;
    Real sf = (s = strtok(nullptr," ")) ? atof(s) : Real(1);
    free(prms);
    return new Parametric(ne1,ne2,ne3,sf);
  }

  Real scaling(1);
  if (func.empty())
    IFEM::cout <<"  ** utl::parseIntFunc: \""<< type
               <<"\" not implemented - returning identity"<< std::endl;
  else
  {
    char* prms = strdup(func.c_str());
    scaling = atof(strtok(prms," "));
    free(prms);
  }
  return new Identity(scaling);
}


ScalarFunc* utl::parseTimeFunc (const char* func, const std::string& type,
                                Real eps)
{
  char* cstr = func ? strdup(func) : nullptr;
  const ScalarFunc* sf = nullptr;
  if (type == "expression")
  {
    IFEM::cout <<"(expression) ";
    sf = parseFunction("expression",cstr,eps);
  }
  else
    sf = parseFunction(type.c_str(),cstr);
  IFEM::cout << std::endl;
  if (cstr) free(cstr);

  return const_cast<ScalarFunc*>(sf);
}


VecTimeFunc* utl::parseVecTimeFunc (const char* func, const std::string& type)
{
  VecTimeFunc* vfunc = nullptr;
  if (strncasecmp(type.c_str(),"PiecewiseLin",12) == 0)
  {
    char* cline = strdup(func);
    char* ctemp = cline;
    char* fname = strtok(cline," ");
    int   colum = (cline = strtok(nullptr," ")) ? atoi(cline) : 2;
    IFEM::cout <<"PiecewiseLin(t,"<< fname <<","<< colum <<")";
    vfunc = new LinVecFunc(fname,colum);
    free(ctemp);
  }
  IFEM::cout << std::endl;

  return vfunc;
}


RealFunc* utl::parseRealFunc (const std::string& func,
                              const std::string& type, bool print)
{
  Real p = Real(0);
  if (func.empty())
    return new ConstFunc(p);

  if (print)
    IFEM::cout <<": ";

  const RealFunc* f = nullptr;
  if (type == "expression")
  {
    if (print)
      IFEM::cout << func;
    EvalFunc::numError = 0;
    f = new EvalFunction(func.c_str());
    if (EvalFunc::numError > 0)
    {
      delete f;
      f = nullptr;
    }
  }
  else if (type == "linear")
  {
    p = atof(func.c_str());
    if (print)
      IFEM::cout << p <<"*t";
    f = new ConstTimeFunc(new LinearFunc(p));
  }
  else if (type == "chebyshev")
  {
    if (print)
      IFEM::cout << func;
    f = new ChebyshevFunc(func, true);
  }
  else if (type == "constant" || func.find_first_of("\t ") == std::string::npos)
  {
    p = atof(func.c_str());
    if (print)
      IFEM::cout << p;
    f = new ConstFunc(p);
  }
  else
  {
    char* funcType = strdup(type.c_str());
    char* funcBody = strdup(func.c_str());
    p = atof(strtok(funcBody," "));
    f = parseRealFunc(funcType,p,print);
    free(funcBody);
    free(funcType);
  }

  return const_cast<RealFunc*>(f);
}


VecFunc* utl::parseVecFunc (const std::string& func, const std::string& type,
                            const std::string& variables)
{
  VecFunc* f = nullptr;

  if (func.empty())
    f = new ConstVecFunc(Vec3());

  else if (type == "expression")
  {
    IFEM::cout <<": "<< func;
    EvalFunc::numError = 0;
    f = new VecFuncExpr(func,variables);
    if (EvalFunc::numError > 0)
    {
      delete f;
      f = nullptr;
    }
  }
  else if (type == "constant")
  {
    Vec3 v;
    char* ff = strdup(func.c_str());
    char* ss = strtok(ff," ");
    for (int i = 0; i < 3 && ss; i++, ss = strtok(nullptr," "))
      v[i] = atof(ss);
    free(ff);
    IFEM::cout <<": "<< v;
    f = new ConstVecFunc(v);
  }
  else if (type == "chebyshev")
    f = new ChebyshevVecFunc(splitString(func),true);

  else if (type.substr(0,5) == "field")
  {
    std::vector<std::string> prms = splitString(func);
    if (prms.size() > 2)
    {
      int level = 0;
      if (prms.size() > 3)
      {
        if (prms[3].find('f') != std::string::npos)
          level = FieldFuncBase::FIXED_LEVEL;
        else
          level = atoi(prms[3].c_str());
      }
      if (type == "fieldgrad")
        f = new ScalarGradFieldFunction(prms[0],prms[1],prms[2],level);
      else if (type == "fieldlaplacian")
        f = new ScalarLaplacianFieldFunction(prms[0],prms[1],prms[2],level);
      else
        f = new VecFieldFunction(prms[0],prms[1],prms[2],level);
    }
  }

  return f;
}


TensorFunc* utl::parseTensorFunc (const std::string& func,
                                  const std::string& type)
{
  TensorFunc* f = nullptr;

  if (type == "chebyshev" && !func.empty())
    f = new ChebyshevTensorFunc(splitString(func),true);

  else if (type.substr(0,5) == "field")
  {
    std::vector<std::string> prms = splitString(func);
    if (prms.size() > 2)
    {
      int level = 0;
      if (prms.size() > 3)
      {
        if (prms[3].find('f') != std::string::npos)
          level = FieldFuncBase::FIXED_LEVEL;
        else
          level = atoi(prms[3].c_str());
      }
      if (type == "fieldgrad")
        f = new VecGradFieldFunction(prms[0],prms[1],prms[2],level);
      else if (type == "fieldlaplacian")
        f = new VecLaplacianFieldFunction(prms[0],prms[1],prms[2],level);
    }
  }

  return f;
}


TractionFunc* utl::parseTracFunc (const std::string& func,
                                  const std::string& type, int dir)
{
  Real p = Real(0);
  RealFunc* f = nullptr;
  if (!func.empty())
  {
    if (type == "constant" || func.find_first_of("\t ") == std::string::npos)
      IFEM::cout <<": "<< (p = atof(func.c_str()));
    else
      f = parseRealFunc(func,type);
  }

  return f ? new PressureField(f,dir) : new PressureField(p,dir);
}


TractionFunc* utl::parseTracFunc (const tinyxml2::XMLElement* elem)
{
  Vec3 X0, Xaxis(1,0,0), Zaxis(0,0,1);
  if (utl::getAttribute(elem,"X0",X0))
    IFEM::cout <<"\n\tX0    = "<< X0;
  if (utl::getAttribute(elem,"Xaxis",Xaxis))
    IFEM::cout <<"\n\tXaxis = "<< Xaxis;
  if (utl::getAttribute(elem,"Zaxis",Zaxis))
    IFEM::cout <<"\n\tZaxis = "<< Zaxis;
  IFEM::cout << std::endl;

  char rotaxis = 'X';
  const ScalarFunc*  force = nullptr;
  const VecTimeFunc* fdir  = nullptr;
  const VecTimeFunc* frot  = nullptr;
  const ScalarFunc*  angle = nullptr;
  const RealFunc*    shape = nullptr;
  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  while (child && child->Value() && child->FirstChild())
  {
    std::string type;
    utl::getAttribute(child,"type",type,true);
    if (strcasecmp(child->Value(),"force") == 0)
    {
      IFEM::cout <<"\tForce resultant: ";
      force = parseTimeFunc(child->FirstChild()->Value(),type);
    }
    else if (strcasecmp(child->Value(),"angle") == 0)
    {
      if (utl::getAttribute(child,"axis",rotaxis,false))
        rotaxis = toupper(rotaxis);
      IFEM::cout <<"\tForce angle ("<< rotaxis <<"-axis): ";
      angle = parseTimeFunc(child->FirstChild()->Value(),type);
    }
    else if (strncasecmp(child->Value(),"orient",6) == 0)
    {
      IFEM::cout <<"\tForce orientation: ";
      frot = parseVecTimeFunc(child->FirstChild()->Value(),type.c_str());
    }
    else if (strcasecmp(child->Value(),"direction") == 0)
    {
      IFEM::cout <<"\tForce direction: ";
      fdir = parseVecTimeFunc(child->FirstChild()->Value(),type.c_str());
    }
    else if (strcasecmp(child->Value(),"shape") == 0)
    {
      IFEM::cout <<"\tShape function";
      shape = parseRealFunc(child->FirstChild()->Value(),type);
      IFEM::cout << std::endl;
    }
    child = child->NextSiblingElement();
  }

  if (angle && rotaxis >= 'X' && rotaxis <= 'Z')
    return new ForceDirField(force,angle,rotaxis,shape,Xaxis,Zaxis,X0);
  else if (fdir)
    return new ForceDirField(force,fdir,shape,Xaxis,Zaxis,X0,true);
  else
    return new ForceDirField(force,frot,shape,Xaxis,Zaxis,X0);
}
