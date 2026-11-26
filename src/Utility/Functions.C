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
#include <iterator>


namespace
{
  constexpr Real zTol = Real(1.0e-12); //!< Zero tolerance on function values

  //! \brief Creates a scalar function by parsing a character string.
  ScalarFunc* parseFunction (const char* type, char* cline, Real C = Real(1),
                             bool print = true, bool parseConstant = true)
  {
    if (strncasecmp(type,"expr",4) == 0)
    {
      ScalarFunc* sf = nullptr;
      if ((cline = strtok(cline,":")))
      {
        if (print)
          IFEM::cout << cline;
        ExprEval::numError = 0;
        sf = new EvalFunc(cline,"t",C);
        if (ExprEval::numError > 0)
        {
          delete sf;
          sf = nullptr;
          if (print)
            IFEM::cout << std::endl;
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
    else if (strncasecmp(type,"Pulse",5) == 0)
    {
      // Get the pulse amplitude
      Real Cval = C;
      char* tmp = cline;
      if (parseConstant)
      {
        if ((cline = strtok(cline," ")))
          Cval = atof(cline);
        tmp = nullptr;
      }

      // Get pulse location (x) and width (w).
      // "Pulse" (same width): w x1 x2 x3 ...
      // "Pulses" (different widths): x1 w1 x2 w2  ...
      utl::RealVec x;
      for (cline = strtok(tmp," "); cline; cline = strtok(nullptr," "))
        x.push_back(atof(cline));
      if (x.size() < 2)
      {
        if (print)
          IFEM::cout << std::endl;
        std::cerr <<"\n *** A pulse function needs at least two values"
                  <<" in addition to "<< Cval <<": "<< x.size() << std::endl;
        return nullptr;
      }

      if (Cval != Real(1) && print)
        IFEM::cout << Cval <<"*";
      if (strncasecmp(type,"Pulses",6) == 0)
      {
        utl::RealVec w(x.size()/2);
        for (size_t i = 0; i < x.size(); i++)
          if (i%2 > 0)
            w[i/2] = x[i];
          else if (i > 0)
            x[i/2] = x[i];
        x.resize(w.size()+x.size()%2);
        if (print)
        {
          IFEM::cout <<"Pulse(t";
          for (size_t i = 0; i < x.size(); i++)
            IFEM::cout <<",("<< x[i] <<","
                       << (i < w.size() ? w[i] : w.back()) <<")";
          IFEM::cout <<")";
        }
        return new PulseFunc(Cval,w,x);
      }
      else
      {
        Real w = x.front();
        x.erase(x.begin());
        if (print)
        {
          IFEM::cout <<"Pulse(t,"<< w <<",("<< x.front();
          for (size_t i = 1; i < x.size(); i++)
            IFEM::cout <<","<< x[i];
          IFEM::cout <<"))";
        }
        return new PulseFunc(Cval,w,x);
      }
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
          if (phase > zTol)
            IFEM::cout <<" + "<< phase;
          else if (phase < -zTol)
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
    else if (cline) // linear in time
    {
      Real scale = atof(strncasecmp(type,"Lin",3) == 0 ? cline : type);

      if (print)
        IFEM::cout << scale <<"*t";
      return new LinearFunc(C*scale);
    }
    else // linear in time
    {
      if (print)
        IFEM::cout <<"t";
      return new LinearFunc(C);
    }
  }

  //! \brief Returns \e true if \a func contains a single real constant.
  bool isConstant (const std::string& func, const std::string& type, Real& p)
  {
    if (type == "constant")
      p = atof(func.c_str());
    else
    {
      char* endPtr = nullptr;
      p = strtod(func.c_str(),&endPtr);
      if (strlen(endPtr) > 0)
        return false; // more data, not a single numerical value
    }
    return true;
  }

  //! \brief Nicely prints an expression over several lines.
  void printExpression (const std::string& func, char sep = ';')
  {
    IFEM::cout <<":";
    if (func.find_first_of(sep) == std::string::npos)
    {
      IFEM::cout <<"\n\t\t"<< func;
      return;
    }

    char* tmp = strdup(func.c_str());
    const char ssep[2] = { sep, '\0' };
    for (char* s = strtok(tmp,ssep); s; s = strtok(nullptr,ssep))
    {
      while (*s && isspace(*s)) ++s;
      IFEM::cout <<"\n\t\t"<< s;
      if (sep == ';')
        IFEM::cout << sep;
    }
    free(tmp);
  }

  //! \brief Helper swapping the domain range if reversed.
  bool reversed (Real& x0, Real& x1)
  {
    if (x1 >= x0)
      return false;

    std::swap(x0,x1);
    return true;
  }
}


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
      x = (fvals.back().first+x)/Real(2);
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


LinearFunc::LinearFunc (const utl::RealVec& x, const utl::RealVec& y)
{
  scale = Real(1);
  fvals.reserve(x.size());
  for (size_t i = 0; i < x.size() && i < y.size(); i++)
    fvals.emplace_back(x[i],y[i]);
}


bool LinearFunc::isZero () const
{
  for (const Point& v : fvals)
    if (v.second > zTol || v.second < -zTol)
      return false;

  return scale <= zTol && scale >= -zTol;
}


size_t LinearFunc::locate (Real x) const
{
  if (fvals.size() < 2)
    return 0;

  return std::distance(fvals.begin(),
                       std::lower_bound(fvals.begin(), fvals.end(), x,
                                        [](const Point& a, Real b)
                                        { return a.first < b; }));
}


Real LinearFunc::interpolate (Real x, size_t i1) const
{
  const size_t i2 = i1 + 1;
  const Real x1 = fvals[i1].first;
  const Real x2 = fvals[i2].first;
  const Real f1 = fvals[i1].second;
  const Real f2 = fvals[i2].second;
  return f1 + (f2-f1)*(x-x1)/(x2-x1);
}


Real LinearFunc::evaluate (const Real& x) const
{
  if (fvals.empty())
    return scale*x;

  size_t ix = this->locate(x);
  if (ix == 0) // flat extrapolation before the first segment
    return scale*fvals.front().second;
  else if (ix >= fvals.size()) // flat extrapolation after the last segment
    return scale*fvals.back().second;
  else // need to interpolate between points ix-1 and ix
    return scale*this->interpolate(x,ix-1);
}


Real LinearFunc::deriv (Real x) const
{
  if (fvals.empty())
    return scale;

  size_t ix = this->locate(x);
  if (ix == 0 || ix >= fvals.size())
    return Real(0); // zero derivative outside the line segments

  // Need to interpolate
  const size_t i0 = ix - 1;
  const Real x1 = fvals[i0].first;
  const Real x2 = fvals[ix].first;
  const Real f1 = fvals[i0].second;
  const Real f2 = fvals[ix].second;
  return scale*(f2-f1)/(x2-x1);
}


Real LinearFunc::integrate (Real x0, Real x1, Real) const
{
  if (fvals.empty())
    return (x1-x0)*scale*(x0+x1)/Real(2);

  const bool negative = reversed(x0,x1);

  Real res = Real(0);
  if (fvals.front().first > x0) // flat section before the first segment
    res += fvals.front().second * (std::min(fvals.front().first,x1)-x0);
  if (fvals.back().first < x1) // flat section after the last segment
    res += fvals.back().second * (x1-std::max(fvals.back().first,x0));

  size_t i0 = this->locate(x0);
  size_t i1 = this->locate(x1);
  for (size_t i = i0; i+1 < i1; i++)
  {
    // Integrate the whole segment
    Real fVal = (fvals[i].second + fvals[i+1].second)/Real(2);
    res += fVal*(fvals[i+1].first - fvals[i].first);
  }

  if (i0 > 0 && i0 < fvals.size() && fvals[i0].first > x0)
  {
    // Integrate the first boundary segment
    Real x2 = std::min(x1,fvals[i0].first);
    Real f2 = x2 < x1 ? fvals[i0].second : this->interpolate(x1,i0-1);
    res += (this->interpolate(x0,i0-1) + f2)*(x2-x0)/Real(2);
  }

  if (i1 > i0 && i1 < fvals.size() && fvals[i1--].first >= x1)
  {
    // Integrate the last boundary segment
    Real x2 = std::max(x0,fvals[i1].first);
    Real f2 = x2 > x0 ? fvals[i1].second : this->interpolate(x0,i1);
    res += (f2 + this->interpolate(x1,i1))*(x1-x2)/Real(2);
  }

  return scale*(negative ? -res : res);
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
      x = (fvals.back().first+x)/Real(2);
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
    ix = std::distance(fvals.begin(),
                       std::lower_bound(fvals.begin(), fvals.end(), x,
                                        [](const Point& a, Real b)
                                        { return a.first < b; }));

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
  return x < xmax && xmax > Real(0) ? fval*x/xmax : fval;
}


Real RampFunc::deriv (Real x) const
{
  return x < xmax && xmax > Real(0) ? fval/xmax : Real(0);
}


Real RampFunc::integrate (Real x0, Real x1, Real) const
{
  if (xmax <= Real(0))
    return fval*(x1-x0);

  const Real A = reversed(x0,x1) ? -fval : fval;

  if (x1 <= xmax)
    return A*(x0+x1)*(x1-x0)/(xmax*Real(2));
  else if (x0 >= xmax)
    return A*(x1-x0);

  return A*(x1-xmax + (x0+xmax)*(xmax-x0)/(xmax*Real(2)));
}


Real DiracFunc::evaluate (const Real& x) const
{
  constexpr Real epsX = Real(1.0e-4);

  return x-xmax < epsX && x-xmax > -epsX ? amp : Real(0);
}


Real StepFunc::evaluate (const Real& x) const
{
  return x >= xmax ? amp : Real(0);
}


Real StepFunc::integrate (Real x0, Real x1, Real) const
{
  const Real A = reversed(x0,x1) ? -amp : amp;

  if (x0 >= xmax)
    return A*(x1-x0);
  else if (x1 <= xmax)
    return Real(0);

  return A*(x1-xmax);
}


Real PulseFunc::evaluate (const Real& x) const
{
  Real res = Real(0);
  for (size_t i = 0; i < xloc.size(); i++)
    if (x >= xloc[i] && x <= xloc[i]+delta[i])
      res += amp;
    else if (x < xloc[i])
      break;

  return res;
}


Real PulseFunc::integrate (Real x0, Real x1, Real) const
{
  const bool negative = reversed(x0,x1);

  Real res = Real(0);
  for (size_t i = 0; i < xloc.size() && xloc[i] < x1; i++)
    if (Real x2 = xloc[i]+delta[i]; x2 > x1)
      res += x1 - std::max(x0,xloc[i]);
    else if (x2 > x0)
      res += x2 - std::max(x0,xloc[i]);

  return amp*(negative ? -res : res);
}


Real SineFunc::evaluate (const Real& x) const
{
  return scale*sin(freq*x+phase);
}


Real SineFunc::deriv (Real x) const
{
  return freq*scale*cos(freq*x+phase);
}


Real SineFunc::integrate (Real x0, Real x1, Real) const
{
  constexpr double zero = Real(1.0e-16);

  if (freq >= -zero && freq <= zero)
    return scale*sin(phase)*(x1-x0);

  const Real freq2 = freq/Real(2);
  // Use identity cos(a) - cos(b) = -2*sin((a+b)/2)*sin(a-b)/2)
  // to avoid round-off cancellation for small frequencies
  return (scale/freq2)*sin(freq2*(x0+x1)+phase)*sin(freq2*(x1-x0));
}


Real ScalarFunc::integrate (Real x0, Real x1, Real dx) const
{
  const bool negative = reversed(x0,x1);
  constexpr double ex = Real(1.0e-8);
  if (dx < ex) dx = ex;

  Real res = Real(0);
  Real dx2 = dx/Real(2);
  Real x   = x0;
  while (x+dx <= x1)
  {
    res += this->evaluate(x+dx2)*dx;
    x += dx;
  }
  if (x < x1)
    res += this->evaluate((x+x1)/Real(2))*(x1-x);

  return negative ? -res : res;
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


Real StepXYZFunc::evaluate (const Vec3& X) const
{
  return X.inside(X0,X1,eps) ? fv : Real(0);
}


/*!
  \class DiracSpaceFunc
  \details This class defines a spatial
  <a href="https://en.wikipedia.org/wiki/Mollifier">mollifier function</a>
  that smooths out a discrete Dirac function around the point #Xmax over a disk
  (or sphere) with radius &epsilon; (&epsilon;<sup>2</sup>=#eps2), or
  (in 1D) over the line segment \f$x\in[x_{max}-\epsilon,x_{max}+\epsilon]\f$.

  The function is defined as follows:
  \f[
    f({\bf x}) = \frac{A}{A_0}\exp(\frac{\epsilon^2}{\epsilon^2-r^2})
  \f]
  where
  - \f$ r=|{\bf x} - {\bf X}_{max}| \f$ and \f${\bf X}_{max}\f$ = #Xmax.
  - \a A is the function amplitude (#amp).
  - \a A<sub>0</sub> is the normalization factor (#fact) ensuring that
    the integral of the function over the entire domain equals \a A.

  The constant \a A<sub>0</sub> (or #fact) is computed numerically as follows:
  \f[
    A_0=\int\limits_{r<\epsilon}\exp{\frac{-\epsilon^2}{\epsilon^2-r^2}}d{\bf x}
    \;\approx\left\{\begin{array}{ll}
      0.44399381616807944\,\epsilon   & \mbox{in 1D} \\
      0.46651239317833007\,\epsilon^2 & \mbox{in 2D} \\
      0.44108888727660440\,\epsilon^3 & \mbox{in 3D}
    \end{array}\right.
  \f]
*/

DiracSpaceFunc::DiracSpaceFunc (Real a, const Vec3& X, Real eps, int n)
  : nsd(n), amp(a), Xmax(X), eps2(eps*eps)
{
  // Define the normalization factor such that the integral equals amp
  if (nsd == 1)
    fact = 0.44399381616807944 * eps;
  else if (nsd == 2)
    fact = 0.46651239317833007 * eps2;
  else
    fact = 0.44108888727660440 * eps2*eps;
}


Real DiracSpaceFunc::evaluate (const Vec3& X) const
{
  Real r2 = Real(0);
  for (int i = 0; i < nsd; i++)
    r2 += (X[i] - Xmax[i]) * (X[i] - Xmax[i]);

  return r2 < eps2 ? amp * exp(-eps2/(eps2-r2)) / fact : Real(0);
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
  Check if the expression is time-dependent by searching for the occurence
  of 't' where neither the next nor the previous characters are alphanumeric.
*/

bool utl::isTimeExpression (const std::string& expr)
{
  auto isAlphaNum = [](char c) -> bool
  {
    return isalnum(static_cast<unsigned char>(c)) || c == '_';
  };

  const size_t n = expr.size();
  for (size_t i = expr.find('t'); i < n; i = expr.find('t',i+1))
  {
    bool ok1 = i   == 0 || !isAlphaNum(expr[i-1]);
    bool ok2 = i+1 == n || !isAlphaNum(expr[i+1]);
    if (ok1 && ok2) return true;
  }

  return false;
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
  RealFunc* f = nullptr;
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
        Real x0 = atof(cline);
        Real x1 = atof(strtok(nullptr," "));
        if (print)
          IFEM::cout <<"Step"<< stepDir <<"("<< x0 <<","<< x1 <<"))";
        f = new StepXFunc(A,x0,x1,stepDir);
      }
      break;
    case 7:
      {
        Real x0 = atof(cline);
        Real y0 = atof(strtok(nullptr," "));
        cline = strtok(nullptr," ");
        if (cline && cline[0] == 't')
        {
          Real x1 = atof(strtok(nullptr," "));
          Real y1 = atof(strtok(nullptr," "));
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
          Real time = atof(t);
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

  if (cline) // check for time variation
  {
    if (print)
      IFEM::cout <<" * ";
    ScalarFunc* s = parseFunction(cline,nullptr,C,print,false);
    if (f && s)
      return new SpaceTimeFunc(f,s);
    else if (s)
      return new ConstTimeFunc(s);
  }

  return f; // constant in time
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
  if (sf) IFEM::cout << std::endl;
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
                              const std::string& type,
                              bool print, bool autodiff)
{
  Real p = Real(0);
  if (func.empty())
    return new ConstFunc(p);

  const RealFunc* f = nullptr;
  if (type == "expression")
  {
    if (print)
    {
      if (func.size() > 75)
        printExpression(func);
      else
        IFEM::cout <<": "<< func;
    }
    f = utl::parseExprRealFunc(func,autodiff);
  }
  else if (type == "linear")
  {
    p = atof(func.c_str());
    if (print)
      IFEM::cout <<": "<< p <<"*t";
    f = new ConstTimeFunc(new LinearFunc(p));
  }
  else if (type == "chebyshev")
  {
    if (print)
      IFEM::cout <<": "<< func;
    f = new ChebyshevFunc(func, true);
  }
  else if (isConstant(func,type,p))
  {
    if (print)
      IFEM::cout <<": "<< p;
    f = new ConstFunc(p);
  }
  else
  {
    if (print)
      IFEM::cout <<": ";
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
                            const std::string& variables, bool autodiff)
{
  VecFunc* f = nullptr;

  if (func.empty())
    f = new ConstVecFunc(Vec3());

  else if (type == "expression")
  {
    if (func.size() > 50)
      printExpression(func,'|');
    else
      IFEM::cout <<": "<< func;
    f = utl::parseExprVecFunc(func,variables,autodiff);
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
  {
    IFEM::cout <<": "<< func;
    f = new ChebyshevVecFunc(splitString(func),true);
  }
  else if (type.substr(0,5) == "field")
    if (std::vector<std::string> prms = splitString(func); prms.size() > 2)
    {
      int level = 0;
      if (prms.size() > 3)
      {
        if (prms[3].find('f') != std::string::npos)
          level = FieldFuncBase::FIXED_LEVEL;
        else
          level = atoi(prms[3].c_str());
      }
      IFEM::cout <<":";
      for (int i = 0; i < 3; i++)
        IFEM::cout <<" "<< prms[i];
      IFEM::cout <<" level="<< level;
      if (type == "fieldgrad")
        f = new ScalarGradFieldFunction(prms[0],prms[1],prms[2],level);
      else if (type == "fieldlaplacian")
        f = new ScalarLaplacianFieldFunction(prms[0],prms[1],prms[2],level);
      else
        f = new VecFieldFunction(prms[0],prms[1],prms[2],level);
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
    if (type == "linear")
      f = parseRealFunc(func,type);
    else if (isConstant(func,type,p))
      IFEM::cout <<": "<< p;
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
  for (const tinyxml2::XMLElement* child = elem->FirstChildElement();
       child && child->Value() && child->FirstChild();
       child = child->NextSiblingElement())
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
  }

  if (angle && rotaxis >= 'X' && rotaxis <= 'Z')
    return new ForceDirField(force,angle,rotaxis,shape,Xaxis,Zaxis,X0);
  else if (fdir)
    return new ForceDirField(force,fdir,shape,Xaxis,Zaxis,X0,true);
  else
    return new ForceDirField(force,frot,shape,Xaxis,Zaxis,X0);
}
