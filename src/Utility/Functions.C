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
#include "Vec3.h"
#include <cstring>
#include <fstream>
#include <algorithm>
#include "Tensor.h"

#include "expreval.h"


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
        std::cout << "Interpolate1D(" << cline << "," << (char)('X'+dir) << ")";
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


Interpolate1D::Interpolate1D(const char* file, int dir_) : dir(dir_)
{
  std::ifstream is(file);

  while (is.good()) {
    double x, v;
    is >> x >> v;
    grid.push_back(x);
    values.push_back(v);
  }
}

real Interpolate1D::evaluate(const Vec3& X) const
{
  double x = X[dir];
  std::vector<double>::const_iterator xb =
    find_if(grid.begin(),grid.end()-1,std::bind2nd(std::greater<double>(),x));

  size_t pos = xb-grid.begin();
  double x1 = *(xb-1);
  double x2 = *xb;
  double val1 = values[pos-1];
  double val2 = values[pos];
  double delta = x2-x1;

  if (fabs(delta) < 1.e-8)
    return val1;

  double w1 = (x-x1)/delta;
  double w2 = (x2-x)/delta;

  return w1*val1+w2*val2;
}


EvalFunction::EvalFunction(const char* function)
{
  try {
    expr = new ExprEval::Expression;
    f = new ExprEval::FunctionList;
    v = new ExprEval::ValueList;
    f->AddDefaultFunctions();
    v->AddDefaultValues();
    v->Add("x",0,false);
    v->Add("y",0,false);
    v->Add("z",0,false);
    expr->SetFunctionList(f);
    expr->SetValueList(v);
    expr->Parse(function);
    x = v->GetAddress("x");
    y = v->GetAddress("y");
    z = v->GetAddress("z");
#ifdef USE_OPENMP
    omp_init_lock(&lock);
#endif
  } catch(...) {
    std::cerr <<" *** Error parsing function: " << function << std::endl;
  }
}

EvalFunction::~EvalFunction()
{
  delete expr;
  delete f;
  delete v;
#ifdef USE_OPENMP
  omp_destroy_lock(&lock);
#endif
}

real EvalFunction::evaluate(const Vec3& X) const
{
#ifdef USE_OPENMP
  omp_set_lock(const_cast<omp_lock_t*>(&lock));
#endif
  double result=0.f;
  try {
    *x = X.x;
    *y = X.y;
    *z = X.z;
    result = expr->Evaluate();
  } catch(...) {
    std::cerr << "Error evaluating function" << std::endl;
  }
#ifdef USE_OPENMP
  omp_unset_lock(const_cast<omp_lock_t*>(&lock));
#endif

  return result;
}

  template<>
Vec3 EvalMultiFunction<VecFunc,Vec3,Vec3>::evaluate(const Vec3& X) const
{
  Vec3 result;
  for (size_t i = 0; i < 3 && i < p.size(); ++i)
    result[i] = (*p[i])(X);

  return result;
}

  template<>
SymmTensor EvalMultiFunction<STensorFunc,Vec3,SymmTensor>::evaluate(const Vec3& X) const
{
  int nsd = p.size() > 5 ? 3 : (p.size() > 2 ? 2 : (p.size() > 0 ? 1 : 0));
  SymmTensor sigma(nsd,p.size()==4);
  int k=0;
  for (int i = 1; i <= nsd; ++i)
    for (int j = i; j<= nsd; ++j)
      sigma(i,j) = (*p[k++])(X);
  if (p.size() == 4)
    sigma(3,3) = (*p[3])(X);

  return sigma;
}

  template<>
Tensor EvalMultiFunction<TensorFunc,Vec3,Tensor>::evaluate(const Vec3& X) const
{
  int nsd = sqrt(p.size());
  if (nsd > 3) nsd = 3;
  Tensor sigma(nsd);
  int k=0;
  for (int i = 1; i <= nsd; ++i)
    for (int j = 1; j<= nsd; ++j)
      sigma(i,j) = (*p[k++])(X);

  return sigma;
}
