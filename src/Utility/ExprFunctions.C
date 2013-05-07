// $Id$
//==============================================================================
//!
//! \file ExprFunctions.C
//!
//! \date Dec 1 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Expression function implementations.
//!
//==============================================================================

#include "ExprFunctions.h"
#include "Vec3.h"
#include "Tensor.h"
#include "expreval.h"


EvalFunc::EvalFunc (const char* function, const char* x)
{
  try {
    expr = new ExprEval::Expression;
    f = new ExprEval::FunctionList;
    v = new ExprEval::ValueList;
    f->AddDefaultFunctions();
    v->AddDefaultValues();
    v->Add(x,0,false);
    expr->SetFunctionList(f);
    expr->SetValueList(v);
    expr->Parse(function);
    arg = v->GetAddress(x);
#ifdef USE_OPENMP
    omp_init_lock(&lock);
#endif
  } catch(...) {
    std::cerr <<" *** Error parsing function: "<< function << std::endl;
  }
}


EvalFunc::~EvalFunc ()
{
  delete expr;
  delete f;
  delete v;
#ifdef USE_OPENMP
  omp_destroy_lock(&lock);
#endif
}


Real EvalFunc::evaluate (const Real& x) const
{
#ifdef USE_OPENMP
  omp_set_lock(const_cast<omp_lock_t*>(&lock));
#endif
  Real result = Real(0);
  try {
    *arg = x;
    result = expr->Evaluate();
  } catch(...) {
    std::cerr <<" *** Error evaluating expression function."<< std::endl;
  }
#ifdef USE_OPENMP
  omp_unset_lock(const_cast<omp_lock_t*>(&lock));
#endif

  return result;
}


EvalFunction::EvalFunction (const char* function)
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
    v->Add("t",0,false);
    expr->SetFunctionList(f);
    expr->SetValueList(v);
    expr->Parse(function);
    x = v->GetAddress("x");
    y = v->GetAddress("y");
    z = v->GetAddress("z");
    t = v->GetAddress("t");
#ifdef USE_OPENMP
    omp_init_lock(&lock);
#endif
  } catch(...) {
    std::cerr <<" *** Error parsing function: "<< function << std::endl;
  }

  // Checking if the expression is time-independent
  // Note, this will also catch tings like tan(x), but...
  std::string expr(function);
  IAmConstant = expr.find_first_of('t') > expr.size();
}


EvalFunction::~EvalFunction ()
{
  delete expr;
  delete f;
  delete v;
#ifdef USE_OPENMP
  omp_destroy_lock(&lock);
#endif
}


Real EvalFunction::evaluate (const Vec3& X) const
{
  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
#ifdef USE_OPENMP
  omp_set_lock(const_cast<omp_lock_t*>(&lock));
#endif
  Real result = Real(0);
  try {
    *x = X.x;
    *y = X.y;
    *z = X.z;
    *t = Xt ? Xt->t : Real(0);
    result = expr->Evaluate();
  } catch(...) {
    std::cerr <<" *** Error evaluating expression function."<< std::endl;
  }
#ifdef USE_OPENMP
  omp_unset_lock(const_cast<omp_lock_t*>(&lock));
#endif

  return result;
}


template<>
Vec3 VecFuncExpr::evaluate (const Vec3& X) const
{
  Vec3 result;
  for (size_t i = 0; i < 3 && i < p.size(); ++i)
    result[i] = (*p[i])(X);

  return result;
}


template<>
Tensor TensorFuncExpr::evaluate (const Vec3& X) const
{
  int nsd = p.size() > 8 ? 3 : (p.size() > 3 ? 2 : (p.size() > 0 ? 1 : 0));
  Tensor sigma(nsd);

  int i, j, k = 0;
  for (i = 1; i <= nsd; ++i)
    for (j = 1; j <= nsd; ++j)
      sigma(i,j) = (*p[k++])(X);

  return sigma;
}


template<>
SymmTensor STensorFuncExpr::evaluate (const Vec3& X) const
{
  int nsd = p.size() > 5 ? 3 : (p.size() > 2 ? 2 : (p.size() > 0 ? 1 : 0));
  SymmTensor sigma(nsd,p.size()==4);

  int i, j, k = 0;
  for (i = 1; i <= nsd; ++i)
    for (j = i; j <= nsd; ++j)
      sigma(i,j) = (*p[k++])(X);

  if (p.size() == 4)
    sigma(3,3) = (*p[3])(X);

  return sigma;
}
