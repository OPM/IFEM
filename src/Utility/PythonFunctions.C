// $Id$
//==============================================================================
//!
//! \file PythonFunctions.C
//!
//! \date Sep 7 2021
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Python function implementations.
//!
//==============================================================================

#include "PythonFunctions.h"

#ifdef HAS_PYTHON

#include <Vec3.h>

#include <pybind11/numpy.h>

std::shared_ptr<InterpreterRAII> pyInterp;


static std::shared_ptr<InterpreterRAII> initInterpreter()
{
  if (!pyInterp)
    pyInterp = std::make_shared<InterpreterRAII>();

  return pyInterp;
}


PythonBaseFunc::PythonBaseFunc(const char* module, const char* params) :
  myInterp(initInterpreter())
{
  myModule = pybind11::module::import(module);
  pybind11::object init = myModule.attr("initialize");
  myInstance = init(params);
  eval = myInstance.attr("eval");
}


Real PythonFunc::evaluate (const Real& x) const
{
  Real d;

  // python code cannot be called on multiple threads - at least not easily.
  // for now serialize it
#pragma omp critical
  {
    auto res = eval(x);
    d = res.cast<double>();
  }
  return d;
}


Real PythonFunction::evaluate (const Vec3& x) const
{
  const Vec4* Xt = dynamic_cast<const Vec4*>(&x);
  std::vector<double> data(x.ptr(), x.ptr()+3);
  if (Xt)
    data.push_back(Xt->t);

  double d;
  // python code cannot be called on multiple threads - at least not easily.
  // for now serialize it
#pragma omp critical
  {
    pybind11::array_t<double> X(data.size(), data.data());
    auto res = eval(X);
    d = res.cast<double>();
  }
  return d;
}


Vec3 PythonVecFunc::evaluate (const Vec3& x) const
{
  const Vec4* Xt = dynamic_cast<const Vec4*>(&x);
  std::vector<double> data(x.ptr(), x.ptr()+3);
  if (Xt)
    data.push_back(Xt->t);

  Vec3 d;

  // python code cannot be called on multiple threads - at least not easily.
  // for now serialize it
#pragma omp critical
  {
    pybind11::array_t<double> X(data.size(), data.data());
    auto res = eval(X);
    size_t i = 0;
    for (auto item : res) {
      d[i++] = item.cast<double>();
      if (i > 2)
        break;
    }
  }
  return d;
}


Tensor PythonTensorFunc::evaluate (const Vec3& x) const
{
  const Vec4* Xt = dynamic_cast<const Vec4*>(&x);
  std::vector<double> data(x.ptr(), x.ptr()+3);
  if (Xt)
    data.push_back(Xt->t);

  Tensor d(3);

  // python code cannot be called on multiple threads - at least not easily.
  // for now serialize it
#pragma omp critical
  {
    pybind11::array_t<double> X(data.size(), data.data());
    auto res = eval(X);
    size_t i = 0;
    size_t nsd;
    size_t cmps = ncmp;
    if (ncmp == 0) {
      pybind11::array_t<double> f = res.cast<pybind11::array_t<double>>();
      cmps = f.size();
      nsd = sqrt(cmps);
    } else
      nsd = sqrt(ncmp);
    d = Tensor(nsd); // funkyness needed to reset dimensionality
    for (auto item : res) {
      d(1 + i / nsd, 1 + (i % nsd)) = item.cast<double>();
      ++i;
      if (i >= cmps)
        break;
    }
  }
  return d;
}


SymmTensor PythonSTensorFunc::evaluate (const Vec3& x) const
{
  const Vec4* Xt = dynamic_cast<const Vec4*>(&x);
  std::vector<double> data(x.ptr(), x.ptr()+3);
  if (Xt)
    data.push_back(Xt->t);

  auto d = this->Result;

  // python code cannot be called on multiple threads - at least not easily.
  // for now serialize it
#pragma omp critical
  {
    pybind11::array_t<double> X(data.size(), data.data());
    auto res = eval(X);
    size_t i = 0;
    std::vector<Real>& svec = d;
    for (auto item : res) {
      svec[i++] = item.cast<double>();
      if (i >= svec.size())
        break;
    }
  }
  return d;
}


#endif
