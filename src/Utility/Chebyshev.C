// $Id$
//==============================================================================
//!
//! \file Chebyshev.C
//!
//! \date Jul 13 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Evaluation of Chebyshev polynomials.
//!
//==============================================================================

#include "Chebyshev.h"
#include "MatVec.h"

#include <fstream>
#include <sstream>
#include <numeric>


Real Chebyshev::evalPol1 (int polnum, Real xi)
{
  if (polnum <= 0)
    return 1.0;
  if (polnum == 1)
    return xi;

  return 2.0*xi*evalPol1(polnum-1, xi) - evalPol1(polnum-2, xi);
}


Real Chebyshev::evalPol2 (int polnum, Real xi)
{
  if (polnum <= 0)
    return 1.0;
  if (polnum == 1)
    return 2.0*xi;

  return 2.0*xi*evalPol2(polnum-1, xi) - evalPol2(polnum-2, xi);
}


Real Chebyshev::evalDer1 (int polnum, Real xi)
{
  if (polnum <= 0)
    return 0.0;

  return polnum * evalPol2(polnum-1, xi);
}


Real Chebyshev::evalDer2 (int polnum, Real xi)
{
  return ((polnum + 1)*evalPol1(polnum + 1, xi) - xi*evalPol2(polnum, xi)) / (xi*xi - 1.0);
}


Real Chebyshev::eval2Der1 (int polnum, Real xi)
{
  if (polnum < 2)
    return 0.0;

  if (std::abs(xi-1.0) < 1e-6)
    return (pow(polnum,4) - pow(polnum, 2)) / 3.0;
  else if (std::abs(1.0+xi) < 1e-6)
    return pow(-1, polnum) * (pow(polnum,4) - pow(polnum, 2)) / 3.0;

  return polnum * ((polnum+1)*evalPol1(polnum, xi) - evalPol2(polnum, xi)) / (xi*xi - 1.0);
}


ChebyshevFunc::ChebyshevFunc (const std::string& input, bool file)
{
  if (file) {
    std::ifstream in(input);
    if (!in.good()) {
      n[0] = n[1] = n[2] = 0;
      return;
    }
    read(in);
  } else {
    std::stringstream in;
    in.str(input);
    read(in);
  }
}


void ChebyshevFunc::read (std::istream& in)
{
  for (size_t i = 0; i < 3; ++i) {
    in >> n[i] >> domain[i][0] >> domain[i][1];
    if (n[i] == 0)
      n[i] = 1;
  }
  coefs.resize(n[0]*n[1]*n[2]);
  for (double& coef : coefs)
    in >> coef;
}


Real ChebyshevFunc::evaluate (const Vec3& X) const
{
  const Func eval{Chebyshev::evalPol1, 1.0};
  return this->evaluateTP(X, {eval, eval, eval});
}


Real ChebyshevFunc::evaluateTP (const Vec3& X,
                                const std::array<Func,3>& funcs) const
{
  auto eval = [&funcs,&X,this](int c)
  {
    Vector V(n[c]);
    const double coord = -1.0 + 2.0 * (X[c] - domain[c][0]) /
                                      (domain[c][1] - domain[c][0]);
    for (int i = 0; i < n[c]; ++i)
      V[i] = funcs[c].f(i, coord) * funcs[c].w;
    return V;
  };

  Vector TX = eval(0);
  if (n[1] == 1 && n[2] == 1)
    return std::inner_product(coefs.begin(), coefs.end(), TX.begin(), 0.0);

  Matrix T;
  Vector TY = eval(1);
  T.outer_product(TX, TY);
  if (n[2] == 1)
    return std::inner_product(coefs.begin(), coefs.end(), T.begin(), 0.0);

  Vector TZ = eval(2);
  Matrix T2;
  T2.outer_product(T, TZ);
  return std::inner_product(coefs.begin(), coefs.end(), T2.begin(), 0.0);
}


Real ChebyshevFunc::deriv (const Vec3& X, int c) const
{
  if (n[c-1] == 1)
    return 0.0;

  const Func der{Chebyshev::evalDer1, 2.0 / (domain[c-1][1]-domain[c-1][0])};
  const Func eval{Chebyshev::evalPol1, 1.0};

  return this->evaluateTP(X, {c == 1 ? der : eval,
                              c == 2 ? der : eval,
                              c == 3 ? der : eval});
}


Real ChebyshevFunc::dderiv (const Vec3& X, int c1, int c2) const
{
  if (n[c1-1] == 1 || n[c2-1] == 1)
    return 0.0;

  const Func eval{Chebyshev::evalPol1, 1.0};

  if (c1 == c2) {
    const Func der2{Chebyshev::eval2Der1, 4.0 / pow(domain[c1-1][1]-domain[c1-1][0], 2.0)};
    return this->evaluateTP(X, {c1 == 1 ? der2 : eval,
                                c1 == 2 ? der2 : eval,
                                c1 == 3 ? der2 : eval});
  }

  const Func der1{Chebyshev::evalDer1, 2.0 / (domain[c1-1][1]-domain[c1-1][0])};
  const Func der2{Chebyshev::evalDer1, 2.0 / (domain[c2-1][1]-domain[c2-1][0])};
  return this->evaluateTP(X, {c1 == 1 ? der1 : (c2 == 1 ? der2 : eval),
                              c1 == 2 ? der1 : (c2 == 2 ? der2 : eval),
                              c1 == 3 ? der1 : (c2 == 3 ? der2 : eval)});
}


ChebyshevVecFunc::ChebyshevVecFunc (const std::vector<std::string>& input,
                                    bool file, bool second)
  : secondDer(second)
{
  f[0] = std::make_unique<ChebyshevFunc>(input[0], file);
  if (input.size() > 1)
    f[1] = std::make_unique<ChebyshevFunc>(input[1], file);
  if (input.size() > 2)
    f[2] = std::make_unique<ChebyshevFunc>(input[2], file);
  if (f[0]->getSize()[2] == 1)
    ncmp = 2;
}


Vec3 ChebyshevVecFunc::evaluate (const Vec3& X) const
{
  const Vec4& X4 = static_cast<const Vec4&>(X);
  Vec3 res;

  // Multiple-components - no derivatives
  if (f[1]) {
    res[0] = (*f[0])(X);
    res[1] = (*f[1])(X);
    if (f[2])
      res[2] = (*f[2])(X);
    return res;
  }

  const std::array<int,3>& n = f[0]->getSize();
  const std::vector<Real>& coefs = f[0]->getCoefs();

  Vector TX(n[0]), dTX(n[0]);
  Matrix T;
  for (int i = 0; i < n[0]; ++i) {
    TX[i] = Chebyshev::evalPol1(i, (-1.0+2.0*X4.u[0]));
    if (secondDer)
      dTX[i] = 4.0*Chebyshev::eval2Der1(i, (-1.0+2.0*X4.u[0])); // 4.0 due to dxi/du twice
    else
      dTX[i] = 2.0*Chebyshev::evalDer1(i, (-1.0+2.0*X4.u[0])); // 2.0 due to dxi/du
  }
  Vector TY(n[1]), dTY(n[1]);
  for (int j = 0; j < n[1]; ++j) {
    TY[j] = Chebyshev::evalPol1(j, (-1.0+2.0*X4.u[1]));
    if (secondDer)
      dTY[j] = 4.0*Chebyshev::eval2Der1(j, (-1.0+2.0*X4.u[1])); // 4.0 due to dxi/du twice
    else
      dTY[j] = 2.0*Chebyshev::evalDer1(j, (-1.0+2.0*X4.u[1])); // 2.0 due to dxi/du
  }
  if (n[2] == 1) {
    T.outer_product(dTX, TY);
    res[0] = std::inner_product(coefs.begin(), coefs.end(), T.begin(), 0.0);
    T.outer_product(TX, dTY);
    res[1] = std::inner_product(coefs.begin(), coefs.end(), T.begin(), 0.0);
  } else {
    Vector TZ(n[2]), dTZ(n[2]);
    for (int k = 0; k < n[2]; ++k) {
      TZ[k] = Chebyshev::evalPol1(k, (-1.0 + 2.0*X4.u[2]));
      if (secondDer)
        dTZ[k] = Chebyshev::eval2Der1(k, (-1.0 + 2.0*X4.u[2])); // 4.0 due to dxi/du twice
      else
        dTZ[k] = Chebyshev::evalDer1(k, (-1.0 + 2.0*X4.u[2])); // 2.0 due to dxi/du
    }
    Matrix T2;
    T.outer_product(dTX, TY);
    T2.outer_product(T, TZ);
    res[0] = std::inner_product(coefs.begin(), coefs.end(), T2.begin(), 0.0);
    T.outer_product(TX, dTY);
    T2.outer_product(T, TZ);
    res[1] = std::inner_product(coefs.begin(), coefs.end(), T2.begin(), 0.0);
    T.outer_product(TX, TY);
    T2.outer_product(T, dTZ);
    res[2] = std::inner_product(coefs.begin(), coefs.end(), T2.begin(), 0.0);
  }

  return res;
}


std::vector<Real> ChebyshevVecFunc::evalGradient (const Vec3& X) const
{
  std::vector<Real> result;
  result.reserve(ncmp*ncmp);
  for (size_t d = 1; d <= ncmp; ++d)
    for (size_t c = 0; c < ncmp; ++c)
      result.push_back(f[c]->deriv(X, d));

  return result;
}


std::vector<Real> ChebyshevVecFunc::evalHessian (const Vec3& X) const
{
  std::vector<Real> result;
  result.reserve(ncmp*ncmp*ncmp);
  for (size_t d2 = 1; d2 <= ncmp; ++d2)
    for (size_t d1 = 1; d1 <= ncmp; ++d1)
      for (size_t c = 0; c < ncmp; ++c)
        result.push_back(f[c]->dderiv(X, d1, d2));

  return result;
}


ChebyshevTensorFunc::ChebyshevTensorFunc (const std::vector<std::string>& input,
                                          bool file, bool second)
{
  if (input.size() < 4) {
    for (size_t i = 0; i < input.size(); ++i)
      f[i].reset(new ChebyshevVecFunc({input[i]}, file, second));
  } else {
    if (input.size() == 4) {
      f[0].reset(new ChebyshevVecFunc({input[0], input[1]}, file));
      f[1].reset(new ChebyshevVecFunc({input[2], input[3]}, file));
    } else {
      f[0].reset(new ChebyshevVecFunc({input[0], input[1], input[2]}, file));
      f[1].reset(new ChebyshevVecFunc({input[3], input[4], input[5]}, file));
      f[2].reset(new ChebyshevVecFunc({input[6], input[7], input[8]}, file));
    }
  }
  ncmp = f[0]->dim();
}


Tensor ChebyshevTensorFunc::evaluate (const Vec3& X) const
{
  Tensor result(ncmp);
  Vec3 r1 = (*f[0])(X);
  Vec3 r2, r3;
  if (f[1])
    r2 = (*f[1])(X);
  if (f[2])
    r3 = (*f[2])(X);
  for (size_t i = 0; i < ncmp; ++i)
    result(1, i+1) = r1[i];
  if (f[1])
    for (size_t i = 0; i < ncmp; ++i)
      result(2, 1+i) = r2[i];
  if (f[2])
    for (size_t i = 0; i < ncmp; ++i)
      result(3, 1+i) = r3[i];

  return result;
}
