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
                                    bool file)
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
  return Vec3((*f[0])(X), (*f[1])(X), f[2] ? (*f[2])(X) : 0.0);
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
                                          bool file)
{
  f.resize(input.size());
  for (size_t i = 0; i < input.size(); ++i)
    f[i] = std::make_unique<ChebyshevFunc>(input[i], file);
  ncmp = input.size();
}


Tensor ChebyshevTensorFunc::evaluate (const Vec3& X) const
{
  const size_t nsd = sqrt(ncmp);
  Tensor result(nsd);
  size_t c = 0;
  for (size_t j = 1; j <= nsd; ++j)
    for (size_t i = 1; i <= nsd; ++i, ++c)
      result(i,j) = (*f[c])(X);

  return result;
}
