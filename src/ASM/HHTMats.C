// $Id$
//==============================================================================
//!
//! \file HHTMats.C
//!
//! \date Nov 13 2014
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of the element matrices for a dynamic FEM problem.
//!
//==============================================================================

#include "HHTMats.h"


HHTMats::HHTMats (double alpha, double a, double b, bool old) : NewmarkMats(a,b)
{
  beta = 0.25*(1.0-alpha)*(1.0-alpha);
  gamma = 0.5 - alpha;
  oldHHT = old;
}


const Matrix& HHTMats::getNewtonMatrix () const
{
  if (A.empty())
    return this->ElmMats::getNewtonMatrix();

  Matrix& N = const_cast<Matrix&>(A.front());
  const double alphaPlus1 = 1.5 - gamma; // = 1.0 + alpha
  const double gammaObh   = gamma/(beta*h);

  // Calculate the Newton matrix of the dynamic problem
  // N = (1+alpha)(1+alpha2*gobh)*K + (1/(beta*h^2)+(1+alpha)*alpha1*gobh)*M
  //   + (1+alpha)*gobh*C
  // where gobh = gamma/(beta*h)

  if (A[1].empty())
    N.fill(0.0);
  else
  {
    N = A[1];
    N.multiply(1.0/(beta*h*h) + alphaPlus1*alpha1*gammaObh);
  }
  if (A.size() > 2)
    N.add(A[2],alphaPlus1 + alphaPlus1*alpha2*gammaObh);
  if (A.size() > 3)
    N.add(A[3],alphaPlus1);
  if (A.size() > 4)
    N.add(A[4],alphaPlus1*gammaObh);

#if SP_DEBUG > 2
  std::cout <<"\nHHTMats::getNewtonMatrix";
  std::cout <<"\nElement mass matrix"<< A[1];
  if (A.size() > 3)
    std::cout <<"Material stiffness matrix"<< A[2]
              <<"Geometric stiffness matrix"<< A[3];
  else if (A.size() > 2)
    std::cout <<"Tangent stiffness matrix"<< A[2];
  if (A.size() > 4)
    std::cout <<"Element damping matrix"<< A[4];
  std::cout <<"scale(M) = "<< 1.0/(beta*h*h) + alphaPlus1*alpha1*gammaObh;
  if (A.size() > 4)
    std::cout <<"  scale(C) = "<< alphaPlus1*gammaObh;
  if (A.size() > 2)
    std::cout <<"  scale(K) = "<< alphaPlus1 + alphaPlus1*alpha2*gammaObh;
  std::cout <<"\nResulting Newton matrix"<< A[0];
#endif

  return A.front();
}


const Vector& HHTMats::getRHSVector () const
{
#if SP_DEBUG > 2
  std::cout <<"\nHHTMats::getRHSVector "
            << (isPredictor ? "(predictor)" : "(corrector)") << std::endl;
#endif

  if (b.empty() || vec.size() < 3)
    return this->ElmMats::getRHSVector();

  const Vector& vel = this->vel(); // current element velocity vector (v)
  const Vector& acc = this->acc(); // current acceleration vector (a)

  bool haveMass = A.size() > 1 && !A[1].empty();
  bool haveStif = A.size() > 2 && !A[2].empty();
  bool haveDamp = A.size() > 4 && !A[4].empty();

  if ((isPredictor || oldHHT) && b.size() > 1)
  {
    Vector& Fia = const_cast<Vector&>(b[1]);

    // Find the actual inertia force from the dynamic equilibrium equation
    // (Fia = Fext - Fs - Fd) and store it in the second RHS vector Fia=b[1]
    Fia = b.front(); // Fia = -Fs
    if (b.size() > 2)
      Fia.add(b[2]); // Fia = Fext - Fs
    if (alpha1 > 0.0 && haveMass)
      Fia.add(A[1]*vel,-alpha1); // Fia -= alpha1*M*v
    if (alpha2 > 0.0 && haveStif)
      Fia.add(A[2]*vel,-alpha2); // Fia -= alpha2*K*v
    if (haveDamp)
      Fia.add(A[4]*vel,-1.0);    // Fia -= C*v
#if SP_DEBUG > 2
    std::cout <<"Element inertia force vector"<< Fia;
#endif
  }

#if SP_DEBUG > 2
  std::cout <<"Element velocity vector"<< vel;
  std::cout <<"Element acceleration vector"<< acc;
  if (b.size() > 2)
    std::cout <<"S_ext"<< b[2] <<"-S_int"<< b.front();
  else
    std::cout <<"S_ext - S_int"<< b.front();
  if (b.size() > 3)
    std::cout <<"S_dmp"<< b[3];
  if (haveMass)
    std::cout <<"S_inert = M*a"<< Vector(A[1]*acc);
#endif

  // Calculate the right-hand-side force vector of the dynamic problem
  // Pred.: RHS = (1+alphaH)*{S_ext + [alpha1*M + alpha2*K + C]*V} + M*(A-a)
  // oldP.: RHS = (1+alphaH)*{S_ext - S_int + [alpha1*M + alpha2*K + C]*V} + M*a
  // Corr.: RHS = (1+alphaH)*{S_ext - S_int - [alpha1*M + alpha2*K + C]*V} - M*a
  // Note: The external load from the previous step is subtracted from
  // the predictor step force vector after the element assembly.
  Vector& RHS = const_cast<Vector&>(b.front());
  double alphaPlus1 = 1.5 - gamma;
  size_t pa = vec.size()-1; // index to predicted acceleration vector (A)
  size_t iv = 1;            // index to current velocity vector (v)
  if (haveMass)
  {
    if (oldHHT)
    {
      RHS *= alphaPlus1; // RHS = -(1+alphaH)*S_int
      RHS.add(A[1]*acc, isPredictor ? 1.0 : -1.0); // RHS (+/-)= M*a
    }
    else if (isPredictor && pa > 3*(vec.size()/3))
    {
      A[1].multiply(vec[pa]-acc,RHS); // RHS = M*(A-a)
      iv = pa-1; // index to predicted element velocity vector (V)
#if SP_DEBUG > 2
      std::cout <<"S_inert = M*(A-a)"<< RHS;
#endif
    }
    else
    {
      RHS *= alphaPlus1; // RHS = -(1+alphaH)*S_int
      RHS.add(A[1]*acc,-1.0); // RHS -= M*a
    }
  }

  if (oldHHT && b.size() > 2)
    RHS.add(b[2],alphaPlus1); // RHS += (1+alphaH)*S_ext

  if (!isPredictor) alphaPlus1 = -alphaPlus1;

  if (alpha1 > 0.0 && haveMass)
    RHS.add(A[1]*vec[iv],alphaPlus1*alpha1); // RHS -= (1+alphaH)*alpha1*M*v
  if (alpha2 > 0.0 && haveStif)
    RHS.add(A[2]*vec[iv],alphaPlus1*alpha2); // RHS -= (1+alphaH)*alpha2*K*v
  if (b.size() > 3)
    RHS.add(b[3]        ,alphaPlus1);        // RHS -= (1+alphaH)*Fd
  else if (haveDamp)
    RHS.add(A[4]*vec[iv],alphaPlus1);        // RHS -= (1+alphaH)*C*v

#if SP_DEBUG > 2
  this->printVec(std::cout);
#endif

  return b.front();
}
