// $Id$
//==============================================================================
//!
//! \file NewmarkMats.C
//!
//! \date Jul 4 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of the element matrices for a dynamic FEM problem.
//! \details See eq. (6.50) and (6.52) in Cottrell et. al. (2009).
//!
//==============================================================================

#include "NewmarkMats.h"


NewmarkMats::NewmarkMats (double a1, double a2, double b, double c,
                          bool generalizedAlpha) : isPredictor(true), h(0.0)
{
  alpha1 = a1;
  alpha2 = a2;

  if (generalizedAlpha)
  {
    alpha_m = fabs(b);
    alpha_f = c;
    double alpha = alpha_f - alpha_m;
    beta = 0.25*(1.0-alpha)*(1.0-alpha);
    gamma = 0.5 - alpha;
  }
  else
  {
    alpha_m = alpha_f = 1.0;
    beta  = fabs(b);
    gamma = c;
  }

  slvDisp = b < 0.0; // Displacement increments are used as primary unknowns
}


const Matrix& NewmarkMats::getNewtonMatrix () const
{
  if (A.empty())
    return this->ElmMats::getNewtonMatrix();

  Matrix& N = const_cast<Matrix&>(A.front());

  N = A[1];
  N.multiply(alpha_m + alpha_f*alpha1*gamma*h);
  if (A.size() > 2)
    N.add(A[2],alpha_f*(alpha2*gamma + beta*h)*h);
  if (A.size() > 3)
    N.add(A[3],alpha_f*gamma*h);
  if (slvDisp)
    N.multiply(1.0/(beta*h*h));
#if SP_DEBUG > 2
  std::cout <<"\nElement mass matrix"<< A[1];
  if (A.size() > 2)
    std::cout <<"Element stiffness matrix"<< A[2];
  if (A.size() > 3)
    std::cout <<"Element damping matrix"<< A[3];
  std::cout <<"Resulting Newton matrix"<< A[0];
#endif

  return A.front();
}


const Vector& NewmarkMats::getRHSVector () const
{
  if (b.empty())
    return this->ElmMats::getRHSVector();

  Vector& dF = const_cast<Vector&>(b.front());

  int ia = vec.size() - 1; // index to element acceleration vector (a)
  int iv = vec.size() - 2; // index to element velocity vector (v)

  bool haveMass = A.size() > 1 && !A[1].empty();
  bool haveStif = A.size() > 2 && !A[2].empty();
  bool haveDamp = A.size() > 3 && !A[3].empty();

#if SP_DEBUG > 2
  if (!vec.empty()) std::cout <<"\n";
  int iu = vec.size() - 3; // index to element displacement vector (u)
  if (iu >= 0) std::cout <<"u:"<< vec[iu];
  if (iv >= 0) std::cout <<"v:"<< vec[iv];
  if (ia >= 0) std::cout <<"a:"<< vec[ia];
  std::cout <<"\nf_ext - f_s"<< dF;
  if (haveMass && ia >= 0)
    std::cout <<"f_i = M*a"<< Vector(A[1]*vec[ia]);
  if (b.size() > 1)
    std::cout <<"f_d"<< b[1];
  else if (haveDamp && iv >= 0)
    std::cout <<"f_d = C*v"<< Vector(A[3]*vec[iv]);
  if (alpha1 > 0.0 && haveMass && iv >= 0)
    std::cout <<"f_d1/alpha1 = M*v (alpha1="<< alpha1 <<")"
              << Vector(A[1]*vec[iv]);
  if (alpha2 > 0.0 && haveStif && iv >= 0)
    std::cout <<"f_d2/alpha2 = K*v (alpha2="<< alpha2 <<")"
              << Vector(A[2]*vec[iv]);
#endif

  if (haveMass && ia >= 0)
    dF.add(A[1]*vec[ia],-1.0);    // dF = Fext - M*a

  if (b.size() > 1)
    dF.add(b[1],-1.0);            // dF -= Fd
  else if (haveDamp && iv >= 0)
    dF.add(A[3]*vec[iv],-1.0);    // dF -= C*v

  if (alpha1 > 0.0 && haveMass && iv >= 0)
    dF.add(A[1]*vec[iv],-alpha1); // dF -= alpha1*M*v

  if (alpha2 > 0.0 && haveStif && iv >= 0)
    dF.add(A[2]*vec[iv],-alpha2); // dF -= alpha2*K*v

#if SP_DEBUG > 2
  this->printVec(std::cout);
#endif

  return b.front();
}
