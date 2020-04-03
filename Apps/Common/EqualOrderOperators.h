//==============================================================================
//!
//! \file EqualOrderOperators.h
//!
//! \date Jul 22 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various equal-ordered discrete operators.
//!
//==============================================================================

#ifndef EQUAL_ORDER_OPERATORS_H
#define EQUAL_ORDER_OPERATORS_H

class Vec3;

#include "FiniteElement.h"
#include "MatVec.h"


namespace WeakOperators
{
  //! \brief Enum for the form of the convection term
  enum ConvectionForm
  {
    CONVECTIVE = 0,   //!<  u_i u_j,i v_j
    CONSERVATIVE = 1, //!< -u_i u_j v_i,j
    SKEWSYMMETRIC = 2 //!< (u_i u_j,i v_j - u_i u_j v_i,j)/2
  };
}


/*! \brief Common discrete operators using equal-ordered discretizations.
 */

class EqualOrderOperators
{
public:
  //! \brief Common weak operators using equal-ordered discretizations.
  class Weak {
  public:
    //! \brief Compute an advection term.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] AC Advecting field
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] form Formulation to use for advection
    //! \param[in] basis Basis to use
    static void Advection(Matrix& EM, const FiniteElement& fe,
                          const Vec3& AC, double scale=1.0,
                          WeakOperators::ConvectionForm form = WeakOperators::CONVECTIVE,
                          int basis=1);

    //! \brief Compute a (nonlinear) convection term.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] U Advecting field
    //! \param[in] dUdX Gradient of advected field
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] form Which form of the convective term to use
    //! \param[in] basis Basis to use
    static void Convection(Matrix& EM, const FiniteElement& fe,
                           const Vec3& U, const Tensor& dUdX, double scale,
                           WeakOperators::ConvectionForm form=WeakOperators::CONVECTIVE,
                           int basis=1);

    //! \brief Compute a divergence term.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] basis Basis for field
    //! \param[in] tbasis Test function basis
    static void Divergence(Matrix& EM, const FiniteElement& fe,
                           double scale=1.0, int basis=1, int tbasis=1);

    //! \brief Compute a divergence term.
    //! \param[out] EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] D Divergence of field
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] basis Test function basis
    static void Divergence(Vector& EV, const FiniteElement& fe,
                           const Vec3& D, double scale=1.0, int basis=1);

    //! \brief Compute a gradient term for a (potentially mixed) vector/scalar field.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] basis Basis for field
    //! \param[in] tbasis Test function basis
    static void Gradient(Matrix& EM, const FiniteElement& fe,
                         double scale=1.0, int basis=1, int tbasis=1);

    //! \brief Compute a gradient term.
    //! \param[out] EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] basis Basis to use
    static void Gradient(Vector& EV, const FiniteElement& fe,
                         double scale=1.0, int basis=1);

    //! \brief Compute a laplacian.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] stress Whether to add extra stress formulation terms
    //! \param[in] basis Basis to use
    static void Laplacian(Matrix& EM, const FiniteElement& fe,
                          double scale=1.0, bool stress=false, int basis=1);

    //! \brief Compute a heteregenous coefficient laplacian.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[out] K The coefficient matrix
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] basis Basis to use
    static void LaplacianCoeff(Matrix& EM, const Matrix& K, const FiniteElement& fe,
                               double scale=1.0, int basis=1);

    //! \brief Compute a mass term.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] basis Basis to use
    static void Mass(Matrix& EM, const FiniteElement& fe,
                     double scale=1.0, int basis=1);

    //! \brief Compute a source term.
    //! \param[out] EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] basis Basis to use
    //! \param[in] cmp Component to add (0 for all)
    static void Source(Vector& EV, const FiniteElement& fe,
                       double scale=1.0, int cmp=1, int basis=1);

    //! \brief Compute a vector-source term.
    //! \param[out] EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] f Vector with contributions
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] basis Basis to use
    static void Source(Vector& EV, const FiniteElement& fe,
                       const Vec3& f, double scale=1.0, int basis=1);
  };

  //! \brief Common weak residual operators using equal-ordered discretizations.
  class Residual {
  public:
    //! \brief Compute an advection term.
    //! \param[out] EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] AC Advecting field
    //! \param[in] g Advected field gradient
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] basis Basis to use
    static void Advection(Vector& EV, const FiniteElement& fe,
                          const Vec3& AC, const Tensor& g,
                          double scale = 1.0, int basis=1);

    //! \brief Compute a convection term in a residual vector.
    //! \param EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] U Advected field
    //! \param[in] dUdX Advected field gradient
    //! \param[in] UC Advecting field
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] form Which form of the convective term to use
    //! \param[in] basis Basis to use
    static void Convection(Vector& EV, const FiniteElement& fe,
                           const Vec3& U, const Tensor& dUdX, const Vec3& UC,
                           double scale,
                           WeakOperators::ConvectionForm form=WeakOperators::CONVECTIVE,
                           int basis=1);

    //! \brief Compute a divergence term in a residual vector.
    //! \param EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] dUdX Gradient of field
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] basis Basis to use
    static void Divergence(Vector& EV, const FiniteElement& fe,
                           const Tensor& dUdX, double scale=1.0,
                           size_t basis=1);

    //! \brief Compute a laplacian term in a residual vector.
    //! \param EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] dUdX Current solution gradient
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] basis Basis to use
    static void Laplacian(Vector& EV, const FiniteElement& fe,
                          const Vec3& dUdX, double scale=1.0,
                          int basis=1);

    //! \brief Compute a laplacian term in a residual vector.
    //! \param EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] dUdX Current solution gradient
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] stress Whether to add extra stress formulation terms
    //! \param[in] basis Basis to use
    template<class T>
    static void Laplacian(Vector& EV, const FiniteElement& fe,
                          const T& dUdX, double scale=1.0,
                          bool stress=false, int basis=1)
    {
      size_t cmp = EV.size() / fe.basis(basis).size();
      for (size_t i = 1;i <= fe.basis(basis).size();i++) {
        for (size_t k = 1;k <= cmp;k++) {
          double diff = 0.0;
          for (size_t l = 1;l <= fe.grad(basis).cols();l++)
            diff += dUdX(k,l)*fe.grad(basis)(i,l);
          diff *= scale*fe.detJxW;

          // Add residual to rhs of momentum equation
          EV((i-1)*cmp + k) += diff;
        }
      }

      // Use stress formulation
      if (stress) {
        for (size_t i = 1;i <= fe.basis(basis).size();i++)
          for (size_t k = 1;k <= cmp;k++)
            for (size_t l = 1;l <= cmp;l++)
              // Diffusion
              EV((i-1)*cmp + k) += scale*dUdX(l,k)*fe.grad(basis)(i,l)*fe.detJxW;
      }
    }
  };
};

#endif
