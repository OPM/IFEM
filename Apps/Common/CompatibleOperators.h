//==============================================================================
//!
//! \file CompatibleOperators.h
//!
//! \date Oct 9 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various weak, div-compatible discrete operators.
//!
//==============================================================================

#ifndef COMPATIBLE_OPERATORS_H_
#define COMPATIBLE_OPERATORS_H_

class FiniteElement;
class Tensor;
class Vec3;

#include "MatVec.h"
#include "EqualOrderOperators.h"

#include <array>

/*! \brief Common operators using div-compatible discretizations.
 *  \details The operators use the block ordering used in the BlockElmMats class.
 */

class CompatibleOperators
{
public:
  //! \brief Common weak operators using div-compatible discretizations.
  class Weak {
  public:
    //! \brief Compute an advection term.
    //! \param[out] EM The element matrices to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] AC Advecting field
    //! \param[in] idx Matrix block indices
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] cnvForm Form of advection operator
    static void Advection(std::vector<Matrix>& EM, const FiniteElement& fe,
                          const Vec3& AC,
                          const std::array<std::array<int,3>,3>& idx,
                          double scale = 1.0,
                          WeakOperators::ConvectionForm cnvForm = WeakOperators::CONVECTIVE);

    //! \brief Compute a (nonlinear) convection term.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] U  Advecting field
    //! \param[in] dUdX Field gradient
    //! \param[in] idx Matrix block indices
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] form Which form of the convective form to use
    static void Convection(std::vector<Matrix>& EM, const FiniteElement& fe,
                           const Vec3& U, const Tensor& dUdX,
                           const std::array<std::array<int,3>,3>& idx,
                           double scale,
                           WeakOperators::ConvectionForm form = WeakOperators::CONVECTIVE);

    //! \brief Compute a gradient term.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] idx Matrix block indices
    //! \param[in] scale Scaling factor for contribution
    static void Gradient(std::vector<Matrix>& EM,
                         const FiniteElement& fe,
                         const std::array<int,3>& idx,
                         double scale = 1.0);

    //! \brief Compute an integration constraint.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[in] idx Matrix block indices
    //! \param[in] fe The finite element to evaluate for
    static void ItgConstraint(std::vector<Matrix>& EM,
                              const FiniteElement& fe,
                              const std::array<int,3>& idx);

    //! \brief Compute a laplacian.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] idx Matrix block indices
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] stress Whether to add extra stress formulation terms
    //! \details Only the upper blocks are added with the stress formulation
    static void Laplacian(std::vector<Matrix>& EM, const FiniteElement& fe,
                          const std::array<std::array<int,3>,3>& idx,
                          double scale = 1.0, bool stress = false);

    //! \brief Compute a mass term.
    //! \param[out] EM The element matrices to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] idx Matrix block indices
    //! \param[in] scale Scaling factor for contribution
    static void Mass(std::vector<Matrix>& EM,
                     const FiniteElement& fe,
                     const std::array<std::array<int,3>,3>& idx,
                     double scale = 1.0);

    //! \brief Compute a scalar source term.
    //! \param[out] EV The element vectors to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param idx Vector block indices
    //! \param[in] scale Scaling factor for contribution
    static void Source(Vectors& EV, const FiniteElement& fe,
                       const std::array<int,3>& idx, double scale = 1.0);

    //! \brief Compute a vector-source term.
    //! \param[out] EV The element vectors to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] f Vector with contributions
    //! \param idx Vector block indices
    //! \param[in] scale Scaling factor for contribution
    static void Source(Vectors& EV, const FiniteElement& fe,
                       const Vec3& f, const std::array<int,3>& idx,
                       double scale = 1.0);
  };

  //! \brief Common weak residual operators using div-compatible discretizations.
  class Residual
  {
  public:
    //! \brief Compute a convection term in a residual vector.
    //! \param EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] U Advected field
    //! \param[in] dUdX Advected field gradient
    //! \param[in] UC Advecting field
    //! \param idx Vector block indices
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] form Which form of the convective form to use
    static void Convection(Vectors& EV, const FiniteElement& fe,
                           const Vec3& U, const Tensor& dUdX,
                           const Vec3& UC,
                           const std::array<int,3>& idx,
                           double scale,
                           WeakOperators::ConvectionForm form = WeakOperators::CONVECTIVE);

    //! \brief Compute a gradient term.
    //! \param[out] EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param idx Vector block indices
    //! \param[in] scale Scaling factor for contribution
    static void Gradient(Vectors& EV, const FiniteElement& fe,
                         const std::array<int,3>& idx, double scale = 1.0);

    //! \brief Compute a laplacian term in a residual vector.
    //! \param[out] EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] dUdX Current solution gradient
    //! \param idx Vector block indices
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] stress Whether to add extra stress formulation terms
    static void Laplacian(Vectors& EV, const FiniteElement& fe,
                          const Tensor& dUdX,
                          const std::array<int,3>& idx,
                          double scale = 1.0,
                          bool stress = false);
  };
};

#endif
