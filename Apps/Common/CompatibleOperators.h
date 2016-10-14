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
    //! \param[in] scale Scaling factor for contribution
    static void Advection(std::vector<Matrix>& EM, const FiniteElement& fe,
                          const Vec3& AC, double scale=1.0);

    //! \brief Compute a gradient term.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    static void Gradient(std::vector<Matrix>& EM,
                         const FiniteElement& fe, double scale=1.0);

    //! \brief Compute a gradient term.
    //! \param[out] EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    static void Gradient(Vectors& EV, const FiniteElement& fe,
                         double scale=1.0);

    //! \brief Compute a laplacian.
    //! \param[out] EM The element matrix to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] stress Whether to add extra stress formulation terms
    //! \details Only the upper blocks are added with the stress formulation
    static void Laplacian(std::vector<Matrix>& EM, const FiniteElement& fe,
                          double scale=1.0, bool stress=false);

    //! \brief Compute a mass term.
    //! \param[out] EM The element matrices to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] scale Scaling factor for contribution
    static void Mass(std::vector<Matrix>& EM,
                     const FiniteElement& fe, double scale=1.0);

    //! \brief Compute a vector-source term.
    //! \param[out] EV The element vectors to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] f Vector with contributions
    //! \param[in] scale Scaling factor for contribution
    static void Source(Vectors& EV, const FiniteElement& fe,
                       const Vec3& f, double scale=1.0);
  };

  //! \brief Common weak residual operators using div-compatible discretizations.
  class Residual
  {
    public:
    //! \brief Compute a laplacian term in a residual vector.
    //! \param[out] EV The element vector to add contribution to
    //! \param[in] fe The finite element to evaluate for
    //! \param[in] dUdX Current solution gradient
    //! \param[in] scale Scaling factor for contribution
    //! \param[in] stress Whether to add extra stress formulation terms
    static void Laplacian(Vectors& EV, const FiniteElement& fe,
                          const Tensor& dUdX, double scale=1.0,
                          bool stress=false);
  };
};

#endif
