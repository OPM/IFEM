//==============================================================================
//!
//! \file ChorinVelPred.h
//!
//! \date Sep 29 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Integrand implementation for velocity prediction in Chorin's method
//!
//==============================================================================

#ifndef _CHORIN_VEL_PRED_H
#define _CHORIN_VEL_PRED_H

#include "Stokes.h"
#include "ElmMats.h"
#include "Vec3.h"
#include "TimeDomain.h"
#include "SIMenums.h"
#include <map>

class LocalSystem;
class VTF;


/*!
  \brief Class representing the integrand of the velocity prediction step
  in Chorin's method.
  \details This class implements the velocity prediction step of Chorin's
  method using NURBS based FEM with equal order elements for velocity and
  pressure. 
*/

class ChorinVelPred : public Stokes
{
 protected:
  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  //! \param[in] form The solution formulation to use (Laplace/Stress)
  //! \param[in] itg The integrandtype to use
  //! \param[in] incPress \e true if incremental pressure formulation
  ChorinVelPred(short int n, SIM::Formulation form, 
		int itg, bool incPress, bool mixed = false);

public:
  //! \brief The destructor frees dynamically allocated objects.
  virtual ~ChorinVelPred();

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  bool initElement(const std::vector<int>& MNPC);

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC1 Nodal point correspondance for the velocity
  //! \param[in] MNPC2 Nodal point correspondance for the pressure
  //! \param[in] n1 Number of nodes in velocity basis on this patch
  bool initElement(const std::vector<int>& MNPC1,
                   const std::vector<int>& MNPC2, size_t n1);

  //! \brief Initializes current element for boundary integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  bool initElementBou(const std::vector<int>& MNPC);
  //! \brief Initializes current element for numerical boundary integration.
  //! \param[in] MNPC1 Nodal point correspondance for the velocity
  //! \param[in] MNPC2 Nodal point correspondance for the pressure
  bool initElementBou(const std::vector<int>& MNPC1,
                      const std::vector<int>& MNPC2, size_t);

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  //!
  //! \details This interface is used when \a getIntegrandType returns 1.
  virtual bool evalInt(LocalIntegral*& elmInt,
		       const TimeDomain& time, double detJW,
		       const Vector& N, const Matrix& dNdX,
		       const Vec3& X) const 
  { return false; }
  
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] d2NdX2 Basis function second derivatives
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] h Characteristic element length
  //!
  //! \details This interface is used when \a getIntegrandType returns 2.
  virtual bool evalInt(LocalIntegral*& elmInt,
		       const TimeDomain& time, double detJW,
		       const Vector& N, const Matrix& dNdX,
		       const Matrix3D& d2NdX2, const Vec3& X,
		       double h = 0.0) const
  { return false; }

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral*& elmInt, const TimeDomain& time,
		       double detJW, const Vector& N, const Matrix& dNdX,
                       const Vec3& X, const Vec3& normal) const
  { return false; }

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N1 Basis function values, velocity
  //! \param[in] dN1dX Basis function gradients, velocity
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  //!
  //! \details The boundary integral is the same as that of the parent class.
  //! It does not depend on the pressure and volumetric change fields.
  //! Thus, this call is forwarded to the single-field parent class method.
  bool evalBou(LocalIntegral*& elmInt, 
               const TimeDomain& time, double detJW,
               const Vector& N1, const Vector&,
               const Matrix& dN1dX, const Matrix&,
               const Vec3& X, const Vec3& normal) const
  {
    return this->evalBou(elmInt,time,detJW,N1,dN1dX,X,normal);
  }

  //! \brief Evaluates the secondary solution at current integration point.
  //! \param[out] s The solution field values
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] MNPC Matrix of nodal point correspondance
  virtual bool evalSol(Vector& s,
                       const Vector& N, const Matrix& dNdX,
                       const Vec3& X, const std::vector<int>& MNPC) const;
  
  //! \brief Evaluates the secondary solution at a result point (mixed problem).
  //! \param[out] s The solution field values at current point
  //! \param[in] N1 Basis function values at current point, velocity
  //! \param[in] N2 Basis function values at current point, pressure
  //! \param[in] dN1dX Basis function gradients at current point, velocity
  //! \param[in] dN2dX Basis function gradients at current point, pressure
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC1 Nodal point correspondance for the velocity
  //! \param[in] MNPC2 Nodal point correspondance for the pressure
  bool evalSol(Vector& s,
               const Vector& N1, const Vector& N2,
               const Matrix& dN1dX, const Matrix& dN2dX, const Vec3& X,
               const std::vector<int>& MNPC1,
               const std::vector<int>& MNPC2) const;

  //! \brief Evaluates the analytical secondary solution at the given point.
  //! \param[out] s The solution field values
  //! \param[in] asol The analytical solution field
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalSol(Vector& s, const TensorFunc& asol, const Vec3& X) const;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer returned.
  //! \param[in] asol Pointer to analytical solution field (optional)
  virtual NormBase* getNormIntegrand(AnaSol* asol = 0) const;

  //! \brief Returns a pointer to an Integrand for boundary force evaluation.
  //! \note The Integrand is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer returned.
  //! \param[in] asol Pointer to analytical solution field (optional)
  virtual NormBase* getForceIntegrand(AnaSol* asol = 0) const;
  
  //! \brief Accesses the velocity solution vectors of current patch
  Vector& getVelocity(int n = 0) { return this->getSolution(n); }

  //! \brief Accesses the pressure solution vectors of the current patch
  Vector& getPressure(int n = 0) { return psol[n]; }

  //! \brief Returns a pointer to current element solution vector.
  const Vector* getElementVelocity() const { return eVs[0]; }

  //! \brief Returns a pointer to current element solution vector.
  const Vector* getElementPressure() const { return ePs[0]; }

  //! \brief If an incremental pressure formulation is used
  bool incrementalPressure () const { return incPressure; }

  //! \brief If a mixed FE formulation is used
  bool mixedFormulation() const { return mixedFEM; }

  //! \brief Get number of velocity solutions
  size_t getNoVelocities() const { return this->getNoSolutions(); }

  //! \brief Get number of pressure solutions
  size_t getNoPressures() const { return psol.size(); }

  //! \brief Calculates viscous part of stress tensor at current point.
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[out] eps Strain tensor at current point
  virtual bool strain(const Matrix& dNdX, Tensor& eps) const;

  //! \brief Calculates the (Cauchy) stress tensor at current point
  //! \param[in] N Basis functions at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[out] sigma Strain tensor at current point
  bool stress(const Vector& N, const Matrix& dNdX, Tensor& sigma) const;

  //! \brief Calculates the (Cauchy) stress tensor at current point
  //! \param[in] N1 Velocity basis functions at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[out] sigma Strain tensor at current point
  bool stress(const Vector& N1, const Vector& N2, 
              const Matrix& dN1dX, const Matrix& dN2dX,
              Tensor& sigma) const;

 protected:
  bool incPressure;           //!< Incremental pressure formulation
  bool mixedFEM;              //!< Mixed FE formulation

  // Finite element quantities
  std::vector<Vector*> ePs;   //!< Pointers to element pressure vectors 
  std::vector<Vector>  psol;  //!< Pressure solution vectors for this patch
};


/*!
  \brief Class representing the integrand of Stokes norms for Chorin method.
*/


class ChorinStokesNorm : public StokesNorm
{
 public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The Stokes problem to evaluate norms for
  //! \param[in] a The analytical velocity and pressure fields (optional)
 ChorinStokesNorm(Stokes& p, AnaSol* a) : StokesNorm(p,a) {}
  //! \brief Empty destructor.
  virtual ~ChorinStokesNorm() {}

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral*& elmInt, double detJW,
                       const Vector& N, const Matrix& dNdX,
                       const Vec3& X) const;

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N1 Basis function values, velocity
  //! \param[in] N2 Basis function values, pressure
  //! \param[in] dN1dX Basis function gradients, velocity
  //! \param[in] dN2dX Basis function gradients, pressure
  //! \param[in] X Cartesian coordinates of current integration point
  //!
  //! \details This interface is used for mixed formulations only.
  //! The default implementation forwards to the stationary version.
  //! Reimplement this method for time-dependent or non-linear problems.
  virtual bool evalInt(LocalIntegral*& elmInt,
                       const TimeDomain& time, double detJW,
                       const Vector& N1, const Vector& N2,
                       const Matrix& dN1dX, const Matrix& dN2dX,
                       const Vec3& X) const;

  //! \brief Returns the number of primary norm quantities.
  virtual size_t getNoFields() const { return anasol ? 4 : 2; }

  //! \brief Returns the number of secondary norm quantities.
  virtual size_t getNoSecFields() const { return anasol ? 4 : 1; }
};


/*!
  \brief Class representing for computing boundary force
*/

class ChorinStokesForce : public StokesForce
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The Stokes problem to evaluate norms for
  //! \param[in] a The analytical velocity and pressure fields (optional)
 ChorinStokesForce(ChorinVelPred& p, AnaSol* a = 0) : StokesForce(p,a) {}
  //! \brief Empty destructor.
  virtual ~ChorinStokesForce() {}

  //! \brief Initializes current element for boundary integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  virtual bool initElementBou(const std::vector<int>& MNPC);
  //! \brief Initializes current element for boundary integration (mixed).
  //! \param[in] MNPC1 Nodal point correspondance for the basis 1
  //! \param[in] MNPC2 Nodal point correspondance for the basis 2
  //! \param[in] n1 Number of nodes in basis 1 on this patch
  virtual bool initElementBou(const std::vector<int>& MNPC1,
                              const std::vector<int>& MNPC2, 
                              size_t n1);

    //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  //!
  //! \details The default implementation forwards to the stationary version.
  //! Reimplement this method for time-dependent or non-linear problems.
  virtual bool evalBou(LocalIntegral*& elmInt,
                       const TimeDomain& time, double detJW,
                       const Vector& N, const Matrix& dNdX,
                       const Vec3& X, const Vec3& normal) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  //! \param[in] N1 Basis function values, field 1
  //! \param[in] N2 Basis function values, field 2
  //! \param[in] dN1dX Basis function gradients, field 1
  //! \param[in] dN2dX Basis function gradients, field 2
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  //!
  //! \details This interface is used for mixed formulations.
  //! The default implementation forwards to the stationary version.
  //! Reimplement this method for time-dependent or non-linear problems.
  virtual bool evalBou(LocalIntegral*& elmInt,
                       const TimeDomain& time, double detJW,
                       const Vector& N1, const Vector& N2,
                       const Matrix& dN1dX, const Matrix& dN2dX,
                       const Vec3& X, const Vec3& normal) const
  { return false; }
};
 
#endif
