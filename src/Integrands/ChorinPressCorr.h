//==============================================================================
//!
//! \file ChorinPressCorr.h
//!
//! \date Sep 30 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Integrand implementations for pressure correction in Chorin's method
//!
//==============================================================================

#ifndef _CHORIN_PRESS_CORR_H
#define _CHORIN_PRESS_CORR_H

#include "IntegrandBase.h"
#include "SIMenums.h"
#include "ElmMats.h"
#include "Vec3.h"
#include "TimeDomain.h"
#include <map>

class LocalSystem;
class VTF;

/*!
  \brief Class representing the integrand of the pressure correction step
  in Chorin's method.
  \details This class implements the pressure correction step of Chorin's
  method using NURBS based FEM with equal order elements for velocity and
  pressure. The resulting system is a Poisson equation for the pressure
  increment.
*/

class ChorinPressCorr : public Integrand
{
 public:
  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  ChorinPressCorr(unsigned short int n = 3, double B0 = 1.0,
		  bool incPress = false, bool mixed = false);
  //! \brief The destructor frees the dynamically allocated data objects.
  virtual ~ChorinPressCorr();

  //! \brief Defines the traction field to use in Neumann boundary conditions.
  void setTraction(VecFunc* tf) { tracFld = tf; }

  //! \brief Clears the integration point traction values.
  void clearTracVal() { tracVal.clear(); }

  //! \brief Defines fluid properties for current volume patch.
  //! \param[in] Density   Mass density
  //! \param[in] Viscosity Dynamic viscosity
  void setFluidProperties(double Density, double Viscosity)
  { rho = Density; mu = Viscosity; }

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  bool initElement(const std::vector<int>& MNPC);

  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC1 Nodal point correspondance for the velocity
  //! \param[in] MNPC2 Nodal point correspondance for the pressure
  //! \param[in] n1 Number of nodes in basis 1 on this patch
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
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral*& elmInt, 
	       const TimeDomain& time, double detJW,
	       const Vector& N, const Matrix& dNdX,
	       const Vec3& X) const;

  //! \brief Evaluates the mixed field problem integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N1 Basis function values, field velocity
  //! \param[in] N2 Basis function values, field pressure
  //! \param[in] dN1dX Basis function gradients, velocity
  //! \param[in] dN2dX Basis function gradients, pressure
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral*& elmInt, 
	       const TimeDomain& time, double detJW,
	       const Vector& N1, const Vector& N2,
	       const Matrix& dN1dX, const Matrix& dN2dX,
	       const Vec3& X) const;
  
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N Basis function values
  //! \param[in] dNdX Basis function gradients
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  bool evalBou(LocalIntegral*& elmInt, 
	       const TimeDomain& time, double detJW,
	       const Vector& N, const Matrix& dNdX,
	       const Vec3& X, const Vec3& normal) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] detJW Jacobian determinant times integration point weight
  //! \param[in] N2 Basis function values, pressure
  //! \param[in] dN1dX Basis function gradients, pressure
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  //!
  //! \details The boundary integral is the same as that of the parent class.
  //! It does not depend on the pressure and volumetric change fields.
  //! Thus, this call is forwarded to the single-field parent class method.
  virtual bool evalBou(LocalIntegral*& elmInt, 
                       const TimeDomain& time, double detJW,
                       const Vector&, const Vector& N2,
                       const Matrix&, const Matrix& dN2dX,
                       const Vec3& X, const Vec3& normal) const
  {
    return this->evalBou(elmInt,time,detJW,N2,dN2dX,X,normal);
  }

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  bool evalSol(Vector& s,
	       const Vector& N, const Matrix& dNdX,
	       const Vec3& X, const std::vector<int>& MNPC) const;

  //! \brief Evaluates the secondary solution at a result point (mixed problem).
  //! \param[out] s The solution field values at current point
  //! \param[in] N2 Basis function values at current point, pressure
  //! \param[in] dN2dX Basis function gradients at current point, pressure
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC2 Nodal point correspondance for the pressure
  virtual bool evalSol(Vector& s,
                       const Vector&, const Vector& N2,
                       const Matrix&, const Matrix& dN2dX, const Vec3& X,
                       const std::vector<int>& MNPC2) const
  {
    return this->ChorinPressCorr::evalSol(s,N2,dN2dX,X,MNPC2);
  }

  //! \brief Evaluates the analytical solution at an integration point.
  //! \param[out] s The analytical secondary solution values at current point
  //! \param[in] asol The analytical solution field
  //! \param[in] X Cartesian coordinates of current point
  bool evalSolScal(Vector& s, const VecFunc& asol, const Vec3& X) const;

  //! \brief Returns whether there are any traction values to write to VTF.
  bool hasTractionValues() const { return !tracVal.empty(); }

  //! \brief Returns whether an incremental pressure formulation is used.
  bool incrementalPressure() const { return incPressure; }

  //! \brief If a mixed FE formulation is used
  bool mixedFormulation() const { return mixedFEM; }

  //! \brief Accesses the velocity solution vectors of current patch
  Vector& getVelocity(int n = 0) { return usol[n]; }

  //! \brief Accesses the pressure solution vectors of current patch
  Vector& getPressure(int n = 0) { return getSolution(n); }

  //! \brief If the element matrix should be assembled
  void assembleMatrix(bool newMat) { myMats->rhsOnly = !newMat; }

  //! \brief Get number of velocity solutions
  size_t getNoVelocities() const { return usol.size(); }

  //! \brief Get number of pressure solutions
  size_t getNoPressures() const { return this->getNoSolutions(); }

 protected:
  //! \brief Utility used by the virtual \a evalInt and \a evalBou methods.
  //! \param elmInt Pointer to the integrated element quantities
  bool getIntegralResult(LocalIntegral*& elmInt) const;

  // Physical properties (constant)
  double mu;   //!< Dynamic viscosity
  double rho;  //!< Fluid density

  unsigned short int nsd;   //!< Number of space dimensions (1, 2, or 3)

  bool   incPressure;       //!< Incremental pressure flag
  bool   mixedFEM;          //!< Mixed FE formulation
  
  double Beta0;             //!< Time integration parameter

  // Finite element quantities
  Matrix* eM;               //!< Pointer to element matrix
  Vector* eS;               //!< Pointer to element rhs vector
  std::vector<Vector*> eVs; //!< Pointer to element velocity vectors
  std::vector<Vector*> ePs; //!< Pointer to element pressure vectors

  std::vector<Vector> usol; //!< Velocity solution vectors for patch

  ElmMats* myMats;          //!< Local element matrices
  
  VecFunc* tracFld;                    //!< Pointer to boundary traction field
  mutable std::map<Vec3,Vec3> tracVal; //!< Traction field point values
};

#endif
