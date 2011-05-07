// $Id$
//==============================================================================
//!
//! \file AxSymmElasticity.h
//!
//! \date Apr 28 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for axial-symmetric elasticity problems.
//!
//==============================================================================

#ifndef _AX_SYMM_ELASTICITY_H
#define _AX_SYMM_ELASTICITY_H

#include "Elasticity.h"


/*!
  \brief Class representing the integrand of axial-symmetric elasticity problem.
  \details Most methods of this class are inherited form the base class.
  Only the \a evalInt, \a evalBou and \a evalSol methods, which are specific
  for axial-symmetric linear elasticity problems are implemented here.
*/

class AxSymmElasticity : public Elasticity
{
public:
  //! \brief Default constructor.
  AxSymmElasticity();
  //! \brief Empty destructor.
  virtual ~AxSymmElasticity() {}

  //! \brief Prints out the problem definition to the given output stream.
  virtual void print(std::ostream& os) const;

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cylindric coordinates of current integration point
  virtual bool evalInt(LocalIntegral*& elmInt, const FiniteElement& fe,
		       const Vec3& X) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cylindric coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral*& elmInt, const FiniteElement& fe,
		       const Vec3& X, const Vec3& normal) const;

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] X Cylindric coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s, const Vector& N, const Matrix& dNdX,
		       const Vec3& X, const std::vector<int>& MNPC) const;

  //! \brief Evaluates the finite element (FE) solution at an integration point.
  //! \param[out] s The FE stress values at current point
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] X Cylindric coordinates of current point
  virtual bool evalSol(Vector& s, const Vector& N, const Matrix& dNdX,
		       const Vec3& X) const;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  virtual size_t getNoFields(int fld = 2) const { return fld > 1 ? 5 : 2; }
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual const char* getField1Name(size_t i, const char* prefix = 0) const;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual const char* getField2Name(size_t i, const char* prefix = 0) const;

  //! \brief Returns \e true if this is an axial-symmetric problem.
  virtual bool isAxiSymmetric() const { return true; }

protected:
  //! \brief Calculates some kinematic quantities at current point.
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] r Radial coordinate of current point
  //! \param[out] eps Strain tensor at current point
  //!
  //! \details The strain displacement matrix \b B is established
  //! and stored in the mutable class member \a Bmat.
  bool kinematics(const Vector& N, const Matrix& dNdX,
		  double r, SymmTensor& eps) const;

private:
  // Work array declared as member to avoid frequent re-allocation
  // within the numerical integration loop (for reduced overhead)
  mutable Matrix CB; //!< Result of the matrix-matrix product C*B
};

#endif
