// $Id$
//==============================================================================
//!
//! \file MortarContact.h
//!
//! \date Jun 05 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Integrand implementations for Mortar-based contact analysis.
//!
//==============================================================================

#ifndef _MORTAR_CONTACT_H
#define _MORTAR_CONTACT_H

#include "GlobalIntegral.h"
#include "IntegrandBase.h"
#include "SparseMatrix.h"
#include "ElmMats.h"

class LocalIntegral;
class RigidBody;
class SAM;

typedef unsigned short int usint; //!< Convenience type declaration


/*!
  \brief Class for storage of system level Mortar matrices.
*/

class MortarMats : public GlobalIntegral
{
public:
  //! \brief The constructor initializes the Mortar matrices to proper size.
  //! \param[in] _sam Reference to the FE assembly management object to use
  //! \param[in] nMast Total number of master nodes
  //! \param[in] n Number of space dimensions
  MortarMats(const SAM& _sam, int nMast, usint n);
  //! \brief Empty destructor.
  virtual ~MortarMats() {}

  //! \brief Initializes the Mortar matrices to zero.
  virtual void initialize(bool);
  //! \brief Finalizes the Mortar matrix assembly.
  virtual bool finalize(bool);

  //! \brief Adds a set of element matrices into the global Mortar matrices.
  //! \param[in] elmObj Pointer to the element matrices to add into \a *this
  //! \param[in] elmId Global number of the element associated with \a *elmObj
  virtual bool assemble(const LocalIntegral* elmObj, int elmId);

  //! \brief Returns whether any nodes currently are in contact or not.
  bool haveContact() const { return contact; }
  //! \brief Returns the number of elements connected to a Lagrange multiplier.
  int getALconnectivity(size_t n) const { return nNoEl[n-1]; }
  //! \brief Returns the weighted area of the given node.
  double weightedArea(size_t n) const { return AA(n); }
  //! \brief Returns the weighted gap for the given node.
  double weightedGap(size_t n) const { return gNA(n); }
  //! \brief Returns an element of the auxiliary Mortar matrix.
  double phi(size_t i, size_t j) const { return phiA(i,j); }

private:
  const SAM&       sam;     //!< Data for FE assembly management
  SparseMatrix     phiA;    //!< Auxiliary constants
  Vector           gNA;     //!< Weighted nodal gaps (non-positive if contact)
  Vector           AA;      //!< Weighted nodal areas (zero if not contact)
  std::vector<int> nNoEl;   //!< Number of elements contributing to each node
  usint            nsd;     //!< Number of space dimensions
  bool             contact; //!< If \e true, at least one node is in contact
};


/*!
  \brief Class representing the integrand for the weighted gap matrices.
*/

class MortarContact : public IntegrandBase
{
public:
  //! \brief The contructor initializes the contact body of this integrand.
  //! \param[in] mst Pointer to the rigid body object to be in contact
  //! \param[in] gi Pointer to the global integrated quantity of this integrand
  //! \param[in] nsd Number of space dimensions
  MortarContact(RigidBody* mst, GlobalIntegral* gi, usint nsd);
  //! \brief The destructor deletes the global integral quantity.
  virtual ~MortarContact() { if (myInt) delete myInt; }

  //! \brief Returns that this integrand has no explicit interior contributions.
  virtual bool hasInteriorTerms() const { return false; }

  //! \brief Defines which FE quantities are needed by the integrand.
  virtual int getIntegrandType() const { return XO_ELEMENTS; }
  //! \brief Returns the number of boundary integration points.
  virtual int getBouIntegrationPoints(int nGP) const { return 2*nGP; }

  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t, bool) const;
  //! \brief Returns the system quantity to be integrated by \a *this.
  virtual GlobalIntegral& getGlobalInt(GlobalIntegral*) const { return *myInt; }

  //! \brief Initializes current element for boundary integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] elmInt Local integral for element
  virtual bool initElementBou(const std::vector<int>& MNPC,
                              LocalIntegral& elmInt);

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

protected:
  //! \brief Accumulates displacement residual and associated tangent stiffness.
  //! \param R Residual force vector
  //! \param Kt Tangent stiffness matrix
  //! \param[out] mNorm Master surface normal at current integration point
  //! \param[out] Nm Interpolation functions of the master surface DOFs
  //! \param[in] Ns Interpolation functions of the slave surface DOFs
  //! \param[in] X Updated Cartesian coordinates of current integration point
  //! \param[in] lambda Weighted nodal gap or Lagrange multiplier values
  //! \param[in] phi Auxiliary mortar matrix for current element
  //! \param[in] detJxW Jacobian determinant times integration weight
  //! \param[in] epsJxW Scaled Jacobian determinant times integration weight
  void accResAndTangent(Vector& R, Matrix& Kt,
		        Vec3& mNorm, Vector& Nm,
		        const Vector& Ns, const Vec3& X,
		        const Vector& lambda, const Matrix& phi,
		        double detJxW, double epsJxW = 0.0) const;

private:
  GlobalIntegral* myInt; //!< Pointer to resulting global integrated quantity

protected:
  RigidBody* master; //!< The rigid body to be in contact
};


/*!
  \brief Class representing the integrand of the Penalty-based Mortar contact.
*/

class MortarPenalty : public MortarContact
{
public:
  //! \brief The contructor initializes the contact body of this integrand.
  //! \param[in] mst Pointer to the rigid body object to be in contact
  //! \param[in] mats Pre-integrated Mortar matrices associated with \a mst.
  //! \param[in] nsd Number of space dimensions
  MortarPenalty(RigidBody* mst, const MortarMats& mats, usint nsd);
  //! \brief Empty destructor.
  virtual ~MortarPenalty() {}

  //! \brief Returns whether this integrand has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return mortar.haveContact(); }

  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t, bool) const;
  //! \brief Returns the system quantity to be integrated by \a *this.
  virtual GlobalIntegral& getGlobalInt(GlobalIntegral* eq) const { return *eq; }

  //! \brief Initializes current element for boundary integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] elmInt Local integral for element
  virtual bool initElementBou(const std::vector<int>& MNPC,
                              LocalIntegral& elmInt);

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

private:
  const MortarMats& mortar; //!< Pre-integrated system-level Mortar matrices
};


/*!
  \brief Class representing the integrand of the Augmented-Lagrange contact.
*/

class MortarAugmentedLag : public MortarContact
{
  //! \brief Class representing the element matrices of the AL formulation.
  class ALElmMats : public ElmMats
  {
  public:
    //! \brief Default constructor.
    ALElmMats(size_t nedof);
    //! \brief Empty destructor.
    virtual ~ALElmMats() {}
    //! \brief Returns the element-level Newton matrix.
    virtual const Matrix& getNewtonMatrix() const;
    //! \brief Returns the element-level right-hand-side vector
    //! associated with the Newton matrix.
    virtual const Vector& getRHSVector() const;

    std::vector<size_t> iLag; //!< Indices for the active Lagrange multipliers
  };

public:
  //! \brief The contructor initializes the contact body of this integrand.
  //! \param[in] mst Pointer to the rigid body object to be in contact
  //! \param[in] mats Pre-integrated Mortar matrices associated with \a mst.
  //! \param[in] nsd Number of space dimensions
  MortarAugmentedLag(RigidBody* mst, const MortarMats& mats, usint nsd);
  //! \brief Empty destructor.
  virtual ~MortarAugmentedLag() {}

  //! \brief Returns a local integral container for the given element.
  virtual LocalIntegral* getLocalIntegral(size_t, size_t, bool) const;
  //! \brief Returns the system quantity to be integrated by \a *this.
  virtual GlobalIntegral& getGlobalInt(GlobalIntegral* eq) const { return *eq; }

  //! \brief Initializes current element for boundary integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] elmInt Local integral for element
  virtual bool initElementBou(const std::vector<int>& MNPC,
                              LocalIntegral& elmInt);

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

private:
  const MortarMats& mortar; //!< Pre-integrated system-level Mortar matrices

  Matrix phiA; //!< Auxiliary constants
  Vector LNA;  //!< Lagrange multiplier values
};

#endif
