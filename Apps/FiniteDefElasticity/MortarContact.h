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
  //! \param[in] nMast Total number of master nodes in the model
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

  //! \brief Returns the weighted area of the given node.
  double weightedArea(size_t n) const { return AA(n); }
  //! \brief Returns the weighted gap for the given node.
  double weightedGap(size_t n) const { return gNA(n); }
  //! \brief Returns an element of the auxiliary Mortar matrix.
  double phi(size_t m, size_t n) const { return phiA(m,n); }

  //! \brief Returns the total number of nodes for the Mortar matrices.
  size_t getNoNodes() const { return phiA.rows() / nsd; }
  //! \brief Returns the number of slave nodes for the Mortar matrices.
  size_t getNoSlaves() const { return phiA.cols(); }

  //! \brief Returns the SAM object associated with these Mortar matrices.
  const SAM& getSAM() const { return sam; }

private:
  const SAM&   sam;  //!< Data for FE assembly management
  SparseMatrix phiA; //!< Matrix of auxiliary constants
  Vector       gNA;  //!< Weighted nodal gaps
  Vector       AA;   //!< Weighted nodal areas
  usint        nsd;  //!< Number of space dimensions
};


/*!
  \brief Class representing the integrand for the weighted gap matrices.
  \details This class actually have two purposes. First of all it implements
  the integrand for computation of the global Mortar matrices (stored in the
  MortarMats class). Secondly, it is a base class for the integrand classes
  of the two formulations for calculation of residual and tangent contributions,
  MortarPenaty and MortarAugmentedLag, where the commonality of these two
  formulations is implemented in this class.
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

  //! \brief Assembles contributions to the tangent stiffness and residual.
  virtual bool assemble(SystemMatrix&, SystemVector&) const { return true; }

protected:
  //! \brief Accumulates displacement residual and associated tangent stiffness.
  //! \param R Residual force vector
  //! \param Kt Tangent stiffness matrix
  //! \param[out] mNorm Master surface normal at current integration point
  //! \param[out] Nm Interpolation functions of the master surface DOFs
  //! \param[in] Ns Interpolation functions of the slave surface DOFs
  //! \param[in] X Updated Cartesian coordinates of current integration point
  //! \param[in] lambda Weighted nodal gap or Lagrange multiplier values
  //! \param[in] epsJxW Scaled Jacobian determinant times integration weight
  void accResAndTangent(Vector& R, Matrix& Kt,
		        Vec3& mNorm, Vector& Nm,
		        const Vector& Ns, const Vec3& X,
		        const Vector& lambda, double epsJxW) const;

  //! \brief Assembles contributions to the tangent stiffness and residual.
  //! \param Ktan System tangent stiffness matrix
  //! \param Res System residual vector
  //! \param[in] mortar Mortar matrices giving the tangent contributions
  //!
  //! \details This method assembles direct nodal contributions to the
  //! system matrices emanating from the Mortar method formulation,
  //! which cannot be assembled through the normal element assembly loop.
  bool assResAndTangent(SystemMatrix& Ktan, SystemVector& Res,
			const MortarMats& mortar) const;

private:
  GlobalIntegral* myInt; //!< Pointer to resulting global integrated quantity

protected:
  RigidBody*        master;      //!< The rigid body to be in contact
  std::vector<bool> activeSlave; //!< Array of slave node status flags
};


/*!
  \brief Class representing the integrand of the Penalty-based contact.
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

  //! \brief Initializes the integrand for a new integration loop.
  //! \param[in] prm Nonlinear solution algorithm parameters
  virtual void initIntegration(const TimeDomain& prm, const Vector&);

  //! \brief Returns whether this integrand has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const;

  //! \brief Initializes the global node number mapping.
  virtual void initNodeMap(const std::vector<int>& nodes) { nodMap = nodes; }

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

  //! \brief Assembles contributions to the tangent stiffness and residual.
  //! \param Ktan System tangent stiffness matrix
  //! \param Res System residual vector
  virtual bool assemble(SystemMatrix& Ktan, SystemVector& Res) const;

private:
  const MortarMats& mortar; //!< Pre-integrated system-level Mortar matrices
  std::vector<int>  nodMap; //!< Nodal map from patch-level to global numbering
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
    ALElmMats(size_t nedof = 0);
    //! \brief Empty destructor.
    virtual ~ALElmMats() {}
    //! \brief Returns the element-level Newton matrix.
    virtual const Matrix& getNewtonMatrix() const;

    std::vector<size_t> iLag; //!< Indices for the active Lagrange multipliers
  };

public:
  //! \brief The contructor initializes the contact body of this integrand.
  //! \param[in] mst Pointer to the rigid body object to be in contact
  //! \param[in] mats Pre-integrated Mortar matrices associated with \a mst.
  //! \param[in] alMap Nodal map for the Lagrange multipliers
  //! \param[in] nsd Number of space dimensions
  MortarAugmentedLag(RigidBody* mst, const MortarMats& mats,
                     const std::map<int,int>& alMap, usint nsd);
  //! \brief Empty destructor.
  virtual ~MortarAugmentedLag() {}

  //! \brief Initializes the integrand for a new integration loop.
  //! \param[in] prm Nonlinear solution algorithm parameters
  //! \param[in] psol Global primary solution vector in DOF-order
  virtual void initIntegration(const TimeDomain& prm, const Vector& psol);

  //! \brief Initializes the global node number mapping.
  virtual void initNodeMap(const std::vector<int>& nodes) { nodMap = nodes; }

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

  //! \brief Assembles contributions to the tangent stiffness and residual.
  //! \param Ktan System tangent stiffness matrix
  //! \param Res System residual vector
  virtual bool assemble(SystemMatrix& Ktan, SystemVector& Res) const;

private:
  const MortarMats&        mortar; //!< Pre-integrated Mortar matrices
  const std::map<int,int>& ALmap;  //!< Nodal map for the Lagrange multipliers
  std::vector<int>         nodMap; //!< Nodal map from patch to global numbering
  Vector                   lambda; //!< Lagrange multiplier values
};

#endif
