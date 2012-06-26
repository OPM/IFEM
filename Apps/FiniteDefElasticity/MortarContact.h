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
  //! \param[in] Reference to the FE assembly management object to use
  //! \param[in] nMast Number of master nodes
  //! \param[in] n Numner of space dimensions
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
  //! \brief Returns whether the given node node currently is in contact or not.
  bool activeNode(size_t n) const { return AA(n) > 0.0; }
  //! \brief Returns the weighter gap for the given node.
  double weightedGap(size_t n) const { return AA(n) > 0.0 ? gNA(n) : 0.0; }
  //! \brief Returns an element of the auxiliary Mortar matrix.
  double phi(size_t i, size_t j) const { return AA(j) > 0.0 ? phiA(i,j) : 0.0; }

private:
  const SAM&   sam;     //!< Data for FE assembly management
  SparseMatrix phiA;    //!< Auxiliary constants
  Vector       gNA;     //!< Weighted nodal gaps (non-positive if in contact)
  Vector       AA;      //!< Weighted nodal areas (zero if not in contact)
  usint        nsd;     //!< Number of space dimensions
  bool         contact; //!< If \e true, at least one node is in contact
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
  //! \param[in] nsd Number of space dimenensions
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

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

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
  //! \param[in] nsd Number of space dimenensions
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
  //! \brief Initializes current element for boundary integration (mixed).
  //! \param[in] MNPC1 Nodal point correspondance for the basis 1
  //! \param[in] elmInt Local integral for element
  virtual bool initElementBou(const std::vector<int>& MNPC1,
                              const std::vector<int>&, size_t,
                              LocalIntegral& elmInt)
  { return this->initElementBou(MNPC1,elmInt); }

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
  Vector gNA;  //!< Weighted nodal gaps
};

#endif
