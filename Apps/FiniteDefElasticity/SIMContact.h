// $Id$
//==============================================================================
//!
//! \file SIMContact.h
//!
//! \date Jun 05 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution drivers for finite deformation contact analysis.
//!
//==============================================================================

#ifndef _SIM_CONTACT_H
#define _SIM_CONTACT_H

#include "TopologySet.h"
#include "Property.h"
#include "Function.h"
#include <vector>
#include <map>

class TiXmlElement;
class RigidBody;
class IntegrandBase;
class ASMbase;
class MPC;
class SAM;


class SIMContact
{
protected:
  //! \brief The default constructor is protected to allow sub-classes only.
  SIMContact() {}
  //! \brief Empty destructor.
  virtual ~SIMContact();

  //! \brief Creates a property set for contact condition on an entity set.
  //! \param[in] slave Name of the slave entity set
  //! \param[out] code Property code associated with the contact set
  virtual bool createContactSet(const std::string& slaveSet, int& code) = 0;

private:
  //! \brief Establishes nodal connectivity for the contact elements.
  //! \param[in] master The master rigid body object
  //! \param[in] model The flexible FE model
  //! \param[in] entitys Set of all named topological entities in the model
  //! \param[in] slave Name of the slave entity set
  bool addContactElms(RigidBody* master,
		      const std::vector<ASMbase*>& model,
		      const TopologySet& entitys,
		      const std::string& slave);

protected:
  //! Property code to integrand map
  typedef std::multimap<int,IntegrandBase*> IntegrandMap;

  //! \brief Parses a subelement of the \a contact XML-tag.
  //! \param[in] elem The XML element to parse
  //! \param[in] model The flexible FE model
  //! \param[in] entitys Set of all named topological entities in the model
  bool parseContactTag(const TiXmlElement* elem,
		       const std::vector<ASMbase*>& model,
		       const TopologySet& entitys);

  //! \brief Performs some pre-processing tasks on the contact bodies.
  //! \param All integrands of the contact problem (including the main one)
  //! \param[in] sam Auxiliary data for FE assembly management
  //! \param[in] nsd Number of space dimensions
  bool preprocessContact(IntegrandMap& itgs, const SAM& sam, size_t nsd);

  //! \brief Updates the time-dependent in-homogeneous Dirichlet coefficients.
  //! \param[in] time Current time
  void updateDirichlet(double time);

  //! \brief Updates the positions of the contact bodies.
  //! \param[in] displ Current total displacement vector in DOF order
  bool updateContactBodies(const std::vector<double>& displ);

private:
  //! Rigid/flexible body contact pair representation
  typedef std::pair<RigidBody*,ASMbase*> ContactPair;

  std::vector<ContactPair>   myContacts; //!< Actual contact pair definitions
  std::vector<RigidBody*>    myBodies;   //!< All rigid bodies of the model
  std::vector<ScalarFunc*>   myFuncs;    //!< Functions for prescribed movements
  std::map<MPC*,ScalarFunc*> dMap;       //!< MPC equation to functions map
};

#endif
