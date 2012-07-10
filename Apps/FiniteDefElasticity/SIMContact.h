// $Id$
//==============================================================================
//!
//! \file SIMContact.h
//!
//! \date Jun 05 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Administration of finite deformation contact analysis.
//!
//==============================================================================

#ifndef _SIM_CONTACT_H
#define _SIM_CONTACT_H

#include "TopologySet.h"
#include "Function.h"
#include <vector>
#include <map>

class TiXmlElement;
class RigidBody;
class ASMbase;
class IntegrandBase;
class SAM;
class MPC;


/*!
  \brief Administration of finite deformation contact analysis.
*/

class SIMContact
{
  //! \brief Enum defining the available contact formulations.
  enum Method { NONE, PENALTY, AUGMENTED_LAGRANGE };

protected:
  //! \brief The default constructor is protected to allow sub-classes only.
  SIMContact() : contactMethod(NONE) {}
  //! \brief The destructor deletes the rigid body objects and scalar functions.
  virtual ~SIMContact();

  //! \brief Creates a property set for contact condition on an entity set.
  //! \param[in] slaveSet Name of the slave boundary entity set
  //! \param[out] code Property code associated with the contact set
  virtual bool createContactSet(const std::string& slaveSet, int& code) = 0;

private:
  //! Rigid/flexible body contact pair representation
  typedef std::pair<RigidBody*,ASMbase*> ContactPair;

  //! \brief Predicate class used when searching in contact pairs.
  class hasBody
  {
    const RigidBody* body; //!< The rigid body to search for
  public:
    //! \brief Constructor initializing the rigid body pointer.
    hasBody(const RigidBody* b) : body(b) {}
    //! \brief Returns \e true if \a p contains \a body.
    bool operator() (const ContactPair& p) const { return p.first == body; }
  };

  //! \brief Establishes nodal connectivity for the contact elements.
  //! \param contacts Array of contact pair definitions
  //! \param[in] master The master rigid body object
  //! \param[in] slave Name of the slave entity set
  //! \param[in] model The flexible FE model
  //! \param[in] entitys Set of all named topological entities in the model
  bool addContactElms(std::vector<ContactPair>& contacts, RigidBody* master,
		      const std::string& slave, const TopologySet& entitys,
		      const std::vector<ASMbase*>& model);

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

  //! \brief Adds Lagrangian multipliers as unknowns to the specified patch.
  //! \param patch The patch to recieve Lagrangian multipliers
  //! \param ngnod Total number of nodes in the model
  bool addLagrangeMultipliers(ASMbase* patch, int& ngnod);

  //! \brief Performs some pre-processing tasks on the contact bodies.
  //! \param itgs All integrands of the contact problem (including the main one)
  //! \param[in] sam Auxiliary data for FE assembly management
  //! \param[in] nsd Number of space dimensions
  bool preprocessContact(IntegrandMap& itgs, const SAM& sam, size_t nsd);

  //! \brief Updates time-dependent in-homogeneous Dirichlet coefficients.
  //! \param[in] time Current time
  void updateDirichlet(double time);

  //! \brief Updates the positions of the contact bodies.
  //! \param[in] displ Current total displacement vector in DOF order
  bool updateContactBodies(const std::vector<double>& displ);

private:
  std::vector<RigidBody*>    myBodies; //!< All rigid bodies of the model
  std::vector<ScalarFunc*>   myFuncs;  //!< Functions for prescribed movements
  std::map<MPC*,ScalarFunc*> dMap;     //!< MPC equation to functions map
  std::map<int,int>          myALp;    //!< Augmented Lagrangian multiplier map

  Method contactMethod; //!< Contact formulation option
};

#endif
