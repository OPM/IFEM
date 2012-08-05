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
#include "MatVec.h"
#include <map>

class TiXmlElement;
class RigidBody;
class ASMbase;
class IntegrandBase;
class SystemMatrix;
class SystemVector;
class SAM;
class MPC;
class VTF;


/*!
  \brief Class for administration of finite deformation contact analysis.
  \details The class incapsulates data and methods for the inclusion of one
  or more rigid bodies that are in contact with the flexible body.
  The class has no public members and can only be used to augment functionality
  to the main FE solution driver class through inheritance.
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

  //! \brief Returns whether contact is included in this simulation or not.
  bool withContact() const { return contactMethod > NONE && !myBodies.empty(); }

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
  //! \param[in] entitys Set of all named topological entities in the model
  //! \param[in] model The flexible FE model
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

  //! \brief Adds Lagrange multipliers as unknowns to the specified patch.
  //! \param patch The patch that will recieve the Lagrange multipliers
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

  //! \brief Renumbers the nodes of the contact bodies.
  //! \param[in] old2new Old-to-new node number mapping
  void renumberContactBodies(const std::map<int,int>& old2new);
  //! \brief Updates the positions of the contact bodies.
  //! \param[in] displ Current total displacement vector in DOF order
  bool updateContactBodies(const RealArray& displ);

  //! \brief Assembles contributions to the tangent stiffness and residual.
  //! \param[in] problem The integrand containing the tangent contributions
  //! \param Ktan System tangent stiffness matrix
  //! \param Res System residual vector
  //!
  //! \details This method assembles direct nodal contributions to the
  //! system matrices emanating from the Mortar method formulation,
  //! which cannot be assembled through the normal element assembly loop.
  bool assembleMortarTangent(const IntegrandBase* problem,
                             SystemMatrix* Ktan, SystemVector* Res);

  //! \brief Writes the geometry of the rigid bodies to the VTF-file.
  //! \param vtf The VTF-file to receive the geometry data
  //! \param nBlock Running result block counter
  bool writeGlvBodies(VTF* vtf, int& nBlock);

  //! \brief Writes the current position of the rigid bodies to the VTF-file.
  //! \param vtf The VTF-file to receive the position data
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  bool writeGlvBodyMovements(VTF* vtf, int iStep, int& nBlock) const;

  //! \brief Dumps total reaction forces for each rigid body to a given stream.
  //! \param[in] sam Auxiliary data for FE assembly management
  //! \param[in] RF Compressed reaction force vector for the entire model
  //! \param os Output stream to write the reaction forces to
  //! \param[in] precision Number of digits after the decimal point
  void printBodyReactions(const SAM& sam, const Vector& RF,
                          std::ostream& os, std::streamsize precision) const;

private:
  std::vector<RigidBody*>    myBodies; //!< All rigid bodies of the model
  std::vector<ScalarFunc*>   myFuncs;  //!< Functions for prescribed movements
  std::map<MPC*,ScalarFunc*> dMap;     //!< MPC equation to functions map
  std::map<int,int>          myALp;    //!< Augmented Lagrange multiplier map

  Method contactMethod; //!< Contact formulation option
};

#endif
