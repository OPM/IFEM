// $Id$
//==============================================================================
//!
//! \file SIMNodalConstraint.h
//!
//! \date Nov 4 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Base class for simulators constraining a topologyset to a given node.
//!
//==============================================================================

#ifndef _SIM_NODALCONSTRAINT_H_
#define _SIM_NODALCONSTRAINT_H_

#include "IFEM.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "Utilities.h"
#include "tinyxml.h"

class ASMbase;
class ASMs1D;
class ASMs2D;
class ASMs3D;
class ASMu2D;


//! \brief Describes a topologyset constrained to a vertex.
struct TopSetToVertex {
  std::string topset; //!< Topologyset to be constrained
  size_t basis;       //!< Basis to constrain
  size_t patch;       //!< Patch number of vertex
  size_t vertex;      //!< Vertex number
  size_t comp;        //!< Component to constrain

  //! \brief Default constructor.
  TopSetToVertex() : basis(1), patch(1), vertex(1), comp(1) {}
};


//! \brief Helper class for constraining entities to a node
template<class Dim>
class NodalConstraintProcessor {
public:
  //! \brief Convenience typedef
  typedef std::vector<TopSetToVertex> ConstraintVec;

  //! \brief Apply nodal constraints
  //! \param sim The simulator to apply the constraints to
  //! \param myEntitys The topologysets of the simulator
  //! \param vertConstraints Constraints to apply
  static void apply(Dim& sim, TopologySet& myEntitys,
                    const ConstraintVec& vertConstraints);
};

//! \brief Specialization for 1D
template<>
void NodalConstraintProcessor<SIM1D>::apply(SIM1D& sim, TopologySet& myEntitys,
                                            const ConstraintVec& vConstr);


//! \brief Specialization for 2D
template<>
void NodalConstraintProcessor<SIM2D>::apply(SIM2D& sim, TopologySet& myEntitys,
                                            const ConstraintVec& vConstr);

//! \brief Specialization for 3D
template<>
void NodalConstraintProcessor<SIM3D>::apply(SIM3D& sim, TopologySet& myEntitys,
                                            const ConstraintVec& vConstr);


//! \brief Inherit this class to equip your SIM with nodal constraints.
template<class Dim>
class SIMNodalConstraint : public Dim {
public:
  //! \brief Default constructor.
  //! \param[in] n1 Dimension of the primary solution field
  //! \param[in] check If \e true, ensure the model is in a right-hand system
  SIMNodalConstraint(const SIMbase::CharVec& unf,
                     bool checkRHS=false) :
    Dim(unf,checkRHS) {}

  //! \brief Empty destructor
  virtual ~SIMNodalConstraint() {}

  //! \copydoc SIMbase::preprocessBeforeAsmInit(int&)
  //! \details Sets up the nodal constraints.
  virtual bool preprocessBeforeAsmInit(int&)
  {
    NodalConstraintProcessor<Dim>::apply(*this, this->myEntitys,
                                         vertConstraints);
    return true;
  }

  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    if (!strcasecmp(elem->Value(),"constraintovertex")) {
      vertConstraints.resize(vertConstraints.size()+1);
      utl::getAttribute(elem,"set",vertConstraints.back().topset);
      utl::getAttribute(elem,"patch",vertConstraints.back().patch);
      utl::getAttribute(elem,"vertex",vertConstraints.back().vertex);
      utl::getAttribute(elem,"comp",vertConstraints.back().comp);
      utl::getAttribute(elem,"basis",vertConstraints.back().basis);
      IFEM::cout << "\tConstraining set \"" << vertConstraints.back().topset
                 << "\" to P"<< vertConstraints.back().patch << "V"
                 << vertConstraints.back().vertex
                 << " in direction " << vertConstraints.back().comp;
      if (vertConstraints.back().basis > 1)
        IFEM::cout << " (basis " << vertConstraints.back().basis << ")";
      IFEM::cout << std::endl;
    } else
      return this->Dim::parse(elem);

    return true;
  }
protected:
  std::vector<TopSetToVertex> vertConstraints; //!< Registered vertex constraints
};

#endif
