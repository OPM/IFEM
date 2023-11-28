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

#ifndef _SIM_NODAL_CONSTRAINT_H_
#define _SIM_NODAL_CONSTRAINT_H_

#include <string>
#include <vector>

class TiXmlElement;


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

using ConstraintVec = std::vector<TopSetToVertex>; //!< Convenience type


//! \brief Inherit this class to equip your SIM with nodal constraints.
template<class Dim> class SIMNodalConstraint : public Dim
{
  ConstraintVec vertConstraints; //!< Registered vertex constraints

public:
  //! \brief The constructor forwards to the parent class constructor.
  explicit SIMNodalConstraint(const std::vector<unsigned char>& unf) : Dim(unf) {}
  //! \brief Empty destructor.
  virtual ~SIMNodalConstraint() {}

protected:
  //! \brief Sets up the nodal constraints.
  virtual bool preprocessBeforeAsmInit(int&) { return this->applyConstraint(); }

  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem);

private:
  //! \brief Applies the nodal constraints on the defined topology sets.
  bool applyConstraint();
};

#endif
