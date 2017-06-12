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

#include "IFEM.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "Utilities.h"
#include "tinyxml.h"


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

typedef std::vector<TopSetToVertex> ConstraintVec; //!< Convenience type


//! \brief Inherit this class to equip your SIM with nodal constraints.
template<class Dim> class SIMNodalConstraint : public Dim
{
  ConstraintVec vertConstraints; //!< Registered vertex constraints

public:
  //! \brief The constructor forwards to the parent class constructor.
  SIMNodalConstraint(const std::vector<unsigned char>& unf) : Dim(unf) {}
  //! \brief Empty destructor.
  virtual ~SIMNodalConstraint() {}

protected:
  //! \brief Sets up the nodal constraints.
  virtual bool preprocessBeforeAsmInit(int&) { return this->applyConstraint(); }

  using Dim::parse;
  //! \brief Parses a data section from an XML element.
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),"constraintovertex"))
      return this->Dim::parse(elem);

    TopSetToVertex topset;
    utl::getAttribute(elem,"set",topset.topset);
    utl::getAttribute(elem,"patch",topset.patch);
    utl::getAttribute(elem,"vertex",topset.vertex);
    utl::getAttribute(elem,"comp",topset.comp);
    utl::getAttribute(elem,"basis",topset.basis);
    vertConstraints.push_back(topset);
    IFEM::cout <<"\tConstraining set \""<< topset.topset
               <<"\" to P"<< topset.patch <<" V"<< topset.vertex
               <<" in direction "<< topset.comp;
    if (topset.basis > 1)
      IFEM::cout <<" (basis "<< topset.basis <<")";
    IFEM::cout << std::endl;
    return true;
  }

private:
  //! \brief Applies the nodal constraints on the defined topology sets.
  bool applyConstraint();
};

//! \brief Specialization for 1D
template<> bool SIMNodalConstraint<SIM1D>::applyConstraint();
//! \brief Specialization for 2D
template<> bool SIMNodalConstraint<SIM2D>::applyConstraint();
//! \brief Specialization for 3D
template<> bool SIMNodalConstraint<SIM3D>::applyConstraint();

#endif
