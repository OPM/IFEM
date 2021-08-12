// $Id$
//==============================================================================
//!
//! \file NodeVecFunc.h
//!
//! \date Apr 2 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Vector function wrapper for a nodal field.
//!
//==============================================================================

#ifndef _NODE_VEC_FUNC_H
#define _NODE_VEC_FUNC_H

#include "Function.h"
#include <map>

class SIMbase;


/*!
  \brief A class that wraps a nodal field as a vector-valued spatial function.
*/

class NodeVecFunc : public VecFunc
{
public:
  //! \brief The constructor initializes the references.
  NodeVecFunc(const SIMbase& m, const std::vector<double>* v = nullptr)
    : idMap(dummy), model(m), value(v) {}
  //! \brief This constructor provides a node number map in addition.
  NodeVecFunc(const SIMbase& m, const std::vector<double>* v,
              const std::map<int,int>& nodeIdMap)
    : idMap(nodeIdMap), model(m), value(v) {}
  //! \brief This constructor provides a node number map in addition.
  NodeVecFunc(const SIMbase& m, const std::vector<double>& v,
              const std::map<int,int>& nodeIdMap)
    : idMap(nodeIdMap), model(m), value(&v) {}
  //! \brief Empty destructor.
  virtual ~NodeVecFunc() {}

  //! \brief Returns whether the function is identically zero or not.
  virtual bool isZero() const;

  //! \brief Returns whether the function is time-independent or not.
  virtual bool isConstant() const { return false; }

protected:
  //! \brief Evaluates the function at the point \a X.
  virtual Vec3 evaluate(const Vec3& X) const;

  //! \brief Returns the node index (if any) matching the given coordinates.
  std::pair<int,int> getPointIndex(const Vec3& Xp) const;

private:
  const std::map<int,int>    dummy; //!< Dummy empty map
  const std::map<int,int>&   idMap; //!< Map of node indices
  mutable std::map<Vec3,int> ptMap; //!< Map of evaluated nodal points

protected:
  const SIMbase&             model; //!< FE model on which the field is defined
  const std::vector<double>* value; //!< The nodal field values
};

#endif
