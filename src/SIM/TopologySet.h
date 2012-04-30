// $Id$
//==============================================================================
//!
//! \file TopologySet.h
//!
//! \date Apr 17 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of topological entities.
//!
//==============================================================================

#ifndef _TOPOLOGY_SET_H
#define _TOPOLOGY_SET_H

#include <iostream>
#include <string>
#include <set>
#include <map>
#include <cmath>


/*!
  \brief Struct for representing a topological item.
*/

struct TopItem
{
  size_t    patch; //!< Patch index (one-based)
  short int item;  //!< Local item index within the patch (one-based)
  short int idim;  //!< Dimension on the local item [-3,3]

  //! \brief Default constructor.
  TopItem(size_t p = 0, short int i = 0, short int d = 0)
  : patch(p), item(i), idim(d) {}

  //! \brief The less-than operator defining the ordering of topological items.
  friend bool operator<(const TopItem& a, const TopItem& b)
  {
    if (std::abs(a.idim) < std::abs(b.idim))
      return true;
    else if (std::abs(a.idim) > std::abs(b.idim))
      return false;
    else if (a.patch < b.patch)
      return true;
    else if (a.patch > b.patch)
      return false;

    return a.item < b.item;
  }

  //! \brief Output stream operator.
  friend std::ostream& operator<<(std::ostream& os, const TopItem& top)
  {
    return os <<" ("<< top.patch <<","<< top.item <<","<< top.idim <<"D)";
  }
};

typedef std::set<TopItem> TopEntity; //!< Items defining a topological entity

typedef std::map<std::string,TopEntity> TopologySet; //!< Named topology sets

#endif
