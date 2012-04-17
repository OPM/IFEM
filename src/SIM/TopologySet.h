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

#include <string>
#include <set>
#include <map>


/*!
  \brief Struct for representing a topological item.
*/

struct TopItem
{
  size_t             patch; //!< Patch index (one-based)
  unsigned short int item;  //!< Local item index within the patch (one-based)
  unsigned short int idim;  //!< Dimension on the local item [0,3]

  //! \brief Default constructor.
  TopItem(size_t p = 1, unsigned short int i = 0, unsigned short int d = 0)
  : patch(p), item(i), idim(d) {}

  //! \brief The less-than oprator defining the ordering of topological items.
  friend bool operator<(const TopItem& a, const TopItem& b)
  {
    if (a.idim < b.idim)
      return true;
    else if (a.idim > b.idim)
      return false;
    else if (a.patch < b.patch)
      return true;
    else if (a.patch > b.patch)
      return true;

    return a.item < b.item;
  }
};

typedef std::set<TopItem> TopEntity; //!< Items defining a topological entity

typedef std::map<std::string,TopEntity> TopologySet; //!< Named topology sets

#endif
