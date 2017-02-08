// $Id$
//==============================================================================
//!
//! \file ThreadGroups.h
//!
//! \date May 15 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Threading group partitioning.
//!
//==============================================================================

#ifndef _THREAD_GROUPS_H
#define _THREAD_GROUPS_H

#include <vector>
#include <cstddef>


/*!
  \brief Class containing threading group partitioning.
*/

class ThreadGroups
{
  typedef std::vector<bool>   BoolVec; //!< List of boolean flags
  typedef std::vector<int>    IntVec;  //!< List of elements on one thread
  typedef std::vector<IntVec> IntMat;  //!< Element lists for all threads

public:
  //! \brief Calculates a 2D thread group partitioning based on strips.
  //! \param[in] el1 Flags non-zero knot spans in first parameter direction
  //! \param[in] el2 Flags non-zero knot spans in second parameter direction
  //! \param[in] p1 Polynomial degree in first parameter direction
  //! \param[in] p2 Polynomial degree in second parameter direction
  void calcGroups(const BoolVec& el1, const BoolVec& el2,
                  int p1, int p2);
  //! \brief Calculates a 2D thread group partitioning based on strips.
  //! \param[in] nel1 Number of elements in the first direction
  //! \param[in] nel2 Number of elements in the second direction
  //! \param[in] minsize Minimum element strip size
  void calcGroups(int nel1, int nel2, int minsize);

  //! \brief Calculates a 3D thread group partitioning based on strips.
  //! \param[in] el1 Flags non-zero knot spans in first parameter direction
  //! \param[in] el2 Flags non-zero knot spans in second parameter direction
  //! \param[in] el3 Flags non-zero knot spans in third parameter direction
  //! \param[in] p1 Polynomial degree in first parameter direction
  //! \param[in] p2 Polynomial degree in second parameter direction
  //! \param[in] p3 Polynomial degree in third parameter direction
  void calcGroups(const BoolVec& el1, const BoolVec& el2, const BoolVec& el3,
                  int p1, int p2, int p3);
  //! \brief Calculates a 3D thread group partitioning based on strips.
  //! \param[in] nel1 Number of elements in the first direction
  //! \param[in] nel2 Number of elements in the second direction
  //! \param[in] nel3 Number of elements in the third direction
  //! \param[in] minsize Minimum element strip size
  void calcGroups(int nel1, int nel2, int nel3, int minsize);

  //! \brief Maps a partitioning through a map.
  //! \details The original entry \a n in the group is mapped onto \a map[n].
  void applyMap(const IntVec& map);

  //! \brief Returns the number of groups.
  size_t size() const { return tg[1].empty() ? 1 : 2; }
  //! \brief Indexing operator.
  const IntMat& operator[](int i) const { return tg[i]; }
  //! \brief Indexing operator.
  IntMat& operator[](int i) { return tg[i]; }

protected:
  //! \brief Calculates the parameter direction of the treading strips in 2D.
  static int getStripDirection(int nel1, int nel2, int parts);
  //! \brief Calculates the parameter direction of the treading strips in 3D.
  static int getStripDirection(int nel1, int nel2, int nel3, int parts);

private:
  IntMat tg[2]; //!< Threading groups (always two, but the second may be empty)
};

#endif
