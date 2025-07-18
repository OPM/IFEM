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
  //! Directions to consider for element stripes.
  enum StripDirection { NONE, U, V, W, ANY };

  //! \brief Default constructor.
  explicit ThreadGroups(StripDirection dir = ANY) : stripDir(dir) {}

  //! \brief Calculates a 2D thread group partitioning based on stripes.
  //! \param[in] el1 Flags non-zero knot spans in first parameter direction
  //! \param[in] el2 Flags non-zero knot spans in second parameter direction
  //! \param[in] p1 Polynomial degree in first parameter direction
  //! \param[in] p2 Polynomial degree in second parameter direction
  void calcGroups(const BoolVec& el1, const BoolVec& el2,
                  int p1, int p2);
  //! \brief Calculates a 2D thread group partitioning based on stripes.
  //! \param[in] nel1 Number of elements in the first direction
  //! \param[in] nel2 Number of elements in the second direction
  //! \param[in] minsize Minimum element strip size
  void calcGroups(int nel1, int nel2, int minsize);

  //! \brief Calculates a 3D thread group partitioning based on stripes.
  //! \param[in] el1 Flags non-zero knot spans in first parameter direction
  //! \param[in] el2 Flags non-zero knot spans in second parameter direction
  //! \param[in] el3 Flags non-zero knot spans in third parameter direction
  //! \param[in] p1 Polynomial degree in first parameter direction
  //! \param[in] p2 Polynomial degree in second parameter direction
  //! \param[in] p3 Polynomial degree in third parameter direction
  void calcGroups(const BoolVec& el1, const BoolVec& el2, const BoolVec& el3,
                  int p1, int p2, int p3);
  //! \brief Calculates a 3D thread group partitioning based on stripes.
  //! \param[in] nel1 Number of elements in the first direction
  //! \param[in] nel2 Number of elements in the second direction
  //! \param[in] nel3 Number of elements in the third direction
  //! \param[in] minsize Minimum element strip size
  void calcGroups(int nel1, int nel2, int nel3, int minsize);
  //! \brief Initializes the threading groups in case of no multi-threading.
  //! \param[in] nel Total number of elements
  void oneGroup(size_t nel);
  //! \brief Initializes the threading groups in case of a single stripe.
  //! \param[in] nel Total number of elements
  void oneStripe(size_t nel);
  //! \brief Initializes the threading groups in case of a single stripe.
  //! \param[in] elms Elements to include in stripe
  void oneStripe(const std::vector<int>& elms);

  //! \brief Maps a partitioning through a map.
  //! \details The original entry \a n in the group is mapped onto \a map[n].
  void applyMap(const IntVec& map);

  //! \brief Returns the number of groups.
  size_t size() const { return tg[1].empty() ? 1 : 2; }
  //! \brief Return true if both groups are empty.
  bool empty() const { return tg[0].empty() && tg[1].empty(); }
  //! \brief Indexing operator.
  const IntMat& operator[](int i) const { return tg[i]; }
  //! \brief Indexing operator.
  IntMat& operator[](int i) { return tg[i]; }

  //! \brief Filters current threading groups through a white-list of elements.
  ThreadGroups filter(const IntVec& elmList) const;

protected:
  //! \brief Calculates the parameter direction of the treading stripes in 2D.
  static StripDirection getStripDirection(int nel1, int nel2,
                                          int parts);
  //! \brief Calculates the parameter direction of the treading stripes in 3D.
  static StripDirection getStripDirection(int nel1, int nel2, int nel3,
                                          int parts);

  //! \brief Prints out a threading group definition.
  static void printGroup(const IntMat& group, int g);

public:
  StripDirection stripDir; //!< Actual direction to split elements

private:
  IntMat tg[2]; //!< Threading groups (always two, but the second may be empty)
};

#endif
