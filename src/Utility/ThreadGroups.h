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

  \details The elements are divided into tasks, which are sequences of
  elements assembled one after the other by one thread, and the tasks are
  divided into colors. The colors are assembled one after the other, while
  the tasks of a color are assembled concurrently. Two tasks of one color must
  therefore never write to the same equation.

  Structured patches use tiles of elements as tasks, for their locality, while
  unstructured patches use one task per element.
*/

class ThreadGroups
{
  typedef std::vector<bool>   BoolVec; //!< List of boolean flags
  typedef std::vector<int>    IntVec;  //!< List of elements of a task
  typedef std::vector<IntVec> IntMat;  //!< Element lists for all tasks

public:
  //! \brief Assembles the elements in sequence, without multi-threading.
  //! \param[in] nel Total number of elements
  //! \param[in] elms Elements to include, all if empty, none if elms[0] < 0
  void sequential(size_t nel, const IntVec& elms = {});
  //! \brief Assembles all elements concurrently.
  //! \details This is used for integrands that are thread safe by themselves.
  //! \param[in] nel Total number of elements
  //! \param[in] elms Elements to include, all if empty, none if elms[0] < 0
  void concurrent(size_t nel, const IntVec& elms = {});

  //! \brief Sets the groups from a coloring of tasks.
  //! \param[in] tasks The element lists of the tasks
  //! \param[in] colors The color of each task
  void setColors(const IntMat& tasks, const IntVec& colors);

  //! \brief Divides a 2D grid of knot spans into tiles.
  //! \param[in] el1 Flags non-zero knot spans in first parameter direction
  //! \param[in] el2 Flags non-zero knot spans in second parameter direction
  //! \param[in] b1 Minimum non-zero knot spans of a tile in first direction
  //! \param[in] b2 Minimum non-zero knot spans of a tile in second direction
  //! \param[out] parity The parity color of each tile, if wanted
  //! \return The elements of each tile, tiles in the natural order of the grid
  //!
  //! \details The widths of the tiles in a direction differ by at most one,
  //! which balances the work of the tiles of a color.
  //! Tiles of the same parity are separated by at least one tile.
  //! Elements of such tiles thus share no basis functions of polynomial degree
  //! \a p, where \a b is at least \a p, in the absence of constraints.
  static IntMat tiles(const BoolVec& el1, const BoolVec& el2, int b1, int b2,
                      IntVec* parity = nullptr);
  //! \brief Divides a 3D grid of knot spans into tiles.
  //! \param[in] el1 Flags non-zero knot spans in first parameter direction
  //! \param[in] el2 Flags non-zero knot spans in second parameter direction
  //! \param[in] el3 Flags non-zero knot spans in third parameter direction
  //! \param[in] b1 Minimum non-zero knot spans of a tile in first direction
  //! \param[in] b2 Minimum non-zero knot spans of a tile in second direction
  //! \param[in] b3 Minimum non-zero knot spans of a tile in third direction
  //! \param[out] parity The parity color of each tile, if wanted
  //! \return The elements of each tile, tiles in the natural order of the grid
  static IntMat tiles(const BoolVec& el1, const BoolVec& el2,
                      const BoolVec& el3, int b1, int b2, int b3,
                      IntVec* parity = nullptr);

  //! \brief Maps a partitioning through a map.
  //! \details The original entry \a n in the group is mapped onto \a map[n].
  void applyMap(const IntVec& map);

  //! \brief Returns the number of colors.
  size_t size() const { return tg.size(); }
  //! \brief Returns true if there are no colors.
  bool empty() const { return tg.empty(); }
  //! \brief Returns the total number of elements in all tasks.
  size_t noElms() const;
  //! \brief Returns the tasks of a color.
  const IntMat& operator[](size_t i) const { return tg[i]; }

  //! \brief Filters current threading groups through a white-list of elements.
  ThreadGroups filter(const IntVec& elmList) const;

  //! \brief Prints statistics of the colors.
  //! \param[in] listAllSizes If \e true, list the element count of each color
  void analyze(bool listAllSizes = false) const;

private:
  std::vector<IntMat> tg; //!< The tasks of each color
};

#endif
