// $Id$
//==============================================================================
//!
//! \file GraphColoring.h
//!
//! \date Sep 28 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Multi-coloring of element conflict graphs.
//!
//==============================================================================

#ifndef _GRAPH_COLORING_H
#define _GRAPH_COLORING_H

#include <cstddef>
#include <vector>


/*!
  \brief Class for multi-coloring the conflict graph of a set of elements.

  \details Each element writes to a set of nodes, and two elements are in
  conflict if they share at least one node. A coloring assigns a color to each
  element such that no two conflicting elements have the same color, i.e.,
  the elements of one color can be processed concurrently.

  The graph is stored as the element-to-node connectivity and its transpose,
  and the neighbors of an element are found through the nodes it shares with
  them. This needs far less memory than an explicit adjacency graph for
  high-order discretizations, where elements share many nodes.
*/

class GraphColoring
{
  using IntVec = std::vector<int>;    //!< General integer vector
  using IntMat = std::vector<IntVec>; //!< General 2D integer matrix

public:
  //! \brief Enum defining the available coloring algorithms.
  enum class Algorithm
  {
    FirstFit,     //!< Greedy coloring in element order
    LargestFirst, //!< Greedy coloring in order of decreasing degree
    DSatur,       //!< Greedy coloring by largest saturation degree first
    RLF           //!< Recursive largest first
  };

  //! \brief The constructor sets up the conflict graph.
  //! \param[in] elms The nodes each element writes to
  //! \param[in] nnod Number of nodes, all node indices must be in [0,nnod)
  GraphColoring(const IntMat& elms, size_t nnod);

  //! \brief Returns the number of elements in the graph.
  size_t size() const { return elmPtr.size() - 1; }

  //! \brief Colors the elements.
  //! \param[in] alg The coloring algorithm to use
  //! \return The color of each element, numbered from zero
  IntVec color(Algorithm alg) const;

  //! \brief Groups the elements by color.
  //! \param[in] colors The color of each element
  //! \return The elements of each color, in ascending order
  static IntMat groups(const IntVec& colors);

private:
  class Neighbors; //!< Helper visiting the neighbors of an element

  //! \brief Greedy coloring visiting the elements in the given order.
  IntVec greedy(const IntVec& order) const;
  //! \brief Colors the elements with the DSatur algorithm.
  IntVec DSatur() const;
  //! \brief Colors the elements with the recursive largest first algorithm.
  IntVec RLF() const;
  //! \brief Returns the number of neighbors of each element.
  IntVec degrees() const;

  IntVec elmPtr; //!< Offsets into \ref elmNod for each element
  IntVec elmNod; //!< The nodes of each element
  IntVec nodPtr; //!< Offsets into \ref nodElm for each node
  IntVec nodElm; //!< The elements of each node
};

#endif
