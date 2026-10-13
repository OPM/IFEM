// $Id$
//==============================================================================
//!
//! \file GraphColoring.C
//!
//! \date Sep 28 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Multi-coloring of element conflict graphs.
//!
//==============================================================================

#include "GraphColoring.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>


/*!
  \brief Helper visiting each neighbor of an element exactly once.
  \details Elements typically share several nodes, so walking the elements of
  each node of an element visits a neighbor more than once. A visit stamp per
  element filters out the repeated visits.
*/

class GraphColoring::Neighbors
{
public:
  //! \brief The constructor initializes the visit stamps.
  explicit Neighbors(const GraphColoring& g) : graph(g), stamp(g.size(),0) {}

  //! \brief Calls \a func for each neighbor of element \a e.
  template<class Func>
  void forEach(int e, Func&& func)
  {
    if (++mark == 0) // wrapped around, reset all stamps
    {
      std::fill(stamp.begin(),stamp.end(),0u);
      mark = 1;
    }

    stamp[e] = mark;
    for (int p = graph.elmPtr[e]; p < graph.elmPtr[e+1]; p++)
    {
      const int node = graph.elmNod[p];
      for (int q = graph.nodPtr[node]; q < graph.nodPtr[node+1]; q++)
        if (int n = graph.nodElm[q]; stamp[n] != mark)
        {
          stamp[n] = mark;
          func(n);
        }
    }
  }

private:
  const GraphColoring& graph; //!< The graph to visit
  std::vector<unsigned int> stamp; //!< Last visit of each element
  unsigned int mark = 0; //!< Stamp of the current visit
};


GraphColoring::GraphColoring (const IntMat& elms, size_t nnod)
{
  elmPtr.resize(elms.size()+1,0);
  nodPtr.resize(nnod+1,0);
  for (size_t e = 0; e < elms.size(); e++)
  {
    elmPtr[e+1] = elmPtr[e] + elms[e].size();
    for (int node : elms[e])
      if (node < 0 || static_cast<size_t>(node) >= nnod)
        throw std::out_of_range("GraphColoring: Node index is out of range");
      else
        ++nodPtr[node+1];
  }

  elmNod.reserve(elmPtr.back());
  for (const IntVec& nodes : elms)
    elmNod.insert(elmNod.end(),nodes.begin(),nodes.end());

  std::partial_sum(nodPtr.begin(),nodPtr.end(),nodPtr.begin());
  nodElm.resize(nodPtr.back());
  IntVec next(nodPtr.begin(),nodPtr.end()-1);
  for (size_t e = 0; e < elms.size(); e++)
    for (int node : elms[e])
      nodElm[next[node]++] = e;
}


std::vector<int> GraphColoring::color (Algorithm alg) const
{
  IntVec order(this->size());
  std::iota(order.begin(),order.end(),0);

  switch (alg) {
  case Algorithm::FirstFit:
    return this->greedy(order);

  case Algorithm::LargestFirst:
  {
    const IntVec degree = this->degrees();
    std::stable_sort(order.begin(),order.end(),
                     [&degree](int a, int b) { return degree[a] > degree[b]; });
    return this->greedy(order);
  }

  case Algorithm::DSatur:
    return this->DSatur();

  case Algorithm::RLF:
    return this->RLF();
  }

  return {};
}


std::vector<std::vector<int>> GraphColoring::groups (const IntVec& colors)
{
  IntMat result;
  for (size_t e = 0; e < colors.size(); e++)
  {
    if (colors[e] >= static_cast<int>(result.size()))
      result.resize(colors[e]+1);
    result[colors[e]].push_back(e);
  }

  return result;
}


std::vector<int> GraphColoring::degrees () const
{
  IntVec degree(this->size(),0);
  Neighbors neighbors(*this);
  for (size_t e = 0; e < this->size(); e++)
    neighbors.forEach(e,[&degree,e](int) { ++degree[e]; });

  return degree;
}


std::vector<int> GraphColoring::greedy (const IntVec& order) const
{
  // Marking the colors of the neighbors is idempotent,
  // so there is no need to filter out repeated visits here
  IntVec colors(this->size(),-1);
  IntVec usedBy; // the element that last saw each color on a neighbor
  for (int e : order)
  {
    for (int p = elmPtr[e]; p < elmPtr[e+1]; p++)
      for (int q = nodPtr[elmNod[p]]; q < nodPtr[elmNod[p]+1]; q++)
        if (int c = colors[nodElm[q]]; c >= 0)
        {
          if (c >= static_cast<int>(usedBy.size()))
            usedBy.resize(c+1,-1);
          usedBy[c] = e;
        }

    int c = 0;
    while (c < static_cast<int>(usedBy.size()) && usedBy[c] == e)
      ++c;
    colors[e] = c;
  }

  return colors;
}


/*!
  The element with the most distinct colors among its neighbors (saturation) is
  colored next, with ties broken by the largest number of uncolored neighbors
  and then by the lowest element index. The candidates are kept in a binary
  heap which is updated whenever the priority of an element changes.
*/

std::vector<int> GraphColoring::DSatur () const
{
  const int nel = this->size();
  IntVec colors(nel,-1);
  IntVec degree = this->degrees(); // number of uncolored neighbors
  IntMat adjColors(nel); // distinct colors of the colored neighbors

  auto&& before = [&adjColors,&degree](int a, int b)
  {
    if (adjColors[a].size() != adjColors[b].size())
      return adjColors[a].size() > adjColors[b].size();
    if (degree[a] != degree[b])
      return degree[a] > degree[b];
    return a < b;
  };

  IntVec heap(nel), pos(nel);
  std::iota(heap.begin(),heap.end(),0);
  std::iota(pos.begin(),pos.end(),0);

  auto&& swapEntries = [&heap,&pos](int i, int j)
  {
    std::swap(heap[i],heap[j]);
    pos[heap[i]] = i;
    pos[heap[j]] = j;
  };
  auto&& siftUp = [&](int i)
  {
    for (int parent = (i-1)/2; i > 0 && before(heap[i],heap[parent]);
         i = parent, parent = (i-1)/2)
      swapEntries(i,parent);
  };
  auto&& siftDown = [&](int i)
  {
    for (int best = i;; i = best)
    {
      const int left = 2*i + 1;
      const int right = left + 1;
      const int size = heap.size();
      if (left < size && before(heap[left],heap[best]))
        best = left;
      if (right < size && before(heap[right],heap[best]))
        best = right;
      if (best == i)
        break;
      swapEntries(i,best);
    }
  };

  for (int i = nel/2-1; i >= 0; i--)
    siftDown(i);

  Neighbors neighbors(*this);
  while (!heap.empty())
  {
    const int e = heap.front();
    swapEntries(0,heap.size()-1);
    heap.pop_back();
    pos[e] = -1;
    siftDown(0);

    int c = 0;
    while (std::find(adjColors[e].begin(),adjColors[e].end(),c) !=
           adjColors[e].end())
      ++c;
    colors[e] = c;
    IntVec().swap(adjColors[e]);

    neighbors.forEach(e,[&](int n)
    {
      if (colors[n] >= 0)
        return;

      --degree[n];
      if (std::find(adjColors[n].begin(),adjColors[n].end(),c) ==
          adjColors[n].end())
      {
        adjColors[n].push_back(c);
        siftUp(pos[n]); // higher saturation outweighs the lower degree
      }
      else
        siftDown(pos[n]);
    });
  }

  return colors;
}


/*!
  The colors are built one at a time. Each color starts with the uncolored
  element that has the most uncolored neighbors. Then the candidate with the
  most neighbors among the elements excluded from the color is added
  repeatedly, until no candidates remain. The candidates are kept in bucket
  lists indexed by that number, so each selection and update is O(1).
*/

std::vector<int> GraphColoring::RLF () const
{
  const int nel = this->size();
  enum State { CANDIDATE, EXCLUDED, COLORED };

  IntVec colors(nel,-1);
  IntVec degree = this->degrees(); // number of uncolored neighbors
  std::vector<State> state(nel,CANDIDATE);

  // Bucket lists of the candidates, keyed by their number of excluded neighbors
  IntVec key(nel,0), next(nel,-1), prev(nel,-1), head;
  int maxKey = 0;
  auto&& remove = [&](int e)
  {
    if (prev[e] >= 0)
      next[prev[e]] = next[e];
    else
      head[key[e]] = next[e];
    if (next[e] >= 0)
      prev[next[e]] = prev[e];
  };
  auto&& insert = [&](int e)
  {
    if (key[e] >= static_cast<int>(head.size()))
      head.resize(key[e]+1,-1);
    prev[e] = -1;
    next[e] = head[key[e]];
    if (next[e] >= 0)
      prev[next[e]] = e;
    head[key[e]] = e;
    maxKey = std::max(maxKey,key[e]);
  };

  Neighbors outer(*this), inner(*this);
  auto&& addToColor = [&](int e, int c)
  {
    remove(e);
    state[e] = COLORED;
    colors[e] = c;
    outer.forEach(e,[&](int n)
    {
      if (state[n] == COLORED)
        return;

      --degree[n];
      if (state[n] == CANDIDATE)
      {
        remove(n);
        state[n] = EXCLUDED;
        inner.forEach(n,[&](int m)
        {
          if (state[m] == CANDIDATE)
          {
            remove(m);
            ++key[m];
            insert(m);
          }
        });
      }
    });
  };

  int remaining = nel;
  for (int c = 0; remaining > 0; c++)
  {
    // Reset the candidates to all uncolored elements
    head.assign(1,-1);
    maxKey = 0;
    int start = -1;
    for (int e = nel-1; e >= 0; e--)
      if (state[e] != COLORED)
      {
        state[e] = CANDIDATE;
        key[e] = 0;
        insert(e);
        if (start < 0 || degree[e] >= degree[start])
          start = e;
      }

    for (int e = start; e >= 0; --remaining)
    {
      addToColor(e,c);
      while (maxKey > 0 && head[maxKey] < 0)
        --maxKey;
      e = head[maxKey];
    }
  }

  return colors;
}
