//==============================================================================
//!
//! \file TestGraphColoring.C
//!
//! \date Sep 28 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for multi-coloring of element conflict graphs.
//!
//==============================================================================

#include "GraphColoring.h"

#include "Catch2Support.h"

#include <algorithm>
#include <random>
#include <set>
#include <stdexcept>

using IntVec = std::vector<int>;
using IntMat = std::vector<IntVec>;
using Algorithm = GraphColoring::Algorithm;


namespace
{

const auto allAlgorithms = {
  Algorithm::FirstFit,
  Algorithm::LargestFirst,
  Algorithm::DSatur,
  Algorithm::RLF
};


//! \brief Returns the number of colors in a coloring.
int numColors(const IntVec& colors)
{
  return colors.empty() ? 0 : 1 + *std::max_element(colors.begin(),colors.end());
}


//! \brief Checks that no two elements sharing a node have the same color.
void requireValid(const IntMat& elms, const IntVec& colors)
{
  REQUIRE(colors.size() == elms.size());
  REQUIRE(std::all_of(colors.begin(),colors.end(),[](int c) { return c >= 0; }));

  size_t conflicts = 0;
  for (size_t e = 0; e < elms.size(); e++)
  {
    const std::set<int> nodes(elms[e].begin(),elms[e].end());
    for (size_t f = e+1; f < elms.size(); f++)
      if (colors[e] == colors[f])
        for (int node : elms[f])
          conflicts += nodes.count(node);
  }
  REQUIRE(conflicts == 0);
}


//! \brief Returns the connectivity of a structured grid of bilinear elements.
IntMat quadGrid(int nx, int ny)
{
  IntMat elms;
  for (int j = 0; j < ny; j++)
    for (int i = 0; i < nx; i++)
      elms.push_back({ i + j*(nx+1), i+1 + j*(nx+1),
                       i+1 + (j+1)*(nx+1), i + (j+1)*(nx+1) });
  return elms;
}


/*!
  \brief Returns the crown graph on 2n vertices as element connectivities.
  \details Element 2i is joined to element 2j+1 for all i != j, each edge by
  a node of its own. The graph is bipartite, but greedy coloring in element
  order needs n colors.
*/

IntMat crownGraph(int n, size_t& nnod)
{
  IntMat elms(2*n);
  nnod = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if (i != j)
      {
        elms[2*i].push_back(nnod);
        elms[2*j+1].push_back(nnod++);
      }
  return elms;
}

}


TEST_CASE("TestGraphColoring.Empty")
{
  for (Algorithm alg : allAlgorithms)
    REQUIRE(GraphColoring({},0).color(alg).empty());
}


TEST_CASE("TestGraphColoring.NoConflicts")
{
  const IntMat elms = {{0}, {1}, {}, {2,2}};
  for (Algorithm alg : allAlgorithms)
    REQUIRE(GraphColoring(elms,3).color(alg) == IntVec(4,0));
}


TEST_CASE("TestGraphColoring.NodeOutOfRange")
{
  REQUIRE_THROWS_AS(GraphColoring({{0,3}},3),std::out_of_range);
  REQUIRE_THROWS_AS(GraphColoring({{-1}},3),std::out_of_range);
}


TEST_CASE("TestGraphColoring.QuadGrid")
{
  const IntMat elms = quadGrid(4,4);
  const GraphColoring graph(elms,25);
  for (Algorithm alg : allAlgorithms)
  {
    const IntVec colors = graph.color(alg);
    requireValid(elms,colors);
    REQUIRE(numColors(colors) == 4);
  }

  // Greedy coloring in element order checkerboards the grid in 2x2 blocks
  const IntMat groups = GraphColoring::groups(graph.color(Algorithm::FirstFit));
  REQUIRE(groups.size() == 4);
  REQUIRE(groups[0] == IntVec{0, 2, 8, 10});
  REQUIRE(groups[1] == IntVec{1, 3, 9, 11});
  REQUIRE(groups[2] == IntVec{4, 6, 12, 14});
  REQUIRE(groups[3] == IntVec{5, 7, 13, 15});
}


TEST_CASE("TestGraphColoring.Crown")
{
  size_t nnod;
  const IntMat elms = crownGraph(6,nnod);
  const GraphColoring graph(elms,nnod);
  for (Algorithm alg : allAlgorithms)
    requireValid(elms,graph.color(alg));

  // All elements have the same degree, so these just follow element order
  REQUIRE(numColors(graph.color(Algorithm::FirstFit)) == 6);
  REQUIRE(numColors(graph.color(Algorithm::LargestFirst)) == 6);

  // DSatur is exact on bipartite graphs, and so is RLF on this one
  REQUIRE(numColors(graph.color(Algorithm::DSatur)) == 2);
  REQUIRE(numColors(graph.color(Algorithm::RLF)) == 2);
}


TEST_CASE("TestGraphColoring.LargestFirst")
{
  // The path a-b-c-d with the elements in the order a, d, b, c.
  // Coloring the two ends first costs a third color in the middle.
  const IntMat elms = {{0}, {2}, {0,1}, {1,2}};
  const GraphColoring graph(elms,3);
  REQUIRE(graph.color(Algorithm::FirstFit) == IntVec{0, 0, 1, 2});
  REQUIRE(graph.color(Algorithm::LargestFirst) == IntVec{1, 0, 0, 1});
}


TEST_CASE("TestGraphColoring.Random")
{
  std::mt19937 rng(1);
  for (int t = 0; t < 50; t++)
  {
    const int nel = 1 + t*7;
    const int nnod = 1 + t*3;
    std::uniform_int_distribution<int> node(0,nnod-1);
    std::uniform_int_distribution<int> count(0,5);
    IntMat elms(nel);
    for (IntVec& nodes : elms)
      for (int n = count(rng); n > 0; n--)
        nodes.push_back(node(rng)); // repeated nodes are allowed

    const GraphColoring graph(elms,nnod);
    for (Algorithm alg : allAlgorithms)
    {
      const IntVec colors = graph.color(alg);
      requireValid(elms,colors);

      // Every color is used, and grouping preserves all elements
      const IntMat groups = GraphColoring::groups(colors);
      REQUIRE(static_cast<int>(groups.size()) == numColors(colors));
      size_t total = 0;
      for (const IntVec& group : groups)
      {
        REQUIRE(!group.empty());
        REQUIRE(std::is_sorted(group.begin(),group.end()));
        total += group.size();
      }
      REQUIRE(total == elms.size());
    }
  }
}
