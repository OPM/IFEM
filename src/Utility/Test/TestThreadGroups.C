//==============================================================================
//!
//! \file TestThreadGroups.C
//!
//! \date Oct 13 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for threading group partitioning.
//!
//==============================================================================

#include "ThreadGroups.h"

#include "Catch2Support.h"

#include <set>
#include <vector>

using IntVec = std::vector<int>;
using IntMat = std::vector<IntVec>;


TEST_CASE("TestThreadGroups.Sequential")
{
  ThreadGroups groups;
  groups.sequential(5);
  REQUIRE(groups.size() == 1);
  REQUIRE(groups[0] == IntMat{{0, 1, 2, 3, 4}});

  groups.sequential(5, {3, 1});
  REQUIRE(groups.size() == 1);
  REQUIRE(groups[0] == IntMat{{3, 1}});

  groups.sequential(5, {-1});
  REQUIRE(groups.empty());
}


TEST_CASE("TestThreadGroups.Concurrent")
{
  ThreadGroups groups;
  groups.concurrent(3);
  REQUIRE(groups.size() == 1);
  REQUIRE(groups[0] == IntMat{{0}, {1}, {2}});

  groups.concurrent(5, {4, 2});
  REQUIRE(groups.size() == 1);
  REQUIRE(groups[0] == IntMat{{4}, {2}});

  groups.concurrent(5, {-1});
  REQUIRE(groups.empty());
}


TEST_CASE("TestThreadGroups.Tiles2D")
{
  const std::vector<bool> el(8, true);
  IntVec parity;
  const IntMat tiles = ThreadGroups::tiles(el, el, 2, 2, &parity);

  REQUIRE(tiles.size() == 16);
  REQUIRE(tiles[0] == IntVec{0, 1, 8, 9});
  REQUIRE(tiles[1] == IntVec{2, 3, 10, 11});
  REQUIRE(tiles[4] == IntVec{16, 17, 24, 25});
  REQUIRE(tiles[15] == IntVec{54, 55, 62, 63});
  REQUIRE(parity == IntVec{0, 1, 0, 1, 2, 3, 2, 3, 0, 1, 0, 1, 2, 3, 2, 3});

  // Every element is in exactly one tile
  std::set<int> elms;
  for (const IntVec& tile : tiles)
    elms.insert(tile.begin(), tile.end());
  REQUIRE(elms.size() == 64);
  REQUIRE(*elms.rbegin() == 63);

  ThreadGroups groups;
  groups.setColors(tiles, parity);
  REQUIRE(groups.size() == 4);
  for (size_t c = 0; c < groups.size(); c++)
    REQUIRE(groups[c].size() == 4);
  REQUIRE(groups[3][0] == IntVec{18, 19, 26, 27});
}


TEST_CASE("TestThreadGroups.TilesBalanced")
{
  // Ten spans in tiles of at least three give widths of three and four
  const std::vector<bool> el1(10, true);
  const std::vector<bool> el2{true};
  const IntMat tiles = ThreadGroups::tiles(el1, el2, 3, 1);
  REQUIRE(tiles == IntMat{{0, 1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
}


TEST_CASE("TestThreadGroups.TilesZeroSpans")
{
  // The zero spans stay with the tile preceding them
  const std::vector<bool> el1{true, false, true, true, false, true};
  const std::vector<bool> el2{true};
  const IntMat tiles = ThreadGroups::tiles(el1, el2, 2, 1);
  REQUIRE(tiles == IntMat{{0, 1, 2}, {3, 4, 5}});
}


TEST_CASE("TestThreadGroups.Tiles3D")
{
  const std::vector<bool> el(4, true);
  IntVec parity;
  const IntMat tiles = ThreadGroups::tiles(el, el, el, 2, 2, 2, &parity);

  REQUIRE(tiles.size() == 8);
  REQUIRE(parity == IntVec{0, 1, 2, 3, 4, 5, 6, 7});
  REQUIRE(tiles[0] == IntVec{0, 1, 4, 5, 16, 17, 20, 21});
  REQUIRE(tiles[7] == IntVec{42, 43, 46, 47, 58, 59, 62, 63});
}


TEST_CASE("TestThreadGroups.SetColors")
{
  // Empty tasks are dropped, and so are colors left without tasks
  ThreadGroups groups;
  groups.setColors({{0, 1}, {}, {2}, {3}}, {2, 0, 2, 5});
  REQUIRE(groups.size() == 2);
  REQUIRE(groups[0] == IntMat{{0, 1}, {2}});
  REQUIRE(groups[1] == IntMat{{3}});
}


TEST_CASE("TestThreadGroups.Filter")
{
  ThreadGroups groups;
  groups.setColors({{0, 1, 2}, {3, 4}, {5}}, {0, 0, 1});

  const ThreadGroups filtered = groups.filter({4, 1, 2});
  REQUIRE(filtered.size() == 1);
  REQUIRE(filtered[0] == IntMat{{1, 2}, {4}});
}


TEST_CASE("TestThreadGroups.ApplyMap")
{
  ThreadGroups groups;
  groups.setColors({{0, 1}, {2}}, {0, 1});
  groups.applyMap({10, 11, 12});
  REQUIRE(groups[0] == IntMat{{10, 11}});
  REQUIRE(groups[1] == IntMat{{12}});
}
