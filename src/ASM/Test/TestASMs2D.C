//==============================================================================
//!
//! \file TestASMs2D.C
//!
//! \date Feb 14 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for structured 2D spline FE models.
//!
//==============================================================================

#include "ASMSquare.h"
#include "ASMTrap.h"
#include "ASM2DTests.h"

#include <memory>
#include <set>
#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace {

class ASMdegenerate2D : public ASMs2D
{
public:
  explicit ASMdegenerate2D(int iedge)
  {
    // Create a degenerated patch where edge "iedge" is collapsed into a vertex
    std::stringstream geo; // -- Xi ----    --- Eta ----
    geo <<"200 1 0 0 2 0"<<" 2 2 0 0 1 1"<<" 2 2 0 0 1 1"<<" 0 0 ";
    switch (iedge) {
    case 1:
      geo <<"3 0 0 0 3 1\n"; // P3=P1
      break;
    case 2:
      geo <<"3 0 0 1 3 0\n"; // P4=P2
      break;
    case 3:
      geo <<"0 0 0 1 3 1\n"; // P2=P1
      break;
    case 4:
      geo <<"3 0 0 1 0 1\n"; // P4=P3
      break;
    default:
      geo <<"3 0 0 1 3 1\n";
      break;
    }
    REQUIRE(this->read(geo));
  }
};

}


TEST_CASE("TestASMs2D.BoundaryElements")
{
  ASM2DTests<ASMSquare>::BoundaryElements();
}


TEST_CASE("TestASMs2D.Collapse")
{
  for (int iedge = 1; iedge <= 4; iedge++)
  {
    ASMbase::resetNumbering();
    ASMdegenerate2D pch(iedge);
    REQUIRE(pch.uniformRefine(0,2));
    REQUIRE(pch.uniformRefine(1,1));
    REQUIRE(pch.generateFEMTopology());
    std::cout <<"Degenerating E"<< iedge << std::endl;
#ifdef SP_DEBUG
    pch.write(std::cout,0);
#endif
    REQUIRE(pch.collapseEdge(iedge));
  }
}


TEST_CASE("TestASMs2D.Connect")
{
  ASM2DTests<ASMSquare>::Connect(ASM::Spline);
}


TEST_CASE("TestASMs2D.ConstrainEdge")
{
  ASM2DTests<ASMSquare>::ConstrainEdge();
}


TEST_CASE("TestASMs2D.ConstrainEdgeOpen")
{
  ASM2DTests<ASMSquare>::ConstrainEdgeOpen();
}


TEST_CASE("TestASMs2D.ElementConnectivities")
{
  const auto ref = std::array{
    std::array{-1,  1, -1,  2},
    std::array{ 0, -1, -1,  3},
    std::array{-1,  3,  0, -1},
    std::array{ 2, -1,  1, -1}
  };
  ASM2DTests<ASMSquare>::GetElementConnectivities(ref);
}


TEST_CASE("TestASMs2D.ElmNodes")
{
  const auto ref = std::array{
    std::array{0,1,3,4},
    std::array{1,2,4,5},
    std::array{3,4,6,7},
    std::array{4,5,7,8},
  };
  const auto ref_proj = std::array{
    std::array{0,1,2,4,5,6,8,9,10},
    std::array{1,2,3,5,6,7,9,10,11},
    std::array{4,5,6,8,9,10,12,13,14},
    std::array{5,6,7,9,10,11,13,14,15},
  };

  ASM2DTests<ASMSquare>::ElmNodes(ref, ref_proj);
}


TEST_CASE("TestASMs2D.GetElementCorners")
{
  ASM2DTests<ASMTrap<ASMs2D>>::GetElementCorners();
}


TEST_CASE("TestASMs2D.Write")
{
  ASM2DTests<ASMSquare>::Write(false);
}

namespace {

//! \brief Structured patch giving access to its thread groups.
class ASMs2DThreads : public ASMSquare
{
public:
  ASMs2DThreads() : ASMSquare(1) {}

  //! \brief Generates the thread groups from tiles of the quadratic degree.
  void genThreadGroups() { this->generateTileGroups(2, 2, true); }

  //! \brief Returns the thread groups.
  const ThreadGroups& getThreadGroups() const { return threadGroups; }

  //! \brief Returns the parity coloring of the tiles, ignoring constraints.
  ThreadGroups parityGroups() const
  {
    const std::vector<bool> el(6,true);
    IntVec parity;
    const IntMat tiles = ThreadGroups::tiles(el,el,2,2,&parity);
    ThreadGroups groups;
    groups.setColors(tiles,parity);
    return groups;
  }

  //! \brief Returns the number of pairs of tasks of a color sharing a node.
  size_t countConflicts() const { return this->countConflicts(threadGroups); }

  //! \brief Returns the number of pairs of tasks of a color sharing a node.
  size_t countConflicts(const ThreadGroups& groups) const
  {
    const IntMat nodes = this->getElmWriteNodes();
    size_t conflicts = 0;
    for (size_t c = 0; c < groups.size(); c++)
    {
      const IntMat& tasks = groups[c];
      std::vector<std::set<int>> taskNodes(tasks.size());
      for (size_t t = 0; t < tasks.size(); t++)
        for (int iel : tasks[t])
          taskNodes[t].insert(nodes[iel].begin(),nodes[iel].end());
      for (size_t t = 0; t < tasks.size(); t++)
        for (size_t u = t+1; u < tasks.size(); u++)
          for (int node : taskNodes[t])
            if (taskNodes[u].count(node))
            {
              ++conflicts;
              break;
            }
    }
    return conflicts;
  }
};

#ifdef USE_OPENMP
//! \brief Returns a quadratic patch of 6x6 elements, i.e., 3x3 tiles.
std::unique_ptr<ASMs2DThreads> quadraticPatch()
{
  ASMbase::resetNumbering();
  auto pch = std::make_unique<ASMs2DThreads>();
  REQUIRE(pch->raiseOrder(1,1));
  REQUIRE(pch->uniformRefine(0,5));
  REQUIRE(pch->uniformRefine(1,5));
  REQUIRE(pch->generateFEMTopology());
  return pch;
}
#endif

}

#ifdef USE_OPENMP
TEST_CASE("TestASMs2D.ThreadGroupsPeriodic")
{
  const int nThreads = omp_get_max_threads();
  omp_set_num_threads(2);

  // Without constraints, the tiles get the four parity colors
  std::unique_ptr<ASMs2DThreads> pch = quadraticPatch();
  pch->genThreadGroups();
  REQUIRE(pch->getThreadGroups().size() == 4);
  REQUIRE(pch->countConflicts() == 0);

  // The first and last column of tiles have the same parity. Coupling the
  // nodes on the edges between them through a seam makes them neighbors.
  SECTION("MPC") {
    for (int j = 0; j < 8; j++)
      REQUIRE(pch->add2PC(8+8*j, 1, 1+8*j));
  }
  SECTION("Collapsed") {
    pch->closeBoundaries(1, 1, 1);
  }

  // The parity colors now conflict, while the colors of the patch do not
  REQUIRE(pch->countConflicts(pch->parityGroups()) > 0);
  pch->genThreadGroups();
  REQUIRE(pch->getThreadGroups().size() > 4);
  REQUIRE(pch->countConflicts() == 0);

  omp_set_num_threads(nThreads);
}
#endif
