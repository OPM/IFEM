//==============================================================================
//!
//! \file ASM2DmxTests.h
//!
//! \date Sep 30 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for 2D mixed FE models.
//!
//==============================================================================

#include "ASMbase.h"
#include "ASMmxBase.h"

#include "Catch2Support.h"

#include <numeric>


template<class Patch>
class ASM2DmxTests
{
public:
  static void AddInterfaceElms()
  {
    ASMmxBase::Type = ASMmxBase::DIV_COMPATIBLE;
    ASMbase::resetNumbering();
    Patch pch({1,1,1});
    pch.raiseOrder(1,1);
    REQUIRE(pch.generateFEMTopology());

    const size_t nel = pch.getNoElms(true,true);
    const size_t nnod = pch.getNoNodes(0); // every basis, no extras yet
    const std::vector<size_t> sizes = pch.getElementSizes(1);
    REQUIRE(sizes.size() == 3);

    // The elements of edge 2, and a made-up neighbour for each of them whose
    // nodes are numbered well clear of this patch's own
    IntVec elms;
    pch.getBoundaryElms(2,elms,-1,true);
    REQUIRE(!elms.empty());

    const size_t nen = std::accumulate(sizes.begin(),sizes.end(),0u);
    std::vector<IntVec> xnodes;
    int next = 1000;
    for (size_t e = 0; e < elms.size(); e++)
    {
      IntVec nodes(nen);
      std::iota(nodes.begin(),nodes.end(),next);
      next += nen;
      xnodes.push_back(nodes);
    }

    std::vector<std::vector<size_t>> xsizes(elms.size(),sizes);

    int nextId = 1000000;
    REQUIRE(pch.addInterfaceElms(2,xnodes,xsizes,nextId));

    // One element for each of the boundary, appended after the ordinary ones
    REQUIRE(pch.getNoElms(true,true) == nel + elms.size());

    for (size_t e = 0; e < elms.size(); e++)
    {
      const IntVec& mine = pch.getElementNodes(elms[e]+1);
      const IntVec& both = pch.getElementNodes(nel+e+1);
      REQUIRE(both.size() == mine.size() + xnodes[e].size());

      // One basis at a time, this patch's nodes and then the neighbour's,
      // which is the order the element matrix of a mixed element is laid out
      // in. Getting this wrong puts the contributions of one basis onto the
      // degrees of freedom of another and is not visible in the answer.
      size_t c = 0, mofs = 0, sofs = 0;
      for (size_t b = 0; b < sizes.size(); b++)
      {
        for (size_t i = 0; i < sizes[b]; i++, c++)
          REQUIRE(both[c] == mine[mofs+i]);
        mofs += sizes[b];
        for (size_t i = 0; i < sizes[b]; i++, c++)
          REQUIRE(pch.getNodeID(both[c]+1) == xnodes[e][sofs+i]);
        sofs += sizes[b];
      }
    }

    // The neighbour's nodes keep the numbers they have where they belong, and
    // one shared by two of the elements is taken in once
    REQUIRE(pch.hasXNodes());
    REQUIRE(nextId == 1000000 + static_cast<int>(elms.size()));

    // And the patch says which of the nodes it holds are not its own. What
    // such a node carries is the owner's to say: this patch would answer
    // with its own first basis, which is the right answer only where every
    // basis has the same number of fields, and the pressure basis of a
    // velocity-only model has none at all.
    for (size_t i = 1; i <= pch.getNoNodes(0); i++)
      REQUIRE(pch.isBorrowedNode(i) == (i > nnod));
  }

};
