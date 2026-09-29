//==============================================================================
//!
//! \file TestASMu2DLag.C
//!
//! \date Oct 22 2025
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for assembly of unstructured 2D %Lagrange FE models.
//!
//==============================================================================

#include "ASMu2DLag.h"

#include "Catch2Support.h"
#include "MPC.h"

#include <sstream>

namespace
{

class ASMu2DLagTest : public ASMu2DLag
{
public:
  ASMu2DLagTest() : ASMu2DLag(2,2,'x')
  {
    ASMbase::resetNumbering();
    ASM::coloring = ASM::FIRST_FIT; // the reference groups are greedy ones
  }

  void genThreadGroups(bool separateGroup1noded = false)
  {
    this->generateThreadGroupsMultiColored(true, separateGroup1noded);
  }
  //! \brief Returns the elements of each color.
  std::vector<std::vector<int>> getColors() const
  {
    std::vector<std::vector<int>> colors(threadGroups.size());
    for (size_t c = 0; c < threadGroups.size(); c++)
      for (const std::vector<int>& task : threadGroups[c])
        colors[c].insert(colors[c].end(),task.begin(),task.end());
    return colors;
  }
  bool collapse(int node1, int node2)
  {
    return ASMbase::collapseNodes(*this,node1,*this,node2);
  }
};


void generateXMLModel(std::stringstream& str,
                      float Lx, float Ly, int nx, int ny, int oneNodeElms = 0)
{
  str << "<patch>\n<nodes>\n";
  const double dx = Lx / nx;
  const double dy = Ly / ny;
  for (int j = 0; j < ny+1; ++j)
    for (int i = 0; i < nx+1; ++i)
      str << i*dx << " " << j*dy << " 0.0\n";
  for (int i = 0; i < oneNodeElms; ++i)
    str << i+2 << " 0 0\n";
  str << "</nodes>\n<elements nenod=\"4\">\n";
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      str << i + j*(nx+1) << " " << i + 1 + j*(nx+1) << " "
            << i + 1 + (j+1)*(nx+1) << " " << i + (j+1)*(nx+1) << '\n';
  str << "</elements>\n";
  if (oneNodeElms > 0) {
    str << "<elements nenod=\"1\">\n";
    for (int i = 0; i < oneNodeElms; ++i)
      str << i+(nx+1)*(ny+1) << '\n';
    str << "</elements>\n";
  }
  str << "</patch>\n";
}

}


TEST_CASE("TestASMu2DLag.GenerateThreadGroups3x3")
{
  std::stringstream str;
  generateXMLModel(str, 1.0, 1.0, 3, 3);

  ASMu2DLagTest pch;
  REQUIRE(pch.read(str));
  REQUIRE(pch.generateFEMTopology());

  pch.genThreadGroups();
  const std::vector<std::vector<int>> groups = pch.getColors();
  REQUIRE(groups.size() == 4);
  REQUIRE(groups[0].size() == 4);
  REQUIRE(groups[1].size() == 2);
  REQUIRE(groups[2].size() == 2);
  REQUIRE(groups[3].size() == 1);
  REQUIRE(groups[0] == std::vector{0, 2, 6, 8});
  REQUIRE(groups[1] == std::vector{1, 7});
  REQUIRE(groups[2] == std::vector{3, 5,});
  REQUIRE(groups[3] == std::vector{4});
}


TEST_CASE("TestASMu2DLag.GenerateThreadGroups4x4")
{
  std::stringstream str;
  generateXMLModel(str, 1.0, 1.0, 4, 4);

  ASMu2DLagTest pch;
  REQUIRE(pch.read(str));
  REQUIRE(pch.generateFEMTopology());

  pch.genThreadGroups();
  const std::vector<std::vector<int>> groups = pch.getColors();
  REQUIRE(groups.size() == 4);
  REQUIRE(groups[0].size() == 4);
  REQUIRE(groups[1].size() == 4);
  REQUIRE(groups[2].size() == 4);
  REQUIRE(groups[3].size() == 4);
  REQUIRE(groups[0] == std::vector{0, 2, 8, 10});
  REQUIRE(groups[1] == std::vector{1, 3, 9, 11});
  REQUIRE(groups[2] == std::vector{4, 6, 12, 14});
  REQUIRE(groups[3] == std::vector{5, 7, 13, 15});
}


TEST_CASE("TestASMu2DLag.GenerateThreadGroups3x3TwoPC")
{
  std::stringstream str;
  generateXMLModel(str, 1.0, 1.0, 3, 3);

  ASMu2DLagTest pch;
  REQUIRE(pch.read(str));
  pch.add2PC(1, 1, 16);
  REQUIRE(pch.generateFEMTopology());

  str.clear();
  str.seekg(0, std::ios_base::beg);
  ASMu2DLagTest pch2;
  REQUIRE(pch2.read(str));
  pch2.add2PC(16, 1, 1);
  REQUIRE(pch2.generateFEMTopology());

  auto checks = [](ASMu2DLagTest& p)
  {
    p.genThreadGroups();
    const std::vector<std::vector<int>> groups = p.getColors();
    REQUIRE(groups.size() == 5);
    REQUIRE(groups[0].size() == 3);
    REQUIRE(groups[1].size() == 2);
    REQUIRE(groups[2].size() == 2);
    REQUIRE(groups[3].size() == 1);
    REQUIRE(groups[4].size() == 1);
    REQUIRE(groups[0] == std::vector{0, 2, 6});
    REQUIRE(groups[1] == std::vector{1, 7});
    REQUIRE(groups[2] == std::vector{3, 5,});
    REQUIRE(groups[3] == std::vector{4});
    REQUIRE(groups[4] == std::vector{8});
  };

  checks(pch);
  checks(pch2);
}


TEST_CASE("TestASMu2DLag.GenerateThreadGroups3x3MPC")
{
  std::stringstream str;
  generateXMLModel(str, 1.0, 1.0, 3, 3);

  ASMu2DLagTest pch;
  REQUIRE(pch.read(str));
  MPC* mpc = new MPC(1, 1);
  mpc->addMaster(16, 1);
  mpc->addMaster(13, 1);
  pch.addMPC(mpc);
  REQUIRE(pch.generateFEMTopology());

  pch.genThreadGroups();

  const std::vector<std::vector<int>> groups = pch.getColors();
  REQUIRE(groups.size() == 4);
  REQUIRE(groups[0].size() == 3);
  REQUIRE(groups[1].size() == 3);
  REQUIRE(groups[2].size() == 2);
  REQUIRE(groups[3].size() == 1);
  REQUIRE(groups[0] == std::vector{0, 2, 7});
  REQUIRE(groups[1] == std::vector{1, 6, 8});
  REQUIRE(groups[2] == std::vector{3, 5,});
  REQUIRE(groups[3] == std::vector{4});
}


TEST_CASE("TestASMu2DLag.GenerateThreadGroups3x3OneNode")
{
  std::stringstream str;
  generateXMLModel(str, 1.0, 1.0, 3, 3, 4);

  ASMu2DLagTest pch;
  REQUIRE(pch.read(str));
  REQUIRE(pch.generateFEMTopology());

  auto checks = [](ASMu2DLagTest& p, bool with1)
  {
    p.genThreadGroups(with1);
    const std::vector<std::vector<int>> groups = p.getColors();
    const auto ref =
          with1 ? std::vector{
                    std::vector{9, 10, 11, 12},
                    std::vector{0, 2, 6, 8},
                    std::vector{1, 7},
                    std::vector{3, 5},
                    std::vector{4}
                  }
                :
                  std::vector{
                    std::vector{0, 2, 6, 8, 9, 10, 11, 12},
                    std::vector{1, 7},
                    std::vector{3, 5},
                    std::vector{4}
                  };
    REQUIRE(groups.size() == ref.size());
    for (size_t i = 0; i < ref.size(); ++i) {
      REQUIRE(groups[i].size() == ref[i].size());
      REQUIRE(groups[i] == ref[i]);
    }
  };

  checks(pch, false);
  checks(pch, true);
}


TEST_CASE("TestASMu2DLag.GenerateThreadGroups3x3OneNodeSPC")
{
  std::stringstream str;
  generateXMLModel(str, 1.0, 1.0, 3, 3, 4);

  ASMu2DLagTest pch;
  REQUIRE(pch.read(str));
  pch.add2PC(17, 1, 1);
  REQUIRE(pch.generateFEMTopology());

  auto checks = [](ASMu2DLagTest& p, bool with1)
  {
    p.genThreadGroups(with1);
    const std::vector<std::vector<int>> groups = p.getColors();
    const auto ref =
          with1 ? std::vector{
                    std::vector{10, 11, 12},
                    std::vector{0, 2, 6, 8},
                    std::vector{1, 7, 9},
                    std::vector{3, 5},
                    std::vector{4}
                  }
                :
                  std::vector{
                    std::vector{0, 2, 6, 8, 10, 11, 12},
                    std::vector{1, 7, 9},
                    std::vector{3, 5},
                    std::vector{4}
                  };
    REQUIRE(groups.size() == ref.size());
    for (size_t i = 0; i < ref.size(); ++i) {
      REQUIRE(groups[i].size() == ref[i].size());
      REQUIRE(groups[i] == ref[i]);
    }
  };

  checks(pch, false);
  checks(pch, true);
}


TEST_CASE("TestASMu2DLag.GenerateThreadGroups3x3RemoteMaster")
{
  std::stringstream str;
  generateXMLModel(str, 1.0, 1.0, 3, 3);

  // The corner nodes 1 and 16 are both coupled to node 100 of another patch,
  // so the corner elements 0 and 8 must not be assembled concurrently
  ASMu2DLagTest pch;
  REQUIRE(pch.read(str));
  pch.add2PC(1, 1, 100);
  pch.add2PC(16, 1, 100);
  REQUIRE(pch.generateFEMTopology());

  pch.genThreadGroups();
  const std::vector<std::vector<int>> groups = pch.getColors();
  REQUIRE(groups.size() == 5);
  REQUIRE(groups[0] == std::vector{0, 2, 6});
  REQUIRE(groups[1] == std::vector{1, 7});
  REQUIRE(groups[2] == std::vector{3, 5});
  REQUIRE(groups[3] == std::vector{4});
  REQUIRE(groups[4] == std::vector{8});
}


TEST_CASE("TestASMu2DLag.GenerateThreadGroups3x3Collapsed")
{
  std::stringstream str;
  generateXMLModel(str, 1.0, 1.0, 3, 3);

  // The corner nodes 1 and 16 get a common global node number,
  // so the corner elements 0 and 8 share an equation
  ASMu2DLagTest pch;
  REQUIRE(pch.read(str));
  REQUIRE(pch.generateFEMTopology());
  REQUIRE(pch.collapse(1, 16));

  pch.genThreadGroups();
  const std::vector<std::vector<int>> groups = pch.getColors();
  REQUIRE(groups.size() == 5);
  REQUIRE(groups[0] == std::vector{0, 2, 6});
  REQUIRE(groups[1] == std::vector{1, 7});
  REQUIRE(groups[2] == std::vector{3, 5});
  REQUIRE(groups[3] == std::vector{4});
  REQUIRE(groups[4] == std::vector{8});
}
