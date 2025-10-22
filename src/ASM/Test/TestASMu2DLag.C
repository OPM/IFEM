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
  ASMu2DLagTest() : ASMu2DLag(2,2,'x') { ASMbase::resetNumbering(); }

  void genThreadGroups(bool separateGroup1noded = false)
  {
    this->generateThreadGroupsMultiColored(true, separateGroup1noded);
  }
  const ThreadGroups& getThreadGroups() const { return threadGroups; }
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
  const ThreadGroups& groups = pch.getThreadGroups();
  REQUIRE(groups[0].size() == 4);
  REQUIRE(groups[0][0].size() == 4);
  REQUIRE(groups[0][1].size() == 2);
  REQUIRE(groups[0][2].size() == 2);
  REQUIRE(groups[0][3].size() == 1);
  REQUIRE(groups[0][0] == std::vector{0, 2, 6, 8});
  REQUIRE(groups[0][1] == std::vector{1, 7});
  REQUIRE(groups[0][2] == std::vector{3, 5,});
  REQUIRE(groups[0][3] == std::vector{4});
}


TEST_CASE("TestASMu2DLag.GenerateThreadGroups4x4")
{
  std::stringstream str;
  generateXMLModel(str, 1.0, 1.0, 4, 4);

  ASMu2DLagTest pch;
  REQUIRE(pch.read(str));
  REQUIRE(pch.generateFEMTopology());

  pch.genThreadGroups();
  const ThreadGroups& groups = pch.getThreadGroups();
  REQUIRE(groups[0].size() == 4);
  REQUIRE(groups[0][0].size() == 4);
  REQUIRE(groups[0][1].size() == 4);
  REQUIRE(groups[0][2].size() == 4);
  REQUIRE(groups[0][3].size() == 4);
  REQUIRE(groups[0][0] == std::vector{0, 2, 8, 10});
  REQUIRE(groups[0][1] == std::vector{1, 3, 9, 11});
  REQUIRE(groups[0][2] == std::vector{4, 6, 12, 14});
  REQUIRE(groups[0][3] == std::vector{5, 7, 13, 15});
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
    const ThreadGroups& groups = p.getThreadGroups();
    REQUIRE(groups[0].size() == 5);
    REQUIRE(groups[0][0].size() == 3);
    REQUIRE(groups[0][1].size() == 2);
    REQUIRE(groups[0][2].size() == 2);
    REQUIRE(groups[0][3].size() == 1);
    REQUIRE(groups[0][4].size() == 1);
    REQUIRE(groups[0][0] == std::vector{0, 2, 6});
    REQUIRE(groups[0][1] == std::vector{1, 7});
    REQUIRE(groups[0][2] == std::vector{3, 5,});
    REQUIRE(groups[0][3] == std::vector{4});
    REQUIRE(groups[0][4] == std::vector{8});
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

  const ThreadGroups& groups = pch.getThreadGroups();
  REQUIRE(groups[0].size() == 4);
  REQUIRE(groups[0][0].size() == 3);
  REQUIRE(groups[0][1].size() == 3);
  REQUIRE(groups[0][2].size() == 2);
  REQUIRE(groups[0][3].size() == 1);
  REQUIRE(groups[0][0] == std::vector{0, 2, 7});
  REQUIRE(groups[0][1] == std::vector{1, 6, 8});
  REQUIRE(groups[0][2] == std::vector{3, 5,});
  REQUIRE(groups[0][3] == std::vector{4});
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
    const ThreadGroups& groups = p.getThreadGroups();
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
    REQUIRE(groups[0].size() == ref.size());
    for (size_t i = 0; i < ref.size(); ++i) {
      REQUIRE(groups[0][i].size() == ref[i].size());
      REQUIRE(groups[0][i] == ref[i]);
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
    const ThreadGroups& groups = p.getThreadGroups();
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
    REQUIRE(groups[0].size() == ref.size());
    for (size_t i = 0; i < ref.size(); ++i) {
      REQUIRE(groups[0][i].size() == ref[i].size());
      REQUIRE(groups[0][i] == ref[i]);
    }
  };

  checks(pch, false);
  checks(pch, true);
}
