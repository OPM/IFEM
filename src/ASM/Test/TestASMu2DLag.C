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

#include "Catch2Support.h"
#include "ASMu2DLag.h"
#include "MPC.h"
#include "Vec3Oper.h"

#include <fstream>
#include <sstream>


namespace
{
  std::istream& generateXMLModel(std::stringstream& str,
                                 float Lx, float Ly, int nx, int ny,
                                 int oneNodeElms = 0)
  {
    str << "<patch>\n<nodes>\n";
    const double dx = Lx / nx;
    const double dy = Ly / ny;
    for (int j = 0; j <= ny; ++j)
      for (int i = 0; i <= nx; ++i)
        str << i*dx << " " << j*dy << " 0.0\n";
    for (int i = 0; i < oneNodeElms; ++i)
      str << i+2 << " 0 0\n";
    str << "</nodes>\n<elements nenod=\"4\">\n";
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
        str << i  +  j   *(nx+1) << " " << i+1 +  j   *(nx+1) << " "
            << i  + (j+1)*(nx+1) << " " << i+1 + (j+1)*(nx+1) << "\n";
    str << "</elements>\n";
    if (oneNodeElms > 0) {
      str << "<elements nenod=\"1\">\n";
      for (int i = 0; i < oneNodeElms; ++i)
        str << (nx+1)*(ny+1) + i << "\n";
      str << "</elements>\n";
    }
    str << "</patch>\n";
    return str;
  }

  class ASMu2DLagTest : public ASMu2DLag
  {
  public:
    ASMu2DLagTest(float Lx, float Ly, int nx, int ny, int oneNodeElms = 0)
      : ASMu2DLag(2,2,'x')
    {
      ASMbase::resetNumbering();
      std::stringstream myStream;
      this->read(generateXMLModel(myStream,Lx,Ly,nx,ny,oneNodeElms));
    }

    const IntMat& genThreadGroups(bool separateGroup1noded = false)
    {
      this->generateThreadGroupsMultiColored(false,separateGroup1noded);
      return threadGroups[0];
    }
  };
}


TEST_CASE("TestASMu2DLag.GenerateThreadGroups3x3")
{
  ASMu2DLagTest pch(1.0, 1.0, 3, 3);
  REQUIRE(pch.generateFEMTopology());

  const IntMat& group = pch.genThreadGroups();

  REQUIRE(group.size() == 4);
  REQUIRE(group[0].size() == 4);
  REQUIRE(group[1].size() == 2);
  REQUIRE(group[2].size() == 2);
  REQUIRE(group[3].size() == 1);
  REQUIRE(group[0] == std::vector{0, 2, 6, 8});
  REQUIRE(group[1] == std::vector{1, 7});
  REQUIRE(group[2] == std::vector{3, 5,});
  REQUIRE(group[3] == std::vector{4});
}


TEST_CASE("TestASMu2DLag.GenerateThreadGroups4x4")
{
  ASMu2DLagTest pch(1.0, 1.0, 4, 4);
  REQUIRE(pch.generateFEMTopology());

  const IntMat& group = pch.genThreadGroups();

  REQUIRE(group.size() == 4);
  REQUIRE(group[0].size() == 4);
  REQUIRE(group[1].size() == 4);
  REQUIRE(group[2].size() == 4);
  REQUIRE(group[3].size() == 4);
  REQUIRE(group[0] == std::vector{0, 2, 8, 10});
  REQUIRE(group[1] == std::vector{1, 3, 9, 11});
  REQUIRE(group[2] == std::vector{4, 6, 12, 14});
  REQUIRE(group[3] == std::vector{5, 7, 13, 15});
}


TEST_CASE("TestASMu2DLag.GenerateThreadGroups3x3TwoPC")
{
  ASMu2DLagTest pch1(1.0, 1.0, 3, 3);
  pch1.add2PC(1, 1, 16);
  REQUIRE(pch1.generateFEMTopology());

  ASMu2DLagTest pch2(1.0, 1.0, 3, 3);
  pch2.add2PC(16, 1, 1);
  REQUIRE(pch2.generateFEMTopology());

  for (ASMu2DLagTest* pch : { &pch1, &pch2 })
  {
    const IntMat& group = pch->genThreadGroups();

    REQUIRE(group.size() == 5);
    REQUIRE(group[0].size() == 3);
    REQUIRE(group[1].size() == 2);
    REQUIRE(group[2].size() == 2);
    REQUIRE(group[3].size() == 1);
    REQUIRE(group[4].size() == 1);
    REQUIRE(group[0] == std::vector{0, 2, 6});
    REQUIRE(group[1] == std::vector{1, 7});
    REQUIRE(group[2] == std::vector{3, 5,});
    REQUIRE(group[3] == std::vector{4});
    REQUIRE(group[4] == std::vector{8});
  }
}


TEST_CASE("TestASMu2DLag.GenerateThreadGroups3x3MPC")
{
  ASMu2DLagTest pch(1.0, 1.0, 3, 3);
  MPC* mpc = new MPC(1, 1);
  mpc->addMaster(16, 1);
  mpc->addMaster(13, 1);
  pch.addMPC(mpc);
  REQUIRE(pch.generateFEMTopology());

  const IntMat& group = pch.genThreadGroups();

  REQUIRE(group.size() == 4);
  REQUIRE(group[0].size() == 3);
  REQUIRE(group[1].size() == 3);
  REQUIRE(group[2].size() == 2);
  REQUIRE(group[3].size() == 1);
  REQUIRE(group[0] == std::vector{0, 2, 7});
  REQUIRE(group[1] == std::vector{1, 6, 8});
  REQUIRE(group[2] == std::vector{3, 5,});
  REQUIRE(group[3] == std::vector{4});
}


TEST_CASE("TestASMu2DLag.GenerateThreadGroups3x3OneNode")
{
  ASMu2DLagTest pch(1.0, 1.0, 3, 3, 4);
  REQUIRE(pch.generateFEMTopology());

  for (bool with1 : { false, true })
  {
    const IntMat& group = pch.genThreadGroups(with1);

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

    REQUIRE(group.size() == ref.size());
    for (size_t i = 0; i < ref.size(); ++i) {
      REQUIRE(group[i].size() == ref[i].size());
      REQUIRE(group[i] == ref[i]);
    }
  }
}


TEST_CASE("TestASMu2DLag.GenerateThreadGroups3x3OneNodeSPC")
{
  ASMu2DLagTest pch(1.0, 1.0, 3, 3, 4);
  pch.add2PC(17, 1, 1);
  REQUIRE(pch.generateFEMTopology());

  for (bool with1 : { false, true })
  {
    const IntMat& group = pch.genThreadGroups(with1);

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

    REQUIRE(group.size() == ref.size());
    for (size_t i = 0; i < ref.size(); ++i) {
      REQUIRE(group[i].size() == ref[i].size());
      REQUIRE(group[i] == ref[i]);
    }
  }
}


TEST_CASE("TestASMu2DLag.findPoints")
{
  struct TestCase
  {
    const char* fname = nullptr;
    int nX = 0;
    int nY = 0;
    size_t npt = 0;
    Vec3 X;
  };

  const TestCase param =
    GENERATE(TestCase{ "5x4 mesh", 5, 4, 3, Vec3(0.123,0.345,0.0) },
             TestCase{ "src/ASM/Test/refdata/Fine2Dmesh.xinp", 0, 0, 17,
                 Vec3(0.1,0.0,0.0) }
             );

  SECTION(param.fname)
  {
    std::cout <<"\nTesting "<< param.fname;
    if (param.nX > 0) std::cout <<" "<< param.nX <<"x"<< param.nY;
    std::cout <<" "<< param.npt << std::endl;

    ASMbase::resetNumbering();
    ASMu2DLag pch(3,3,'x');
    if (param.nX > 0 && param.nY > 0)
    {
      // Generate a simple mesh on the bi-unit square
      std::stringstream str;
      REQUIRE(pch.read(generateXMLModel(str,1.0,1.0,param.nX,param.nY)));
    }
    else
    {
      // Read mesh from file
      std::ifstream fs(param.fname);
      REQUIRE(pch.read(fs));
    }
    REQUIRE(pch.generateFEMTopology());

    // Pick some random elements and search for their center points
    const size_t nel = pch.getNoElms();
    const size_t npt = param.npt;
    std::vector<size_t> elements; elements.reserve(npt);
    std::vector<Vec3> points; points.reserve(npt);
    for (size_t iel = 1+nel/(2*npt); iel <= nel; iel += nel/npt)
    {
      elements.push_back(iel);
      points.emplace_back(pch.getElementCenter(iel));
      std::cout <<"Element "<< iel <<": "<< points.back() << std::endl;
    }
    points.emplace_back(param.X);

    std::vector<ASMbase::PointParams> elms;
    REQUIRE(pch.findPoints(points,elms));

    size_t ipt = 0;
    std::cout <<"Search results:";
    for (const ASMbase::PointParams& elm : elms)
    {
      std::cout <<"\nPoint "<< ++ipt <<": iel = "<< elm.iel
                <<" xi,eta = "<< elm.u[0] <<","<< elm.u[1]
                <<" dist = "<< elm.dist;
      if (ipt-1 < elements.size())
        REQUIRE(elements[ipt-1] == elm.iel);
      REQUIRE(elm.dist <= 1.0e-12);
    }
    std::cout << std::endl;
  }
}
