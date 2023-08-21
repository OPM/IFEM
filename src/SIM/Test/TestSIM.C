//==============================================================================
//!
//! \file TestSIM.C
//!
//! \date Apr 27 2016
//!
//! \author Arne Morten Kvarving / NTNU
//!
//! \brief Tests for various expected SIM behavior.
//!
//==============================================================================

#include "SIM2D.h"
#include "SIM3D.h"
#include "ASMmxBase.h"
#include "IntegrandBase.h"

#include "gtest/gtest.h"


template<class Dim> class TestProjectSIM : public Dim
{
public:
  TestProjectSIM(const SIMinput::CharVec& nf) : Dim(nf)
  {
    Dim::myProblem = new TestProjectIntegrand(Dim::dimension);
    EXPECT_TRUE(this->createDefaultModel());
    EXPECT_TRUE(this->preprocess());
  }
  virtual ~TestProjectSIM() {}

  bool addMixedMADOF(unsigned char basis, unsigned char nndof)
  {
    return this->addMADOF(basis, nndof);
  }

private:
  class TestProjectIntegrand : public IntegrandBase
  {
  public:
    TestProjectIntegrand(int dim) : IntegrandBase(dim) {}

    using IntegrandBase::evalSol;
    virtual bool evalSol(Vector& s, const FiniteElement&, const Vec3& X,
                         const std::vector<int>&) const
    {
      s.resize(1);
      s(1) = X.sum();
      return true;
    }

    virtual size_t getNoFields(int) const { return 1; }
  };
};


TEST(TestSIM2D, UniqueBoundaryNodes)
{
  const char* boundary_nodes = "<geometry>"
    "<patchfile>src/ASM/Test/refdata/square-4-orient0.g2</patchfile>"
    "<topology>"
    "  <connection master='1' medge='4' slave='2' sedge='3'/>"
    "  <connection master='1' medge='2' slave='3' sedge='1'/>"
    "  <connection master='2' medge='2' slave='4' sedge='1'/>"
    "  <connection master='3' medge='4' slave='4' sedge='3'/>"
    "</topology>"
    "<topologysets>"
    "  <set name='dir' type='edge'>"
    "    <item patch='1'>3</item>"
    "    <item patch='3'>3</item>"
    "  </set>"
    "</topologysets>"
    "</geometry>";

  SIM2D sim(1);
  ASSERT_TRUE(sim.loadXML(boundary_nodes));
  ASSERT_TRUE(sim.createFEMmodel());

  std::vector<int> vec;
  sim.getBoundaryNodes(sim.getUniquePropertyCode("dir",0),vec);

  std::sort(vec.begin(), vec.end());
  ASSERT_TRUE(std::unique(vec.begin(), vec.end()) == vec.end());
}


TEST(TestSIM2D, Periodic)
{
  const char* geometry = "<geometry>"
    "<patchfile>src/SIM/Test/cylinder-shell.g2</patchfile>"
    "<refine patch='1' u='3' v='3'/>"
    "<periodic patch='1' dir='2'/>"
    "</geometry>";

  SIM2D sim(1);
  ASSERT_TRUE(sim.loadXML(geometry));
  ASSERT_TRUE(sim.createFEMmodel());
}


TEST(TestSIM2D, ProjectSolution)
{
  TestProjectSIM<SIM2D> sim({1});

  Matrix ssol;
  ASSERT_TRUE(sim.project(ssol, Vector(sim.getNoDOFs())));

  size_t n = 1;
  for (size_t j = 0; j < 2; ++j)
    for (size_t i = 0; i < 2; ++i)
      EXPECT_FLOAT_EQ(ssol(1, n++), i + j);
}


TEST(TestSIM2D, ProjectSolutionMixed)
{
  TestProjectSIM<SIM2D> sim({1,1});

  Matrix ssol;
  ASSERT_TRUE(sim.project(ssol, Vector(sim.getNoDOFs())));

  size_t n = 1;
  for (size_t j = 0; j < 3; ++j)
    for (size_t i = 0; i < 3; ++i)
      EXPECT_FLOAT_EQ(ssol(1, n++), i/2.0 + j/2.0);
}


TEST(TestSIM3D, ProjectSolution)
{
  TestProjectSIM<SIM3D> sim({1});

  Matrix ssol;
  ASSERT_TRUE(sim.project(ssol, Vector(sim.getNoDOFs())));

  size_t n = 1;
  for (size_t k = 0; k < 2; ++k)
    for (size_t j = 0; j < 2; ++j)
      for (size_t i = 0; i < 2; ++i)
        EXPECT_FLOAT_EQ(ssol(1, n++), i + j + k);
}


TEST(TestSIM3D, ProjectSolutionMixed)
{
  TestProjectSIM<SIM3D> sim({1,1});

  Matrix ssol;
  ASSERT_TRUE(sim.project(ssol, Vector(sim.getNoDOFs())));

  size_t n = 1;
  for (size_t k = 0; k < 3; ++k)
    for (size_t j = 0; j < 3; ++j)
      for (size_t i = 0; i < 3; ++i)
        EXPECT_FLOAT_EQ(ssol(1, n++), i/2.0 + j/2.0 + k/2.0);
}


TEST(TestSIM2D, InjectPatchSolution)
{
  ASMmxBase::Type = ASMmxBase::REDUCED_CONT_RAISE_BASIS1;
  ASMmxBase::itgBasis = 2;
  TestProjectSIM<SIM2D> sim({1,1});
  ASMbase* pch = sim.getPatch(1);

  Vector sol(2*sim.getNoNodes(1) + sim.getNoNodes(2));
  Vector lsol(2*sim.getNoNodes(1));
  size_t i, ofs;
  for (i = 0; i < sim.getNoNodes(1); i++)
    lsol[2*i] = lsol[2*i+1] = i+1;

  ASSERT_TRUE(sim.addMixedMADOF(1, 2));
  sim.injectPatchSolution(sol, lsol, pch, 2, 1);
  for (i = ofs = 0; i < sim.getNoNodes(1); i++, ofs += 2) {
    EXPECT_FLOAT_EQ(sol[ofs], i+1);
    EXPECT_FLOAT_EQ(sol[ofs+1], i+1);
  }
  for (i = 0; i < sim.getNoNodes(2); i++, ofs++)
    EXPECT_FLOAT_EQ(sol[ofs], 0);

  ASSERT_TRUE(sim.addMixedMADOF(2, 2));
  Vector sol2(sim.getNoNodes(1) + 2*sim.getNoNodes(2));
  Vector lsol2(2*sim.getNoNodes(2));
  for (i = 0; i < sim.getNoNodes(2); i++)
    lsol2[2*i] = lsol2[2*i+1] = i+1;

  sim.injectPatchSolution(sol2, lsol2, pch, 2, 2);
  for (i = ofs = 0; i < sim.getNoNodes(1); i++, ofs++)
    EXPECT_FLOAT_EQ(sol2[ofs], 0);

  for (i = 0; i < sim.getNoNodes(2); i++, ofs += 2) {
    EXPECT_FLOAT_EQ(sol2[ofs], i+1);
    EXPECT_FLOAT_EQ(sol2[ofs+1], i+1);
  }
}


class TestSIM2D : public testing::Test,
                  public testing::WithParamInterface<std::pair<int,ASM::Discretization>>
{
};


TEST_P(TestSIM2D, ElmConnectivites)
{
  SIM2D sim;
  std::stringstream str;
  sim.opt.discretization = GetParam().second;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_2D_4_orient";
  str << GetParam().first << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.preprocess());

  auto neighs = sim.getElmConnectivities();

  using IntMat = std::vector<std::vector<int>>;
  const std::vector<IntMat> ref = {{{-1,  1, -1,  2},
                                    { 0,  8, -1,  3},
                                    {-1,  3,  0,  4},
                                    { 2, 10,  1,  5},
                                    {-1,  5,  2,  6},
                                    { 4, 12,  3,  7},
                                    {-1,  7,  4, -1},
                                    { 6, 14,  5, -1},
                                    { 1,  9, -1, 10},
                                    { 8, -1, -1, 11},
                                    { 3, 11,  8, 12},
                                    {10, -1,  9, 13},
                                    { 5, 13, 10, 14},
                                    {12, -1, 11, 15},
                                    { 7, 15, 12, -1},
                                    {14, -1, 13, -1}},

                                   {{-1,  1, -1,  2},
                                    { 0,  9, -1,  3},
                                    {-1,  3,  0,  4},
                                    { 2,  8,  1,  5},
                                    {-1,  5,  2,  6},
                                    { 4, 12,  3,  7},
                                    {-1,  7,  4, -1},
                                    { 6, 14,  5, -1},
                                    {12,  9,  3, 10},
                                    { 8, -1,  1, 11},
                                    {13, 11,  8, -1},
                                    {10, -1,  9, -1},
                                    { 5, 13,  8, 14},
                                    {12, -1, 10, 15},
                                    { 7, 15, 12, -1},
                                    {14, -1, 13, -1}}};

  ASSERT_EQ(neighs.size(), 16U);
  IntMat r = ref[GetParam().first];
  for (size_t e = 0; e < 16; ++e) {
    if (GetParam().second == ASM::LRSpline) {
      std::sort(r[e].begin(), r[e].end());
      auto it = r[e].begin();
      while (*it == -1) ++it;
      r[e].erase(r[e].begin(), it);
      std::sort(neighs[e].begin(), neighs[e].end());
    }
    ASSERT_EQ(neighs[e].size(), r[e].size());
    for (size_t n = 0; n < neighs[e].size(); ++n) {
      EXPECT_EQ(neighs[e][n], r[e][n]);
    }
  }
}


std::vector<std::pair<int,ASM::Discretization>> orientations2D = {{0, ASM::Spline},
                                                                  {1, ASM::Spline}
#ifdef HAS_LRSPLINE
                                                                 ,{0, ASM::LRSpline},
                                                                  {1, ASM::LRSpline}
#endif
                                                                 };
INSTANTIATE_TEST_SUITE_P(TestSIM2D, TestSIM2D, testing::ValuesIn(orientations2D));


class TestSIM3D : public testing::Test,
                  public testing::WithParamInterface<std::pair<int,ASM::Discretization>>
{
};


TEST_P(TestSIM3D, ElmConnectivities)
{
  SIM3D sim;
  std::stringstream str;
  sim.opt.discretization = GetParam().second;
  str << "src/ASM/Test/refdata/DomainDecomposition_MPI_3D_4_orient";
  str << GetParam().first << ".xinp";
  ASSERT_TRUE(sim.read(str.str().c_str()));
  ASSERT_TRUE(sim.preprocess());

  auto neighs = sim.getElmConnectivities();

  using IntMat = std::vector<std::vector<int>>;
  static std::vector<IntMat> ref = {{{-1,  1, -1,  2, -1,  4},
                                     { 0,  8, -1,  3, -1,  5},
                                     {-1,  3,  0, 16, -1,  6},
                                     { 2, 10,  1, 17, -1,  7},
                                     {-1,  5, -1,  6,  0, 32},
                                     { 4, 12, -1,  7,  1, 33},
                                     {-1,  7,  4, 20,  2, 34},
                                     { 6, 14,  5, 21,  3, 35},
                                     { 1,  9, -1, 10, -1, 12},
                                     { 8, -1, -1, 11, -1, 13},
                                     { 3, 11,  8, 24, -1, 14},
                                     {10, -1,  9, 25, -1, 15},
                                     { 5, 13, -1, 14,  8, 40},
                                     {12, -1, -1, 15,  9, 41},
                                     { 7, 15, 12, 28, 10, 42},
                                     {14, -1, 13, 29, 11, 43},
                                     {-1, 17,  2, 18, -1, 20},
                                     {16, 24,  3, 19, -1, 21},
                                     {-1, 19, 16, -1, -1, 22},
                                     {18, 26, 17, -1, -1, 23},
                                     {-1, 21,  6, 22, 16, 48},
                                     {20, 28,  7, 23, 17, 49},
                                     {-1, 23, 20, -1, 18, 50},
                                     {22, 30, 21, -1, 19, 51},
                                     {17, 25, 10, 26, -1, 28},
                                     {24, -1, 11, 27, -1, 29},
                                     {19, 27, 24, -1, -1, 30},
                                     {26, -1, 25, -1, -1, 31},
                                     {21, 29, 14, 30, 24, 56},
                                     {28, -1, 15, 31, 25, 57},
                                     {23, 31, 28, -1, 26, 58},
                                     {30, -1, 29, -1, 27, 59},
                                     {-1, 33, -1, 34,  4, 36},
                                     {32, 40, -1, 35,  5, 37},
                                     {-1, 35, 32, 48,  6, 38},
                                     {34, 42, 33, 49,  7, 39},
                                     {-1, 37, -1, 38, 32, -1},
                                     {36, 44, -1, 39, 33, -1},
                                     {-1, 39, 36, 52, 34, -1},
                                     {38, 46, 37, 53, 35, -1},
                                     {33, 41, -1, 42, 12, 44},
                                     {40, -1, -1, 43, 13, 45},
                                     {35, 43, 40, 56, 14, 46},
                                     {42, -1, 41, 57, 15, 47},
                                     {37, 45, -1, 46, 40, -1},
                                     {44, -1, -1, 47, 41, -1},
                                     {39, 47, 44, 60, 42, -1},
                                     {46, -1, 45, 61, 43, -1},
                                     {-1, 49, 34, 50, 20, 52},
                                     {48, 56, 35, 51, 21, 53},
                                     {-1, 51, 48, -1, 22, 54},
                                     {50, 58, 49, -1, 23, 55},
                                     {-1, 53, 38, 54, 48, -1},
                                     {52, 60, 39, 55, 49, -1},
                                     {-1, 55, 52, -1, 50, -1},
                                     {54, 62, 53, -1, 51, -1},
                                     {49, 57, 42, 58, 28, 60},
                                     {56, -1, 43, 59, 29, 61},
                                     {51, 59, 56, -1, 30, 62},
                                     {58, -1, 57, -1, 31, 63},
                                     {53, 61, 46, 62, 56, -1},
                                     {60, -1, 47, 63, 57, -1},
                                     {55, 63, 60, -1, 58, -1},
                                     {62, -1, 61, -1, 59, -1}},

                                    {{-1,  1, -1,  2, -1,  4}, // 0
                                     { 0, 14, -1,  3, -1,  5}, // 1
                                     {-1,  3,  0, 16, -1,  6}, // 2
                                     { 2, 12,  1, 17, -1,  7}, // 3
                                     {-1,  5, -1,  6,  0, 32}, // 4
                                     { 4, 10, -1,  7,  1, 33}, // 5
                                     {-1,  7,  4, 20,  2, 34}, // 6
                                     { 6,  8,  5, 21,  3, 35}, // 7
                                     { 7,  9, 28, 10, 42, 12}, // 8
                                     { 8, -1, 29, 11, 43, 13}, // 9
                                     { 5, 11,  8, -1, 40, 14}, // 10
                                     {10, -1,  9, -1, 41, 15}, // 11
                                     { 3, 13, 24, 14,  8, -1}, // 12
                                     {12, -1, 25, 15,  9, -1}, // 13
                                     { 1, 15, 12, -1, 10, -1}, // 14
                                     {14, -1, 13, -1, 11, -1}, // 15
                                     {-1, 17,  2, 18, -1, 20}, // 16
                                     {16, 24,  3, 19, -1, 21}, // 17
                                     {-1, 19, 16, -1, -1, 22}, // 18
                                     {18, 26, 17, -1, -1, 23}, // 19
                                     {-1, 21,  6, 22, 16, 48}, // 20
                                     {20, 28,  7, 23, 17, 49}, // 21
                                     {-1, 23, 20, -1, 18, 50}, // 22
                                     {22, 30, 21, -1, 19, 51}, // 23
                                     {17, 25, 12, 26, -1, 28}, // 24
                                     {24, -1, 13, 27, -1, 29}, // 25

                                     {19, 27, 24, -1, -1, 30}, // 26
                                     {26, -1, 25, -1, -1, 31}, // 27
                                     {21, 29,  8, 30, 24, 56}, // 28
                                     {28, -1,  9, 31, 25, 57}, // 29
                                     {23, 31, 28, -1, 26, 58}, // 30
                                     {30, -1, 29, -1, 27, 59}, // 31
                                     {-1, 33, -1, 34,  4, 36}, // 32
                                     {32, 40, -1, 35,  5, 37}, // 33
                                     {-1, 35, 32, 48,  6, 38}, // 34
                                     {34, 42, 33, 49,  7, 39}, // 35
                                     {-1, 37, -1, 38, 32, -1}, // 36
                                     {36, 44, -1, 39, 33, -1}, // 37
                                     {-1, 39, 36, 52, 34, -1}, // 38
                                     {38, 46, 37, 53, 35, -1}, // 39
                                     {33, 41, -1, 42, 10 ,44}, // 40
                                     {40, -1, -1, 43, 11, 45}, // 41
                                     {35, 43, 40, 56,  8, 46}, // 42
                                     {42, -1, 41, 57,  9, 47}, // 43
                                     {37, 45, -1, 46, 40, -1}, // 44
                                     {44, -1, -1, 47, 41, -1}, // 45
                                     {39, 47, 44, 60, 42, -1}, // 46
                                     {46, -1, 45, 61, 43, -1}, // 47
                                     {-1, 49, 34, 50, 20, 52}, // 48
                                     {48, 56, 35, 51, 21, 53}, // 49
                                     {-1, 51, 48, -1, 22, 54}, // 50
                                     {50, 58, 49, -1, 23, 55}, // 51
                                     {-1, 53, 38, 54, 48, -1}, // 52
                                     {52, 60, 39, 55, 49, -1}, // 53
                                     {-1, 55, 52, -1, 50, -1}, // 54
                                     {54, 62, 53, -1, 51, -1}, // 55
                                     {49, 57, 42, 58, 28, 60}, // 56
                                     {56, -1, 43, 59, 29, 61}, // 57
                                     {51, 59, 56, -1, 30, 62}, // 58
                                     {58, -1, 57, -1, 31, 63}, // 59
                                     {53, 61, 46, 62, 56, -1}, // 60
                                     {60, -1, 47, 63, 57, -1}, // 61
                                     {55, 63, 60, -1, 58, -1}, // 62
                                     {62, -1, 61, -1, 59, -1}}}; // 63

  ASSERT_EQ(neighs.size(), 64U);
  IntMat r = ref[GetParam().first];
  for (size_t e = 0; e < 64; ++e) {
    if (GetParam().second == ASM::LRSpline) {
      std::sort(r[e].begin(), r[e].end());
      auto it = r[e].begin();
      while (*it == -1) ++it;
      r[e].erase(r[e].begin(), it);
      std::sort(neighs[e].begin(), neighs[e].end());
    }
    ASSERT_EQ(neighs[e].size(), r[e].size());
    for (size_t n = 0; n < neighs[e].size(); ++n) {
      EXPECT_EQ(neighs[e][n], r[e][n]);
    }
  }
}


const std::vector<std::pair<int,ASM::Discretization>> orientations3D = {{0, ASM::Spline},
                                                                        {1, ASM::Spline}
#ifdef HAS_LRSPLINE
                                                                       ,{0, ASM::LRSpline},
                                                                        {1, ASM::LRSpline}
#endif
                                                                       };

INSTANTIATE_TEST_SUITE_P(TestSIM3D, TestSIM3D, testing::ValuesIn(orientations3D));
