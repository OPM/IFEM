//==============================================================================
//!
//! \file TestMultiPatchModelGenerator.C
//!
//! \date Sep 2 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for multi-patch model generators.
//!
//==============================================================================

#include "MultiPatchModelGenerator.h"
#include "SIMMultiPatchModelGen.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "ASMs1D.h"
#include "ASMs2D.h"
#include "ASMs3D.h"
#include "Functions.h"
#include "SplineUtils.h"
#include "Utilities.h"

#include "gtest/gtest.h"
#include "tinyxml.h"

#include <GoTools/geometry/SplineCurve.h>
#include <GoTools/geometry/SplineSurface.h>
#include <GoTools/utils/Point.h>


auto&& check_vector_int_equals_range = [](const std::vector<int>& arr,
                                          const std::vector<int>& range)
{
  int el = range[0] - 1;
  size_t len = range[1] - el;
  ASSERT_EQ(arr.size(), len);
  for (int i : arr)
    EXPECT_EQ(i, ++el);
};


auto&& check_vector_int_equals = [](const std::vector<int>& arr1,
                                    const std::vector<int>& arr2)
{
  ASSERT_EQ(arr1.size(), arr2.size());
  auto it2 = arr2.begin();
  for (int i : arr1) {
    EXPECT_EQ(i, *it2);
    it2++;
  }
};


auto&& check_vector_double_near = [](const std::vector<double>& arr1,
                                     const std::vector<double>& arr2,
                                     const double tol=1e-12)
{
  ASSERT_EQ(arr1.size(), arr2.size());
  auto it2 = arr2.begin();
  for (double d : arr1) {
    EXPECT_NEAR(d, *it2, tol);
    it2++;
  }
};


auto&& check_point_near = [](const Go::Point p1,
                             const Go::Point p2,
                             const double tol=1e-12)
{
  ASSERT_EQ(p1.dimension(), p2.dimension());
  auto it2 = p2.begin();
  for (double d : p1) {
    EXPECT_NEAR(d, *it2, tol);
    it2++;
  }
};


template<class Generator>
class TestModelGeneratorWrapper : public Generator {
public:
  TestModelGeneratorWrapper(const TiXmlElement* geo) : Generator(geo) {}
  virtual std::string createG2(int nsd, bool = false) const
  {
    bool rational = false;
    utl::getAttribute(Generator::geo,"rational",rational);
    return this->Generator::createG2(nsd,rational);
  }
};

struct GeomTest {
  std::string xml;
  int dim;
  std::string g2;
  std::string sets;
};


class TestMultiPatchModelGenerator1D :
  public testing::Test,
  public testing::WithParamInterface<GeomTest>
{
};


class TestMultiPatchModelGenerator2D :
  public testing::Test,
  public testing::WithParamInterface<GeomTest>
{
};


class TestMultiPatchModelGenerator3D :
  public testing::Test,
  public testing::WithParamInterface<GeomTest>
{
};


auto&& DoTest = [](const GeomTest& ref, const std::string& gen,
                   const TopologySet& sets)
{
  EXPECT_STREQ(gen.c_str(), ref.g2.c_str());

  if (!ref.sets.empty()) {
    std::string gsets;
    for (auto& it : sets) {
      gsets += it.first + ": ";
      for (auto& it2 : it.second) {
        std::stringstream str;
        str << it2.patch << " " << it2.item << " " << it2.idim << " ";
        gsets += str.str();
      }
      gsets += "\n";
    }
    EXPECT_STREQ(gsets.c_str(), ref.sets.c_str());
  }
};



TEST_P(TestMultiPatchModelGenerator2D, Generate)
{
  TiXmlDocument doc;
  doc.Parse(GetParam().xml.c_str());
  TestModelGeneratorWrapper<MultiPatchModelGenerator2D> gen(doc.RootElement());
  std::string g2 = gen.createG2(GetParam().dim);
  SIM2D sim;
  gen.createTopologySets(sim);
  DoTest(GetParam(), g2, sim.getTopology());
}


TEST(TestMultiPatchModelGenerator2D, GenerateLR)
{
  SIMMultiPatchModelGen<SIM2D> sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("refdata/modelgen2d_lr.xinp"));
  sim.preprocess();
  ASSERT_EQ(sim.getNoNodes(), 28U);
}


TEST(TestMultiPatchModelGenerator2D, GenerateLRmx)
{
  SIMMultiPatchModelGen<SIM2D> sim({2,1});
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("refdata/modelgen2d_lr.xinp"));
  sim.preprocess();
  ASSERT_EQ(sim.getNoNodes(), 73U);
}


TEST(TestMultiPatchModelGenerator2D, Subdivisions)
{
  SIMMultiPatchModelGen<SIM2D> sim(1);
  ASSERT_TRUE(sim.read("refdata/modelgen2d_subdivision.xinp"));

  // check FEM topology
  std::vector<std::vector<int>> mlgn =
      {{ 1, 2, 3, 4, 5,
         6, 7, 8, 9,10,
        11,12,13,14,15,
        16,17,18,19,20,
        21,22,23,24,25},
       { 4, 5,26,27,
         9,10,28,29,
        14,15,30,31,
        19,20,32,33,
        24,25,34,35},
       {16,17,18,19,20,
        21,22,23,24,25,
        36,37,38,39,40,
        41,42,43,44,45},
       {19,20,32,33,
        24,25,34,35,
        39,40,46,47,
        44,45,48,49}};
  int i = 0, ngnod = 0;
  std::map<int,int> g2l;
  for (ASMbase* pch : sim.getFEModel()) {
    pch->renumberNodes(g2l,ngnod);
    check_vector_int_equals(mlgn[i++],pch->getMyNodeNums());
  }
}


TEST(TestMultiPatchModelGenerator2D, InnerPatches)
{
  TiXmlDocument doc;
  doc.Parse("<geometry nx=\"3\" ny=\"3\" sets=\"true\"/>");
  TestModelGeneratorWrapper<MultiPatchModelGenerator2D> gen(doc.RootElement());
  SIM2D sim;
  gen.createTopologySets(sim);
  ASSERT_EQ(sim.topology("InnerPatches").size(), 1u);
  ASSERT_EQ(sim.topology("InnerPatches").begin()->patch, 5u);
}


TEST_P(TestMultiPatchModelGenerator3D, Generate)
{
  TiXmlDocument doc;
  doc.Parse(GetParam().xml.c_str());
  TestModelGeneratorWrapper<MultiPatchModelGenerator3D> gen(doc.RootElement());
  std::string g2 = gen.createG2(GetParam().dim);
  SIM3D sim;
  gen.createTopologySets(sim);
  DoTest(GetParam(), g2, sim.getTopology());
}


TEST(TestMultiPatchModelGenerator3D, GenerateLR)
{
  SIMMultiPatchModelGen<SIM3D> sim(1);
  sim.opt.discretization = ASM::LRSpline;
  ASSERT_TRUE(sim.read("refdata/modelgen3d_lr.xinp"));
  sim.preprocess();
  ASSERT_EQ(sim.getNoNodes(), 112U);
}


TEST(TestMultiPatchModelGenerator3D, Subdivisions)
{
  SIMMultiPatchModelGen<SIM3D> sim(1);
  ASSERT_TRUE(sim.read("refdata/modelgen3d_subdivision.xinp"));

  // check FEM topology
  std::vector<std::vector<int>> mlgn = {
      {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,
       28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,
       52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,
       76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,
       100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,
       118,119,120,121,122,123,124,125},
      {4,5,126,127,9,10,128,129,14,15,130,131,19,20,132,133,24,25,134,135,29,30,
       136,137,34,35,138,139,39,40,140,141,44,45,142,143,49,50,144,145,54,55,
       146,147,59,60,148,149,64,65,150,151,69,70,152,153,74,75,154,155,79,80,
       156,157,84,85,158,159,89,90,160,161,94,95,162,163,99,100,164,165,104,105,
       166,167,109,110,168,169,114,115,170,171,119,120,172,173,124,125,174,175},
      {16,17,18,19,20,21,22,23,24,25,176,177,178,179,180,181,182,183,184,185,41,
       42,43,44,45,46,47,48,49,50,186,187,188,189,190,191,192,193,194,195,66,67,
       68,69,70,71,72,73,74,75,196,197,198,199,200,201,202,203,204,205,91,92,93,
       94,95,96,97,98,99,100,206,207,208,209,210,211,212,213,214,215,116,117,
       118,119,120,121,122,123,124,125,216,217,218,219,220,221,222,223,224,225},
      {19,20,132,133,24,25,134,135,179,180,226,227,184,185,228,229,44,45,142,
       143,49,50,144,145,189,190,230,231,194,195,232,233,69,70,152,153,74,75,
       154,155,199,200,234,235,204,205,236,237,94,95,162,163,99,100,164,165,209,
       210,238,239,214,215,240,241,119,120,172,173,124,125,174,175,219,220,242,
       243,224,225,244,245},
      {76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,
       100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,
       118,119,120,121,122,123,124,125,246,247,248,249,250,251,252,253,254,255,
       256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,
       274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,
       292,293,294,295},
      {79,80,156,157,84,85,158,159,89,90,160,161,94,95,162,163,99,100,164,165,
       104,105,166,167,109,110,168,169,114,115,170,171,119,120,172,173,124,125,
       174,175,249,250,296,297,254,255,298,299,259,260,300,301,264,265,302,303,
       269,270,304,305,274,275,306,307,279,280,308,309,284,285,310,311,289,290,
       312,313,294,295,314,315},
      {91,92,93,94,95,96,97,98,99,100,206,207,208,209,210,211,212,213,214,215,
       116,117,118,119,120,121,122,123,124,125,216,217,218,219,220,221,222,223,
       224,225,261,262,263,264,265,266,267,268,269,270,316,317,318,319,320,321,
       322,323,324,325,286,287,288,289,290,291,292,293,294,295,326,327,328,329,
       330,331,332,333,334,335},
      {94,95,162,163,99,100,164,165,209,210,238,239,214,215,240,241,119,120,
       172,173,124,125,174,175,219,220,242,243,224,225,244,245,264,265,302,303,
       269,270,304,305,319,320,336,337,324,325,338,339,289,290,312,313,294,295,
       314,315,329,330,340,341,334,335,342,343}};
  int i = 0, ngnod = 0;
  std::map<int,int> g2l;
  for (ASMbase* pch : sim.getFEModel()) {
    pch->renumberNodes(g2l,ngnod);
    check_vector_int_equals(mlgn[i++],pch->getMyNodeNums());
  }
}


struct SubPatchTest {
  size_t n; // number of coefs in total model
  size_t p; // polynomial order of basis
  std::vector<double> knots; // knot vector of total model
  size_t lknots0; // length of first subknot
  size_t lknots1; // length of second subknot
  std::vector<double> coefs; // control net of total model
  std::vector<double> coefs0; // control net of first subdivision
  std::vector<double> coefs1; // control net of second subdivision
  int dim; // dimensionality of embedding space
  bool rational; // rational or standard splines
  std::vector<int> mlge1; // global element number range owned by patch 1
  std::vector<int> mlgn1; // global node number range owned by patch 1
};


class TestGetSubPatch1D :
  public testing::Test,
  public testing::WithParamInterface<SubPatchTest>
{
};


TEST_P(TestGetSubPatch1D, SubPatch)
{
  Go::SplineCurve cur(GetParam().n, GetParam().p+1, GetParam().knots.begin(),
      GetParam().coefs.begin(), GetParam().dim, GetParam().rational);
  size_t numcoefs0 = cur.basis().numCoefs()/2 + GetParam().p;
  Go::SplineCurve cur0 = MultiPatchModelGenerator1D::getSubPatch(&cur,
      0, numcoefs0, GetParam().p+1);
  size_t numcoefs1 = cur.basis().numCoefs() - numcoefs0 + GetParam().p;
  size_t start1 = numcoefs0 - GetParam().p;
  Go::SplineCurve cur1 = MultiPatchModelGenerator1D::getSubPatch(&cur,
      start1, numcoefs1, GetParam().p+1);

  // Check first sub-knot-vector
  std::vector<double> arr0(GetParam().knots.begin(), GetParam().knots.begin()+GetParam().lknots0);
  std::vector<double> arr1(cur0.basis().begin(), cur0.basis().end());
  std::cout << "sub-knot-vector 0" << std::endl;
  check_vector_double_near(arr0, arr1);

  // Check second sub-knot-vector
  std::vector<double> arr2(GetParam().knots.end()-GetParam().lknots1, GetParam().knots.end());
  std::vector<double> arr3(cur1.basis().begin(), cur1.basis().end());
  std::cout << "sub-knot-vector 1" << std::endl;
  check_vector_double_near(arr2, arr3);

  // Check first sub-control net
  std::vector<double>::const_iterator i0 = GetParam().rational ? cur0.rcoefs_begin() : cur0.coefs_begin();
  std::vector<double>::const_iterator i1 = GetParam().rational ? cur0.rcoefs_end() : cur0.coefs_end();
  std::vector<double> arr4(i0, i1);
  std::cout << "sub-node-vector 0" << std::endl;
  check_vector_double_near(GetParam().coefs0, arr4);

  // Check second sub-control net
  std::vector<double>::const_iterator i2 = GetParam().rational ? cur1.rcoefs_begin() : cur1.coefs_begin();
  std::vector<double>::const_iterator i3 = GetParam().rational ? cur1.rcoefs_end() : cur1.coefs_end();
  std::vector<double> arr5(i2, i3);
  std::cout << "sub-node-vector 1" << std::endl;
  check_vector_double_near(GetParam().coefs1, arr5);

  // Evaluate function at points
  double xiA(0.3333), xiB(GetParam().p+1), xiC(4.2);
  std::string fstr = (GetParam().p == 2) ? "(1-x)*(3.14159-x)" : "x*(1-x)*(3.14159-x)";
  RealFunc* f = utl::parseRealFunc(fstr, "expression");
  Go::SplineCurve* fh = SplineUtils::project(&cur, *f);
  Go::SplineCurve fh0 = MultiPatchModelGenerator1D::getSubPatch(fh, 0, numcoefs0, GetParam().p+1);
  Go::SplineCurve fh1 = MultiPatchModelGenerator1D::getSubPatch(fh, start1, numcoefs1, GetParam().p+1);
  Go::Point fA, fB, fC, fA0, fB0, fB1, fC1;
  fh->point(fA, xiA);
  fh->point(fB, xiB);
  fh->point(fC, xiC);

  std::cout << "point evaluation" << std::endl;
  fh0.point(fA0, xiA);
  check_point_near(fA, fA0);
  fh0.point(fB0, xiB);
  check_point_near(fB, fB0);
  fh1.point(fB1, xiB);
  check_point_near(fB, fB1);
  fh1.point(fC1, xiC);
  check_point_near(fC, fC1);

  // Check FEM topology
  std::stringstream str0, str1;
  str0 << "100 1 0 0\n" << cur0;
  str1 << "100 1 0 0\n" << cur1;
  ASMs1D pch0, pch1;
  pch0.resetNumbering();
  pch0.read(str0);
  pch1.read(str1);
  std::cout << "element/node numbers" << std::endl;
  ASSERT_TRUE(pch0.generateFEMTopology());
  ASSERT_TRUE(pch1.generateFEMTopology());
  std::vector<int> myMLGE = pch1.getMyElementNums();
  check_vector_int_equals_range(myMLGE, GetParam().mlge1);
  std::vector<int> myMLGN = pch1.getMyNodeNums();
  check_vector_int_equals_range(myMLGN, GetParam().mlgn1);
}


const std::vector<SubPatchTest> SubPatch1D =
  {
   {7, 2, // non-rational, embedded in 1D, quadratic
    {0,0,0,1,2,3,4,5,5,5},
     8, 7,
    {0,1,2,3,4,5,6},
    {0,1,2,3,4},
    {3,4,5,6},
     1, false,
    {4,5},
    {6,9}},

   {8, 3, // non-rational, embedded in 1D, cubic
    {0,0,0,0,1,2,3,4,5,5,5,5},
     11, 8,
    {0,1,2,3,4,5,6,7},
    {0,1,2,3,4,5,6},
    {4,5,6,7},
     1, false,
    {5,5},
    {8,11}},

   {7, 2, // rational, embedded in 1D, quadratic
    {0,0,0,1,2,3,4,5,5,5},
     8, 7,
    {0,1,1,1,2,1,3,1,4,1,5,1,6,1},
    {0,1,1,1,2,1,3,1,4,1},
    {3,1,4,1,5,1,6,1},
     1, true,
    {4,5},
    {6,9}},

   {7, 2, // rational, embedded in 2D, quadratic
    {0,0,0,1,2,3,4,5,5,5},
     8, 7,
    {0,9,1,1,9,1,2,9,1,3,9,1,4,9,1,5,9,1,6,9,1},
    {0,9,1,1,9,1,2,9,1,3,9,1,4,9,1},
    {3,9,1,4,9,1,5,9,1,6,9,1},
     2, true,
    {4,5},
    {6,9}},
  };


INSTANTIATE_TEST_CASE_P(TestGetSubPatch1D,
                        TestGetSubPatch1D,
                        testing::ValuesIn(SubPatch1D));


class TestGetSubPatch2D :
  public testing::Test,
  public testing::WithParamInterface<SubPatchTest>
{
};


TEST_P(TestGetSubPatch2D, SubPatch)
{
  Go::SplineSurface srf(GetParam().n, GetParam().n, GetParam().p+1, GetParam().p+1,
      GetParam().knots.begin(), GetParam().knots.begin(),
      GetParam().coefs.begin(), GetParam().dim, GetParam().rational);
  size_t numcoefs0 = srf.basis_u().numCoefs()/2 + GetParam().p;
  Go::SplineSurface srf0 =
    MultiPatchModelGenerator2D::getSubPatch(&srf, {{0, 0}},
                                            {{numcoefs0, numcoefs0}},
                                            {{GetParam().p+1, GetParam().p+1}});
  size_t numcoefs1 = srf.basis_u().numCoefs() - numcoefs0 + GetParam().p;
  size_t start1 = numcoefs0 - GetParam().p;
  Go::SplineSurface srf1 =
    MultiPatchModelGenerator2D::getSubPatch(&srf, {{start1, 0}},
                                            {{numcoefs1, numcoefs0}},
                                            {{GetParam().p+1, GetParam().p+1}});

  // Check first sub-knot-vector
  std::vector<double> arr0(GetParam().knots.begin(), GetParam().knots.begin()+GetParam().lknots0);
  std::vector<double> arr1(srf0.basis_u().begin(), srf0.basis_u().end());
  std::cout << "sub-knot-vector 0" << std::endl;
  check_vector_double_near(arr0, arr1);

  // Check second sub-knot-vector
  std::vector<double> arr2(GetParam().knots.end()-GetParam().lknots1, GetParam().knots.end());
  std::vector<double> arr3(srf1.basis_u().begin(), srf1.basis_u().end());
  std::cout << "sub-knot-vector 1" << std::endl;
  check_vector_double_near(arr2, arr3);

  // Check first sub-control net
  std::vector<double>::const_iterator i0 = GetParam().rational ? srf0.rcoefs_begin() : srf0.coefs_begin();
  std::vector<double>::const_iterator i1 = GetParam().rational ? srf0.rcoefs_end() : srf0.coefs_end();
  std::vector<double> arr4(i0, i1);
  std::cout << "sub-node-vector 0" << std::endl;
  check_vector_double_near(GetParam().coefs0, arr4);

  // Check second sub-control net
  std::vector<double>::const_iterator i2 = GetParam().rational ? srf1.rcoefs_begin() : srf1.coefs_begin();
  std::vector<double>::const_iterator i3 = GetParam().rational ? srf1.rcoefs_end() : srf1.coefs_end();
  std::vector<double> arr5(i2, i3);
  std::cout << "sub-node-vector 1" << std::endl;
  check_vector_double_near(GetParam().coefs1, arr5);

  // Evaluate function at points
  double xiA(0.3333), xiB(2), xiC(4.2);
  RealFunc* f = utl::parseRealFunc("(1-y*y)*(3.14159-x)", "expression");
  Go::SplineSurface* fh = SplineUtils::project(&srf, *f);
  Go::SplineSurface fh0 =
    MultiPatchModelGenerator2D::getSubPatch(fh, {{0, 0}},
                                            {{numcoefs0, numcoefs0}},
                                            {{GetParam().p+1, GetParam().p+1}});
  Go::SplineSurface fh1 =
    MultiPatchModelGenerator2D::getSubPatch(fh, {{start1, 0}},
                                            {{numcoefs1, numcoefs0}},
                                            {{GetParam().p+1, GetParam().p+1}});
  Go::Point fA, fB, fC, fA0, fB0, fB1, fC1;
  fh->point(fA, xiA, xiA);
  fh->point(fB, xiB, xiA);
  fh->point(fC, xiC, xiA);

  std::cout << "point evaluation" << std::endl;
  fh0.point(fA0, xiA, xiA);
  check_point_near(fA, fA0);
  fh0.point(fB0, xiB, xiA);
  check_point_near(fB, fB0);
  fh1.point(fB1, xiB, xiA);
  check_point_near(fB, fB1);
  fh1.point(fC1, xiC, xiA);
  check_point_near(fC, fC1);

  // Check FEM topology
  std::stringstream str0, str1;
  str0 << "200 1 0 0\n" << srf0;
  str1 << "200 1 0 0\n" << srf1;
  ASMs2D pch0, pch1;
  pch0.resetNumbering();
  pch0.read(str0);
  pch1.read(str1);
  std::cout << "element/node numbers" << std::endl;
  ASSERT_TRUE(pch0.generateFEMTopology());
  ASSERT_TRUE(pch1.generateFEMTopology());
  std::vector<int> myMLGE = pch1.getMyElementNums();
  check_vector_int_equals_range(myMLGE, GetParam().mlge1);
  std::vector<int> myMLGN = pch1.getMyNodeNums();
  check_vector_int_equals_range(myMLGN, GetParam().mlgn1);
}


const std::vector<SubPatchTest> SubPatch2D =
  {
   {7, 2, // non-rational, embedded in 2D
    {0,0,0,1,2,3,4,5,5,5},
     8, 7,
    {0,0,1,0,2,0,3,0,4,0,5,0,6,0,
     0,1,1,1,2,1,3,1,4,1,5,1,6,1,
     0,2,1,2,2,2,3,2,4,2,5,2,6,2,
     0,3,1,3,2,3,3,3,4,3,5,3,6,3,
     0,4,1,4,2,4,3,4,4,4,5,4,6,4,
     0,5,1,5,2,5,3,5,4,5,5,5,6,5,
     0,6,1,6,2,6,3,6,4,6,5,6,6,6},
    {0,0,1,0,2,0,3,0,4,0,
     0,1,1,1,2,1,3,1,4,1,
     0,2,1,2,2,2,3,2,4,2,
     0,3,1,3,2,3,3,3,4,3,
     0,4,1,4,2,4,3,4,4,4},
    {3,0,4,0,5,0,6,0,
     3,1,4,1,5,1,6,1,
     3,2,4,2,5,2,6,2,
     3,3,4,3,5,3,6,3,
     3,4,4,4,5,4,6,4},
     2, false,
    {10,15},
    {26,45}},

   {7, 2, // rational, embedded in 2D
    {0,0,0,1,2,3,4,5,5,5},
     8, 7,
    {0,0,1,1,0,1,2,0,1,3,0,1,4,0,1,5,0,1,6,0,1,
     0,1,1,1,1,1,2,1,1,3,1,1,4,1,1,5,1,1,6,1,1,
     0,2,1,1,2,1,2,2,1,3,2,1,4,2,1,5,2,1,6,2,1,
     0,3,1,1,3,1,2,3,1,3,3,1,4,3,1,5,3,1,6,3,1,
     0,4,1,1,4,1,2,4,1,3,4,1,4,4,1,5,4,1,6,4,1,
     0,5,1,1,5,1,2,5,1,3,5,1,4,5,1,5,5,1,6,5,1,
     0,6,1,1,6,1,2,6,1,3,6,1,4,6,1,5,6,1,6,6,1},
    {0,0,1,1,0,1,2,0,1,3,0,1,4,0,1,
     0,1,1,1,1,1,2,1,1,3,1,1,4,1,1,
     0,2,1,1,2,1,2,2,1,3,2,1,4,2,1,
     0,3,1,1,3,1,2,3,1,3,3,1,4,3,1,
     0,4,1,1,4,1,2,4,1,3,4,1,4,4,1},
    {3,0,1,4,0,1,5,0,1,6,0,1,
     3,1,1,4,1,1,5,1,1,6,1,1,
     3,2,1,4,2,1,5,2,1,6,2,1,
     3,3,1,4,3,1,5,3,1,6,3,1,
     3,4,1,4,4,1,5,4,1,6,4,1},
     2, true,
    {10,15},
    {26,45}},

   {7, 2, // rational, embedded in 3D
    {0,0,0,1,2,3,4,5,5,5},
     8, 7,
    {0,0,9,1,1,0,9,1,2,0,9,1,3,0,9,1,4,0,9,1,5,0,9,1,6,0,9,1,
     0,1,9,1,1,1,9,1,2,1,9,1,3,1,9,1,4,1,9,1,5,1,9,1,6,1,9,1,
     0,2,9,1,1,2,9,1,2,2,9,1,3,2,9,1,4,2,9,1,5,2,9,1,6,2,9,1,
     0,3,9,1,1,3,9,1,2,3,9,1,3,3,9,1,4,3,9,1,5,3,9,1,6,3,9,1,
     0,4,9,1,1,4,9,1,2,4,9,1,3,4,9,1,4,4,9,1,5,4,9,1,6,4,9,1,
     0,5,9,1,1,5,9,1,2,5,9,1,3,5,9,1,4,5,9,1,5,5,9,1,6,5,9,1,
     0,6,9,1,1,6,9,1,2,6,9,1,3,6,9,1,4,6,9,1,5,6,9,1,6,6,9,1},
    {0,0,9,1,1,0,9,1,2,0,9,1,3,0,9,1,4,0,9,1,
     0,1,9,1,1,1,9,1,2,1,9,1,3,1,9,1,4,1,9,1,
     0,2,9,1,1,2,9,1,2,2,9,1,3,2,9,1,4,2,9,1,
     0,3,9,1,1,3,9,1,2,3,9,1,3,3,9,1,4,3,9,1,
     0,4,9,1,1,4,9,1,2,4,9,1,3,4,9,1,4,4,9,1},
    {3,0,9,1,4,0,9,1,5,0,9,1,6,0,9,1,
     3,1,9,1,4,1,9,1,5,1,9,1,6,1,9,1,
     3,2,9,1,4,2,9,1,5,2,9,1,6,2,9,1,
     3,3,9,1,4,3,9,1,5,3,9,1,6,3,9,1,
     3,4,9,1,4,4,9,1,5,4,9,1,6,4,9,1},
     3, true,
    {10,15},
    {26,45}},
  };


INSTANTIATE_TEST_CASE_P(TestGetSubPatch2D,
                        TestGetSubPatch2D,
                        testing::ValuesIn(SubPatch2D));


class TestGetSubPatch3D :
  public testing::Test,
  public testing::WithParamInterface<SubPatchTest>
{
};


TEST_P(TestGetSubPatch3D, SubPatch)
{
  Go::SplineVolume vol(GetParam().n, GetParam().n, GetParam().n,
                       GetParam().p+1, GetParam().p+1, GetParam().p+1,
                       GetParam().knots.begin(), GetParam().knots.begin(),
                       GetParam().knots.begin(),
                       GetParam().coefs.begin(), GetParam().dim, GetParam().rational);
  size_t numcoefs0 = vol.basis(0).numCoefs()/2 + GetParam().p;
  Go::SplineVolume vol0 =
    MultiPatchModelGenerator3D::getSubPatch(&vol, {{0, 0, 0}},
                                            {{numcoefs0, numcoefs0, numcoefs0}},
                                            {{GetParam().p+1, GetParam().p+1,
                                              GetParam().p+1}});
  size_t numcoefs1 = vol.basis(0).numCoefs() - numcoefs0 + GetParam().p;
  size_t start1 = numcoefs0 - GetParam().p;
  Go::SplineVolume vol1 =
    MultiPatchModelGenerator3D::getSubPatch(&vol, {{start1, 0, 0}},
                                            {{numcoefs1, numcoefs0, numcoefs0}},
                                            {{GetParam().p+1, GetParam().p+1,
                                              GetParam().p+1}});

  // Check first sub-knot-vector
  std::vector<double> arr0(GetParam().knots.begin(), GetParam().knots.begin()+GetParam().lknots0);
  std::vector<double> arr1(vol0.basis(0).begin(), vol0.basis(0).end());
  std::cout << "sub-knot-vector 0" << std::endl;
  check_vector_double_near(arr0, arr1);

  // Check second sub-knot-vector
  std::vector<double> arr2(GetParam().knots.end()-GetParam().lknots1, GetParam().knots.end());
  std::vector<double> arr3(vol1.basis(0).begin(), vol1.basis(0).end());
  std::cout << "sub-knot-vector 1" << std::endl;
  check_vector_double_near(arr2, arr3);
  std::vector<double> arrX(vol.coefs_begin(), vol.coefs_end());

  // Check first sub-control net
  std::vector<double>::const_iterator i0 = GetParam().rational ? vol0.rcoefs_begin() : vol0.coefs_begin();
  std::vector<double>::const_iterator i1 = GetParam().rational ? vol0.rcoefs_end() : vol0.coefs_end();
  std::vector<double> arr4(i0, i1);
  std::cout << "sub-node-vector 0" << std::endl;
  check_vector_double_near(GetParam().coefs0, arr4);

  // Check second sub-control net
  std::vector<double>::const_iterator i2 = GetParam().rational ? vol1.rcoefs_begin() : vol1.coefs_begin();
  std::vector<double>::const_iterator i3 = GetParam().rational ? vol1.rcoefs_end() : vol1.coefs_end();
  std::vector<double> arr5(i2, i3);
  std::cout << "sub-node-vector 1" << std::endl;
  check_vector_double_near(GetParam().coefs1, arr5);

  // Evaluate function at points
  double xiA(0.3333), xiB(2), xiC(4.2);
  RealFunc* f = utl::parseRealFunc("(1-y*y)*(3.14159-x)*(1-z)", "expression");
  Go::SplineVolume* fh = SplineUtils::project(&vol, *f);
  Go::SplineVolume fh0 =
    MultiPatchModelGenerator3D::getSubPatch(fh, {{0, 0, 0}},
                                            {{numcoefs0, numcoefs0, numcoefs0}},
                                            {{GetParam().p+1, GetParam().p+1,
                                              GetParam().p+1}});
  Go::SplineVolume fh1 =
    MultiPatchModelGenerator3D::getSubPatch(fh, {{start1, 0, 0}},
                                            {{numcoefs1, numcoefs0, numcoefs0}},
                                            {{GetParam().p+1, GetParam().p+1,
                                              GetParam().p+1}});
  Go::Point fA, fB, fC, fA0, fB0, fB1, fC1;
  fh->point(fA, xiA, xiA, xiA);
  fh->point(fB, xiB, xiA, xiA);
  fh->point(fC, xiC, xiA, xiA);

  std::cout << "point evaluation" << std::endl;
  fh0.point(fA0, xiA, xiA, xiA);
  check_point_near(fA, fA0);
  fh0.point(fB0, xiB, xiA, xiA);
  check_point_near(fB, fB0);
  fh1.point(fB1, xiB, xiA, xiA);
  check_point_near(fB, fB1);
  fh1.point(fC1, xiC, xiA, xiA);
  check_point_near(fC, fC1);

  // Check FEM topology
  std::stringstream str0, str1;
  str0 << "700 1 0 0\n" << vol0;
  str1 << "700 1 0 0\n" << vol1;
  ASMs3D pch0, pch1;
  pch0.resetNumbering();
  pch0.read(str0);
  pch1.read(str1);
  std::cout << "element/node numbers" << std::endl;
  ASSERT_TRUE(pch0.generateFEMTopology());
  ASSERT_TRUE(pch1.generateFEMTopology());
  std::vector<int> myMLGE = pch1.getMyElementNums();
  check_vector_int_equals_range(myMLGE, GetParam().mlge1);
  std::vector<int> myMLGN = pch1.getMyNodeNums();
  check_vector_int_equals_range(myMLGN, GetParam().mlgn1);
}


const std::vector<SubPatchTest> SubPatch3D =
  {
   {7, 2, // non-rational, embedded in 3D
    {0,0,0,1,2,3,4,5,5,5},
     8, 7,
    {0,0,0,1,0,0,2,0,0,3,0,0,4,0,0,5,0,0,6,0,0,
     0,1,0,1,1,0,2,1,0,3,1,0,4,1,0,5,1,0,6,1,0,
     0,2,0,1,2,0,2,2,0,3,2,0,4,2,0,5,2,0,6,2,0,
     0,3,0,1,3,0,2,3,0,3,3,0,4,3,0,5,3,0,6,3,0,
     0,4,0,1,4,0,2,4,0,3,4,0,4,4,0,5,4,0,6,4,0,
     0,5,0,1,5,0,2,5,0,3,5,0,4,5,0,5,5,0,6,5,0,
     0,6,0,1,6,0,2,6,0,3,6,0,4,6,0,5,6,0,6,6,0,
     0,0,1,1,0,1,2,0,1,3,0,1,4,0,1,5,0,1,6,0,1,
     0,1,1,1,1,1,2,1,1,3,1,1,4,1,1,5,1,1,6,1,1,
     0,2,1,1,2,1,2,2,1,3,2,1,4,2,1,5,2,1,6,2,1,
     0,3,1,1,3,1,2,3,1,3,3,1,4,3,1,5,3,1,6,3,1,
     0,4,1,1,4,1,2,4,1,3,4,1,4,4,1,5,4,1,6,4,1,
     0,5,1,1,5,1,2,5,1,3,5,1,4,5,1,5,5,1,6,5,1,
     0,6,1,1,6,1,2,6,1,3,6,1,4,6,1,5,6,1,6,6,1,
     0,0,2,1,0,2,2,0,2,3,0,2,4,0,2,5,0,2,6,0,2,
     0,1,2,1,1,2,2,1,2,3,1,2,4,1,2,5,1,2,6,1,2,
     0,2,2,1,2,2,2,2,2,3,2,2,4,2,2,5,2,2,6,2,2,
     0,3,2,1,3,2,2,3,2,3,3,2,4,3,2,5,3,2,6,3,2,
     0,4,2,1,4,2,2,4,2,3,4,2,4,4,2,5,4,2,6,4,2,
     0,5,2,1,5,2,2,5,2,3,5,2,4,5,2,5,5,2,6,5,2,
     0,6,2,1,6,2,2,6,2,3,6,2,4,6,2,5,6,2,6,6,2,
     0,0,3,1,0,3,2,0,3,3,0,3,4,0,3,5,0,3,6,0,3,
     0,1,3,1,1,3,2,1,3,3,1,3,4,1,3,5,1,3,6,1,3,
     0,2,3,1,2,3,2,2,3,3,2,3,4,2,3,5,2,3,6,2,3,
     0,3,3,1,3,3,2,3,3,3,3,3,4,3,3,5,3,3,6,3,3,
     0,4,3,1,4,3,2,4,3,3,4,3,4,4,3,5,4,3,6,4,3,
     0,5,3,1,5,3,2,5,3,3,5,3,4,5,3,5,5,3,6,5,3,
     0,6,3,1,6,3,2,6,3,3,6,3,4,6,3,5,6,3,6,6,3,
     0,0,4,1,0,4,2,0,4,3,0,4,4,0,4,5,0,4,6,0,4,
     0,1,4,1,1,4,2,1,4,3,1,4,4,1,4,5,1,4,6,1,4,
     0,2,4,1,2,4,2,2,4,3,2,4,4,2,4,5,2,4,6,2,4,
     0,3,4,1,3,4,2,3,4,3,3,4,4,3,4,5,3,4,6,3,4,
     0,4,4,1,4,4,2,4,4,3,4,4,4,4,4,5,4,4,6,4,4,
     0,5,4,1,5,4,2,5,4,3,5,4,4,5,4,5,5,4,6,5,4,
     0,6,4,1,6,4,2,6,4,3,6,4,4,6,4,5,6,4,6,6,4,
     0,0,5,1,0,5,2,0,5,3,0,5,4,0,5,5,0,5,6,0,5,
     0,1,5,1,1,5,2,1,5,3,1,5,4,1,5,5,1,5,6,1,5,
     0,2,5,1,2,5,2,2,5,3,2,5,4,2,5,5,2,5,6,2,5,
     0,3,5,1,3,5,2,3,5,3,3,5,4,3,5,5,3,5,6,3,5,
     0,4,5,1,4,5,2,4,5,3,4,5,4,4,5,5,4,5,6,4,5,
     0,5,5,1,5,5,2,5,5,3,5,5,4,5,5,5,5,5,6,5,5,
     0,6,5,1,6,5,2,6,5,3,6,5,4,6,5,5,6,5,6,6,5,
     0,0,6,1,0,6,2,0,6,3,0,6,4,0,6,5,0,6,6,0,6,
     0,1,6,1,1,6,2,1,6,3,1,6,4,1,6,5,1,6,6,1,6,
     0,2,6,1,2,6,2,2,6,3,2,6,4,2,6,5,2,6,6,2,6,
     0,3,6,1,3,6,2,3,6,3,3,6,4,3,6,5,3,6,6,3,6,
     0,4,6,1,4,6,2,4,6,3,4,6,4,4,6,5,4,6,6,4,6,
     0,5,6,1,5,6,2,5,6,3,5,6,4,5,6,5,5,6,6,5,6,
     0,6,6,1,6,6,2,6,6,3,6,6,4,6,6,5,6,6,6,6,6},
    {0,0,0,1,0,0,2,0,0,3,0,0,4,0,0,
     0,1,0,1,1,0,2,1,0,3,1,0,4,1,0,
     0,2,0,1,2,0,2,2,0,3,2,0,4,2,0,
     0,3,0,1,3,0,2,3,0,3,3,0,4,3,0,
     0,4,0,1,4,0,2,4,0,3,4,0,4,4,0,
     0,0,1,1,0,1,2,0,1,3,0,1,4,0,1,
     0,1,1,1,1,1,2,1,1,3,1,1,4,1,1,
     0,2,1,1,2,1,2,2,1,3,2,1,4,2,1,
     0,3,1,1,3,1,2,3,1,3,3,1,4,3,1,
     0,4,1,1,4,1,2,4,1,3,4,1,4,4,1,
     0,0,2,1,0,2,2,0,2,3,0,2,4,0,2,
     0,1,2,1,1,2,2,1,2,3,1,2,4,1,2,
     0,2,2,1,2,2,2,2,2,3,2,2,4,2,2,
     0,3,2,1,3,2,2,3,2,3,3,2,4,3,2,
     0,4,2,1,4,2,2,4,2,3,4,2,4,4,2,
     0,0,3,1,0,3,2,0,3,3,0,3,4,0,3,
     0,1,3,1,1,3,2,1,3,3,1,3,4,1,3,
     0,2,3,1,2,3,2,2,3,3,2,3,4,2,3,
     0,3,3,1,3,3,2,3,3,3,3,3,4,3,3,
     0,4,3,1,4,3,2,4,3,3,4,3,4,4,3,
     0,0,4,1,0,4,2,0,4,3,0,4,4,0,4,
     0,1,4,1,1,4,2,1,4,3,1,4,4,1,4,
     0,2,4,1,2,4,2,2,4,3,2,4,4,2,4,
     0,3,4,1,3,4,2,3,4,3,3,4,4,3,4,
     0,4,4,1,4,4,2,4,4,3,4,4,4,4,4},
    {3,0,0,4,0,0,5,0,0,6,0,0,
     3,1,0,4,1,0,5,1,0,6,1,0,
     3,2,0,4,2,0,5,2,0,6,2,0,
     3,3,0,4,3,0,5,3,0,6,3,0,
     3,4,0,4,4,0,5,4,0,6,4,0,
     3,0,1,4,0,1,5,0,1,6,0,1,
     3,1,1,4,1,1,5,1,1,6,1,1,
     3,2,1,4,2,1,5,2,1,6,2,1,
     3,3,1,4,3,1,5,3,1,6,3,1,
     3,4,1,4,4,1,5,4,1,6,4,1,
     3,0,2,4,0,2,5,0,2,6,0,2,
     3,1,2,4,1,2,5,1,2,6,1,2,
     3,2,2,4,2,2,5,2,2,6,2,2,
     3,3,2,4,3,2,5,3,2,6,3,2,
     3,4,2,4,4,2,5,4,2,6,4,2,
     3,0,3,4,0,3,5,0,3,6,0,3,
     3,1,3,4,1,3,5,1,3,6,1,3,
     3,2,3,4,2,3,5,2,3,6,2,3,
     3,3,3,4,3,3,5,3,3,6,3,3,
     3,4,3,4,4,3,5,4,3,6,4,3,
     3,0,4,4,0,4,5,0,4,6,0,4,
     3,1,4,4,1,4,5,1,4,6,1,4,
     3,2,4,4,2,4,5,2,4,6,2,4,
     3,3,4,4,3,4,5,3,4,6,3,4,
     3,4,4,4,4,4,5,4,4,6,4,4},
     3, false,
    {28,45},
    {126,225}},
  };


INSTANTIATE_TEST_CASE_P(TestGetSubPatch3D,
                        TestGetSubPatch3D,
                        testing::ValuesIn(SubPatch3D));

TEST(TestMultiPatchModelGenerator3D, InnerPatches)
{
  TiXmlDocument doc;
  doc.Parse("<geometry nx=\"3\" ny=\"3\" nz=\"3\"  sets=\"true\"/>");
  TestModelGeneratorWrapper<MultiPatchModelGenerator3D> gen(doc.RootElement());
  SIM3D sim;
  gen.createTopologySets(sim);
  ASSERT_EQ(sim.topology("InnerPatches").size(), 1u);
  ASSERT_EQ(sim.topology("InnerPatches").begin()->patch, 14u);
}


const std::vector<GeomTest> geometry2D =
  {{"<geometry sets=\"true\"/>", 2,
    "200 1 0 0\n"
    "2 0\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "0 0\n"
    "1 0\n"
    "0 1\n"
    "1 1\n",
    "Boundary: 1 1 1 1 2 1 1 3 1 1 4 1 \n"
    "Corners: 1 1 0 1 2 0 1 3 0 1 4 0 \n"
    "Edge1: 1 1 1 \n"
    "Edge1Patches: 1 0 2 \n"
    "Edge2: 1 2 1 \n"
    "Edge2Patches: 1 0 2 \n"
    "Edge3: 1 3 1 \n"
    "Edge3Patches: 1 0 2 \n"
    "Edge4: 1 4 1 \n"
    "Edge4Patches: 1 0 2 \n"
    "InnerPatches: \n"
    "Vertex1: 1 1 0 \n"
    "Vertex2: 1 2 0 \n"
    "Vertex3: 1 3 0 \n"
    "Vertex4: 1 4 0 \n"},

   {"<geometry rational=\"1\"/>", 2,
    "200 1 0 0\n"
    "2 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "0 0 1.0\n"
    "1 0 1.0\n"
    "0 1 1.0\n"
    "1 1 1.0\n", ""},

   {"<geometry scale=\"2\"/>", 2,
     "200 1 0 0\n"
     "2 0\n"
     "2 2\n"
     "0 0 1 1\n"
     "2 2\n"
     "0 0 1 1\n"
     "0 0\n"
     "2 0\n"
     "0 2\n"
     "2 2\n", ""},

   {"<geometry X0=\"2 0\"/>", 2,
    "200 1 0 0\n"
     "2 0\n"
     "2 2\n"
     "0 0 1 1\n"
     "2 2\n"
     "0 0 1 1\n"
     "2 0\n"
     "3 0\n"
     "2 1\n"
     "3 1\n"},

    {"<geometry X0=\"0 2\"/>", 2,
     "200 1 0 0\n"
     "2 0\n"
     "2 2\n"
     "0 0 1 1\n"
     "2 2\n"
     "0 0 1 1\n"
     "0 2\n"
     "1 2\n"
     "0 3\n"
     "1 3\n", ""},

    {"<geometry Lx=\"2\"/>", 2,
     "200 1 0 0\n"
     "2 0\n"
     "2 2\n"
     "0 0 1 1\n"
     "2 2\n"
     "0 0 1 1\n"
     "0 0\n"
     "2 0\n"
     "0 1\n"
     "2 1\n", ""},

    {"<geometry Ly=\"2\"/>", 2,
     "200 1 0 0\n"
     "2 0\n"
     "2 2\n"
     "0 0 1 1\n"
     "2 2\n"
     "0 0 1 1\n"
     "0 0\n"
     "1 0\n"
     "0 2\n"
     "1 2\n", ""},

    {"<geometry sets=\"true\" nx=\"2\"/>", 2,
     "200 1 0 0\n"
     "2 0\n"
     "2 2\n"
     "0 0 1 1\n"
     "2 2\n"
     "0 0 1 1\n"
     "0 0\n"
     "0.5 0\n"
     "0 1\n"
     "0.5 1\n"
     "200 1 0 0\n"
     "2 0\n"
     "2 2\n"
     "0 0 1 1\n"
     "2 2\n"
     "0 0 1 1\n"
     "0.5 0\n"
     "1 0\n"
     "0.5 1\n"
     "1 1\n",
     "Boundary: 1 1 1 1 3 1 1 4 1 2 2 1 2 3 1 2 4 1 \n"
     "Corners: 1 1 0 1 3 0 2 2 0 2 4 0 \n"
     "Edge1: 1 1 1 \n"
     "Edge1Patches: 1 0 2 \n"
     "Edge2: 2 2 1 \n"
     "Edge2Patches: 2 0 2 \n"
     "Edge3: 1 3 1 2 3 1 \n"
     "Edge3Patches: 1 0 2 2 0 2 \n"
     "Edge4: 1 4 1 2 4 1 \n"
     "Edge4Patches: 1 0 2 2 0 2 \n"
     "InnerPatches: \n"
     "Vertex1: 1 1 0 \n"
     "Vertex2: 2 2 0 \n"
     "Vertex3: 1 3 0 \n"
     "Vertex4: 2 4 0 \n"},

    {"<geometry sets=\"true\" ny=\"2\"/>", 2,
     "200 1 0 0\n"
     "2 0\n"
     "2 2\n"
     "0 0 1 1\n"
     "2 2\n"
     "0 0 1 1\n"
     "0 0\n"
     "1 0\n"
     "0 0.5\n"
     "1 0.5\n"
     "200 1 0 0\n"
     "2 0\n"
     "2 2\n"
     "0 0 1 1\n"
     "2 2\n"
     "0 0 1 1\n"
     "0 0.5\n"
     "1 0.5\n"
     "0 1\n"
     "1 1\n",
     "Boundary: 1 1 1 1 2 1 1 3 1 2 1 1 2 2 1 2 4 1 \n"
     "Corners: 1 1 0 1 2 0 2 3 0 2 4 0 \n"
     "Edge1: 1 1 1 2 1 1 \n"
     "Edge1Patches: 1 0 2 2 0 2 \n"
     "Edge2: 1 2 1 2 2 1 \n"
     "Edge2Patches: 1 0 2 2 0 2 \n"
     "Edge3: 1 3 1 \n"
     "Edge3Patches: 1 0 2 \n"
     "Edge4: 2 4 1 \n"
     "Edge4Patches: 2 0 2 \n"
     "InnerPatches: \n"
     "Vertex1: 1 1 0 \n"
     "Vertex2: 1 2 0 \n"
     "Vertex3: 2 3 0 \n"
     "Vertex4: 2 4 0 \n"},

    {"<geometry sets=\"true\" nx=\"2\" ny=\"2\"/>", 2,
     "200 1 0 0\n"
     "2 0\n"
     "2 2\n"
     "0 0 1 1\n"
     "2 2\n"
     "0 0 1 1\n"
     "0 0\n"
     "0.5 0\n"
     "0 0.5\n"
     "0.5 0.5\n"
     "200 1 0 0\n"
     "2 0\n"
     "2 2\n"
     "0 0 1 1\n"
     "2 2\n"
     "0 0 1 1\n"
     "0.5 0\n"
     "1 0\n"
     "0.5 0.5\n"
     "1 0.5\n"
     "200 1 0 0\n"
     "2 0\n"
     "2 2\n"
     "0 0 1 1\n"
     "2 2\n"
     "0 0 1 1\n"
     "0 0.5\n"
     "0.5 0.5\n"
     "0 1\n"
     "0.5 1\n"
     "200 1 0 0\n"
     "2 0\n"
     "2 2\n"
     "0 0 1 1\n"
     "2 2\n"
     "0 0 1 1\n"
     "0.5 0.5\n"
     "1 0.5\n"
     "0.5 1\n"
     "1 1\n",
     "Boundary: 1 1 1 1 3 1 2 2 1 2 3 1 3 1 1 3 4 1 4 2 1 4 4 1 \n"
     "Corners: 1 1 0 2 2 0 3 3 0 4 4 0 \n"
     "Edge1: 1 1 1 3 1 1 \n"
     "Edge1Patches: 1 0 2 3 0 2 \n"
     "Edge2: 2 2 1 4 2 1 \n"
     "Edge2Patches: 2 0 2 4 0 2 \n"
     "Edge3: 1 3 1 2 3 1 \n"
     "Edge3Patches: 1 0 2 2 0 2 \n"
     "Edge4: 3 4 1 4 4 1 \n"
     "Edge4Patches: 3 0 2 4 0 2 \n"
     "InnerPatches: \n"
     "Vertex1: 1 1 0 \n"
     "Vertex2: 2 2 0 \n"
     "Vertex3: 3 3 0 \n"
     "Vertex4: 4 4 0 \n"}};


INSTANTIATE_TEST_CASE_P(TestMultiPatchModelGenerator2D,
                        TestMultiPatchModelGenerator2D,
                        testing::ValuesIn(geometry2D));


const std::vector<GeomTest> geometry3D =
  {{"<geometry sets=\"true\"/>", 3,
    "700 1 0 0\n"
    "3 0\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "0 0 0\n"
    "1 0 0\n"
    "0 1 0\n"
    "1 1 0\n"
    "0 0 1\n"
    "1 0 1\n"
    "0 1 1\n"
    "1 1 1\n",
    "Boundary: 1 1 2 1 2 2 1 3 2 1 4 2 1 5 2 1 6 2 \n"
    "Corners: 1 1 0 1 2 0 1 3 0 1 4 0 1 5 0 1 6 0 1 7 0 1 8 0 \n"
    "Edge1: 1 1 1 \n"
    "Edge10: 1 10 1 \n"
    "Edge11: 1 11 1 \n"
    "Edge12: 1 12 1 \n"
    "Edge2: 1 2 1 \n"
    "Edge3: 1 3 1 \n"
    "Edge4: 1 4 1 \n"
    "Edge5: 1 5 1 \n"
    "Edge6: 1 6 1 \n"
    "Edge7: 1 7 1 \n"
    "Edge8: 1 8 1 \n"
    "Edge9: 1 9 1 \n"
    "Face1: 1 1 2 \n"
    "Face1Patches: 1 0 3 \n"
    "Face2: 1 2 2 \n"
    "Face2Patches: 1 0 3 \n"
    "Face3: 1 3 2 \n"
    "Face3Patches: 1 0 3 \n"
    "Face4: 1 4 2 \n"
    "Face4Patches: 1 0 3 \n"
    "Face5: 1 5 2 \n"
    "Face5Patches: 1 0 3 \n"
    "Face6: 1 6 2 \n"
    "Face6Patches: 1 0 3 \n"
    "Frame: 1 1 1 1 2 1 1 3 1 1 4 1 1 5 1 1 6 1 1 7 1 1 8 1 1 9 1 1 10 1 1 11 1 1 12 1 \n"
    "InnerPatches: \n"
    "Vertex1: 1 1 0 \n"
    "Vertex2: 1 2 0 \n"
    "Vertex3: 1 3 0 \n"
    "Vertex4: 1 4 0 \n"
    "Vertex5: 1 5 0 \n"
    "Vertex6: 1 6 0 \n"
    "Vertex7: 1 7 0 \n"
    "Vertex8: 1 8 0 \n"},

  {"<geometry rational=\"1\"/>", 3,
   "700 1 0 0\n"
   "3 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "0 0 0 1.0\n"
   "1 0 0 1.0\n"
   "0 1 0 1.0\n"
   "1 1 0 1.0\n"
   "0 0 1 1.0\n"
   "1 0 1 1.0\n"
   "0 1 1 1.0\n"
   "1 1 1 1.0\n", ""},

  {"<geometry scale=\"2\"/>", 3,
   "700 1 0 0\n"
   "3 0\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "0 0 0\n"
   "2 0 0\n"
   "0 2 0\n"
   "2 2 0\n"
   "0 0 2\n"
   "2 0 2\n"
   "0 2 2\n"
   "2 2 2\n", ""},

  {"<geometry X0=\"2 0 0\"/>", 3,
   "700 1 0 0\n"
   "3 0\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 0 0\n"
   "3 0 0\n"
   "2 1 0\n"
   "3 1 0\n"
   "2 0 1\n"
   "3 0 1\n"
   "2 1 1\n"
   "3 1 1\n", ""},

  {"<geometry X0=\"0 2 0\"/>", 3,
   "700 1 0 0\n"
   "3 0\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "0 2 0\n"
   "1 2 0\n"
   "0 3 0\n"
   "1 3 0\n"
   "0 2 1\n"
   "1 2 1\n"
   "0 3 1\n"
   "1 3 1\n", ""},

  {"<geometry X0=\"0 0 2\"/>", 3,
   "700 1 0 0\n"
   "3 0\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "0 0 2\n"
   "1 0 2\n"
   "0 1 2\n"
   "1 1 2\n"
   "0 0 3\n"
   "1 0 3\n"
   "0 1 3\n"
   "1 1 3\n", ""},

  {"<geometry Lx=\"2\"/>", 3,
   "700 1 0 0\n"
   "3 0\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "0 0 0\n"
   "2 0 0\n"
   "0 1 0\n"
   "2 1 0\n"
   "0 0 1\n"
   "2 0 1\n"
   "0 1 1\n"
   "2 1 1\n", ""},

  {"<geometry Ly=\"2\"/>", 3,
   "700 1 0 0\n"
   "3 0\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "0 0 0\n"
   "1 0 0\n"
   "0 2 0\n"
   "1 2 0\n"
   "0 0 1\n"
   "1 0 1\n"
   "0 2 1\n"
   "1 2 1\n", ""},

  {"<geometry Lz=\"2\"/>", 3,
   "700 1 0 0\n"
     "3 0\n"
     "2 2\n"
     "0 0 1 1\n"
     "2 2\n"
     "0 0 1 1\n"
     "2 2\n"
     "0 0 1 1\n"
     "0 0 0\n"
     "1 0 0\n"
     "0 1 0\n"
     "1 1 0\n"
     "0 0 2\n"
     "1 0 2\n"
     "0 1 2\n"
     "1 1 2\n", ""},

  {"<geometry sets=\"true\" nx=\"2\"/>", 3,
   "700 1 0 0\n"
   "3 0\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "0 0 0\n"
   "0.5 0 0\n"
   "0 1 0\n"
   "0.5 1 0\n"
   "0 0 1\n"
   "0.5 0 1\n"
   "0 1 1\n"
   "0.5 1 1\n"
   "700 1 0 0\n"
   "3 0\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "0.5 0 0\n"
   "1 0 0\n"
   "0.5 1 0\n"
   "1 1 0\n"
   "0.5 0 1\n"
   "1 0 1\n"
   "0.5 1 1\n"
   "1 1 1\n",
   "Boundary: 1 1 2 1 3 2 1 4 2 1 5 2 1 6 2 2 2 2 2 3 2 2 4 2 2 5 2 2 6 2 \n"
   "Corners: 1 1 0 1 3 0 1 5 0 1 7 0 2 2 0 2 4 0 2 6 0 2 8 0 \n"
   "Edge1: 1 1 1 \n"
   "Edge10: 1 10 1 2 10 1 \n"
   "Edge11: 1 11 1 2 11 1 \n"
   "Edge12: 1 12 1 2 12 1 \n"
   "Edge2: 2 2 1 \n"
   "Edge3: 1 3 1 \n"
   "Edge4: 2 4 1 \n"
   "Edge5: 1 5 1 \n"
   "Edge6: 2 6 1 \n"
   "Edge7: 1 7 1 \n"
   "Edge8: 2 8 1 \n"
   "Edge9: 1 9 1 2 9 1 \n"
   "Face1: 1 1 2 \n"
   "Face1Patches: 1 0 3 \n"
   "Face2: 2 2 2 \n"
   "Face2Patches: 2 0 3 \n"
   "Face3: 1 3 2 2 3 2 \n"
   "Face3Patches: 1 0 3 2 0 3 \n"
   "Face4: 1 4 2 2 4 2 \n"
   "Face4Patches: 1 0 3 2 0 3 \n"
   "Face5: 1 5 2 2 5 2 \n"
   "Face5Patches: 1 0 3 2 0 3 \n"
   "Face6: 1 6 2 2 6 2 \n"
   "Face6Patches: 1 0 3 2 0 3 \n"
   "Frame: 1 1 1 1 3 1 1 5 1 1 7 1 1 9 1 1 10 1 1 11 1 1 12 1 2 2 1 2 4 1 2 6 1 2 8 1 2 9 1 2 10 1 2 11 1 2 12 1 \n"
   "InnerPatches: \n"
   "Vertex1: 1 1 0 \n"
   "Vertex2: 2 2 0 \n"
   "Vertex3: 1 3 0 \n"
   "Vertex4: 2 4 0 \n"
   "Vertex5: 1 5 0 \n"
   "Vertex6: 2 6 0 \n"
   "Vertex7: 1 7 0 \n"
   "Vertex8: 2 8 0 \n"},

  {"<geometry sets=\"true\" ny=\"2\"/>", 3,
   "700 1 0 0\n"
   "3 0\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "0 0 0\n"
   "1 0 0\n"
   "0 0.5 0\n"
   "1 0.5 0\n"
   "0 0 1\n"
   "1 0 1\n"
   "0 0.5 1\n"
   "1 0.5 1\n"
   "700 1 0 0\n"
   "3 0\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "2 2\n"
   "0 0 1 1\n"
   "0 0.5 0\n"
   "1 0.5 0\n"
   "0 1 0\n"
   "1 1 0\n"
   "0 0.5 1\n"
   "1 0.5 1\n"
   "0 1 1\n"
   "1 1 1\n",
   "Boundary: 1 1 2 1 2 2 1 3 2 1 5 2 1 6 2 2 1 2 2 2 2 2 4 2 2 5 2 2 6 2 \n"
   "Corners: 1 1 0 1 2 0 1 5 0 1 6 0 2 3 0 2 4 0 2 7 0 2 8 0 \n"
   "Edge1: 1 1 1 2 1 1 \n"
   "Edge10: 2 10 1 \n"
   "Edge11: 1 11 1 \n"
   "Edge12: 2 12 1 \n"
   "Edge2: 1 2 1 2 2 1 \n"
   "Edge3: 1 3 1 2 3 1 \n"
   "Edge4: 1 4 1 2 4 1 \n"
   "Edge5: 1 5 1 \n"
   "Edge6: 1 6 1 \n"
   "Edge7: 2 7 1 \n"
   "Edge8: 2 8 1 \n"
   "Edge9: 1 9 1 \n"
   "Face1: 1 1 2 2 1 2 \n"
   "Face1Patches: 1 0 3 2 0 3 \n"
   "Face2: 1 2 2 2 2 2 \n"
   "Face2Patches: 1 0 3 2 0 3 \n"
   "Face3: 1 3 2 \n"
   "Face3Patches: 1 0 3 \n"
   "Face4: 2 4 2 \n"
   "Face4Patches: 2 0 3 \n"
   "Face5: 1 5 2 2 5 2 \n"
   "Face5Patches: 1 0 3 2 0 3 \n"
   "Face6: 1 6 2 2 6 2 \n"
   "Face6Patches: 1 0 3 2 0 3 \n"
   "Frame: 1 1 1 1 2 1 1 3 1 1 4 1 1 5 1 1 6 1 1 9 1 1 11 1 2 1 1 2 2 1 2 3 1 2 4 1 2 7 1 2 8 1 2 10 1 2 12 1 \n"
   "InnerPatches: \n"
   "Vertex1: 1 1 0 \n"
   "Vertex2: 1 2 0 \n"
   "Vertex3: 2 3 0 \n"
   "Vertex4: 2 4 0 \n"
   "Vertex5: 1 5 0 \n"
   "Vertex6: 1 6 0 \n"
   "Vertex7: 2 7 0 \n"
   "Vertex8: 2 8 0 \n"},

   {"<geometry sets=\"true\" nz=\"2\"/>", 3,
    "700 1 0 0\n"
    "3 0\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "0 0 0\n"
    "1 0 0\n"
    "0 1 0\n"
    "1 1 0\n"
    "0 0 0.5\n"
    "1 0 0.5\n"
    "0 1 0.5\n"
    "1 1 0.5\n"
    "700 1 0 0\n"
    "3 0\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "0 0 0.5\n"
    "1 0 0.5\n"
    "0 1 0.5\n"
    "1 1 0.5\n"
    "0 0 1\n"
    "1 0 1\n"
    "0 1 1\n"
    "1 1 1\n",
    "Boundary: 1 1 2 1 2 2 1 3 2 1 4 2 1 5 2 2 1 2 2 2 2 2 3 2 2 4 2 2 6 2 \n"
    "Corners: 1 1 0 1 2 0 1 3 0 1 4 0 2 5 0 2 6 0 2 7 0 2 8 0 \n"
    "Edge1: 1 1 1 \n"
    "Edge10: 1 10 1 \n"
    "Edge11: 2 11 1 \n"
    "Edge12: 2 12 1 \n"
    "Edge2: 1 2 1 \n"
    "Edge3: 2 3 1 \n"
    "Edge4: 2 4 1 \n"
    "Edge5: 1 5 1 2 5 1 \n"
    "Edge6: 1 6 1 2 6 1 \n"
    "Edge7: 1 7 1 2 7 1 \n"
    "Edge8: 1 8 1 2 8 1 \n"
    "Edge9: 1 9 1 \n"
    "Face1: 1 1 2 2 1 2 \n"
    "Face1Patches: 1 0 3 2 0 3 \n"
    "Face2: 1 2 2 2 2 2 \n"
    "Face2Patches: 1 0 3 2 0 3 \n"
    "Face3: 1 3 2 2 3 2 \n"
    "Face3Patches: 1 0 3 2 0 3 \n"
    "Face4: 1 4 2 2 4 2 \n"
    "Face4Patches: 1 0 3 2 0 3 \n"
    "Face5: 1 5 2 \n"
    "Face5Patches: 1 0 3 \n"
    "Face6: 2 6 2 \n"
    "Face6Patches: 2 0 3 \n"
    "Frame: 1 1 1 1 2 1 1 5 1 1 6 1 1 7 1 1 8 1 1 9 1 1 10 1 2 3 1 2 4 1 2 5 1 2 6 1 2 7 1 2 8 1 2 11 1 2 12 1 \n"
    "InnerPatches: \n"
    "Vertex1: 1 1 0 \n"
    "Vertex2: 1 2 0 \n"
    "Vertex3: 1 3 0 \n"
    "Vertex4: 1 4 0 \n"
    "Vertex5: 2 5 0 \n"
    "Vertex6: 2 6 0 \n"
    "Vertex7: 2 7 0 \n"
    "Vertex8: 2 8 0 \n"},

  {"<geometry sets=\"true\" nx=\"2\" ny=\"2\" nz=\"2\"/>", 3,
    "700 1 0 0\n"
    "3 0\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "0 0 0\n"
    "0.5 0 0\n"
    "0 0.5 0\n"
    "0.5 0.5 0\n"
    "0 0 0.5\n"
    "0.5 0 0.5\n"
    "0 0.5 0.5\n"
    "0.5 0.5 0.5\n"
    "700 1 0 0\n"
    "3 0\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "0.5 0 0\n"
    "1 0 0\n"
    "0.5 0.5 0\n"
    "1 0.5 0\n"
    "0.5 0 0.5\n"
    "1 0 0.5\n"
    "0.5 0.5 0.5\n"
    "1 0.5 0.5\n"
    "700 1 0 0\n"
    "3 0\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "0 0.5 0\n"
    "0.5 0.5 0\n"
    "0 1 0\n"
    "0.5 1 0\n"
    "0 0.5 0.5\n"
    "0.5 0.5 0.5\n"
    "0 1 0.5\n"
    "0.5 1 0.5\n"
    "700 1 0 0\n"
    "3 0\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "0.5 0.5 0\n"
    "1 0.5 0\n"
    "0.5 1 0\n"
    "1 1 0\n"
    "0.5 0.5 0.5\n"
    "1 0.5 0.5\n"
    "0.5 1 0.5\n"
    "1 1 0.5\n"
    "700 1 0 0\n"
    "3 0\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "0 0 0.5\n"
    "0.5 0 0.5\n"
    "0 0.5 0.5\n"
    "0.5 0.5 0.5\n"
    "0 0 1\n"
    "0.5 0 1\n"
    "0 0.5 1\n"
    "0.5 0.5 1\n"
    "700 1 0 0\n"
    "3 0\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "0.5 0 0.5\n"
    "1 0 0.5\n"
    "0.5 0.5 0.5\n"
    "1 0.5 0.5\n"
    "0.5 0 1\n"
    "1 0 1\n"
    "0.5 0.5 1\n"
    "1 0.5 1\n"
    "700 1 0 0\n"
    "3 0\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "0 0.5 0.5\n"
    "0.5 0.5 0.5\n"
    "0 1 0.5\n"
    "0.5 1 0.5\n"
    "0 0.5 1\n"
    "0.5 0.5 1\n"
    "0 1 1\n"
    "0.5 1 1\n"
    "700 1 0 0\n"
    "3 0\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "2 2\n"
    "0 0 1 1\n"
    "0.5 0.5 0.5\n"
    "1 0.5 0.5\n"
    "0.5 1 0.5\n"
    "1 1 0.5\n"
    "0.5 0.5 1\n"
    "1 0.5 1\n"
    "0.5 1 1\n"
    "1 1 1\n",
    "Boundary: 1 1 2 1 3 2 1 5 2 "
    "2 2 2 2 3 2 2 5 2 "
    "3 1 2 3 4 2 3 5 2 "
    "4 2 2 4 4 2 4 5 2 "
    "5 1 2 5 3 2 5 6 2 "
    "6 2 2 6 3 2 6 6 2 "
    "7 1 2 7 4 2 7 6 2 "
    "8 2 2 8 4 2 8 6 2 \n"
    "Corners: 1 1 0 2 2 0 3 3 0 4 4 0 5 5 0 6 6 0 7 7 0 8 8 0 \n"
    "Edge1: 1 1 1 3 1 1 \n"
    "Edge10: 3 10 1 4 10 1 \n"
    "Edge11: 5 11 1 6 11 1 \n"
    "Edge12: 7 12 1 8 12 1 \n"
    "Edge2: 2 2 1 4 2 1 \n"
    "Edge3: 5 3 1 7 3 1 \n"
    "Edge4: 6 4 1 8 4 1 \n"
    "Edge5: 1 5 1 5 5 1 \n"
    "Edge6: 2 6 1 6 6 1 \n"
    "Edge7: 3 7 1 7 7 1 \n"
    "Edge8: 4 8 1 8 8 1 \n"
    "Edge9: 1 9 1 2 9 1 \n"
    "Face1: 1 1 2 3 1 2 5 1 2 7 1 2 \n"
    "Face1Patches: 1 0 3 3 0 3 5 0 3 7 0 3 \n"
    "Face2: 2 2 2 4 2 2 6 2 2 8 2 2 \n"
    "Face2Patches: 2 0 3 4 0 3 6 0 3 8 0 3 \n"
    "Face3: 1 3 2 2 3 2 5 3 2 6 3 2 \n"
    "Face3Patches: 1 0 3 2 0 3 5 0 3 6 0 3 \n"
    "Face4: 3 4 2 4 4 2 7 4 2 8 4 2 \n"
    "Face4Patches: 3 0 3 4 0 3 7 0 3 8 0 3 \n"
    "Face5: 1 5 2 2 5 2 3 5 2 4 5 2 \n"
    "Face5Patches: 1 0 3 2 0 3 3 0 3 4 0 3 \n"
    "Face6: 5 6 2 6 6 2 7 6 2 8 6 2 \n"
    "Face6Patches: 5 0 3 6 0 3 7 0 3 8 0 3 \n"
    "Frame: 1 1 1 1 5 1 1 9 1 "
    "2 2 1 2 6 1 2 9 1 "
    "3 1 1 3 7 1 3 10 1 "
    "4 2 1 4 8 1 4 10 1 "
    "5 3 1 5 5 1 5 11 1 "
    "6 4 1 6 6 1 6 11 1 "
    "7 3 1 7 7 1 7 12 1 "
    "8 4 1 8 8 1 8 12 1 \n"
    "InnerPatches: \n"
    "Vertex1: 1 1 0 \n"
    "Vertex2: 2 2 0 \n"
    "Vertex3: 3 3 0 \n"
    "Vertex4: 4 4 0 \n"
    "Vertex5: 5 5 0 \n"
    "Vertex6: 6 6 0 \n"
    "Vertex7: 7 7 0 \n"
    "Vertex8: 8 8 0 \n"}};


INSTANTIATE_TEST_CASE_P(TestMultiPatchModelGenerator3D,
                        TestMultiPatchModelGenerator3D,
                        testing::ValuesIn(geometry3D));
