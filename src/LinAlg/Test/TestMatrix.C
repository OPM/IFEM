//==============================================================================
//!
//! \file TestMatrix.C
//!
//! \date Apr 11 2016
//!
//! \author Eivind Fonn / SINTEF
//!
//! \brief Unit tests for matrix
//!
//==============================================================================

#include "matrix.h"

#include "gtest/gtest.h"
#include <numeric>
#include <iomanip>
#include <fstream>


TEST(TestMatrix, AddBlock)
{
  utl::matrix<int> a(3,3),b(2,2);
  std::iota(a.begin(),a.end(),1);
  std::iota(b.begin(),b.end(),1);

  a.addBlock(b, 2, 2, 2, false);
  EXPECT_EQ(a(2,2), 7);
  EXPECT_EQ(a(3,2), 10);
  EXPECT_EQ(a(2,3), 14);
  EXPECT_EQ(a(3,3), 17);

  a.addBlock(b, 1, 1, 1, true);
  EXPECT_EQ(a(1,1), 2);
  EXPECT_EQ(a(2,1), 5);
  EXPECT_EQ(a(1,2), 6);
  EXPECT_EQ(a(2,2), 11);
}


TEST(TestMatrix, AddRows)
{
  utl::matrix<int> a(3,5);
  std::iota(a.begin(),a.end(),1);
  std::cout <<"A:"<< a;

  a.expandRows(1);
  std::cout <<"B:"<< a;
  int fasit = 1;
  for (size_t j = 1; j <= a.cols(); j++)
  {
    for (size_t i = 1; i <= 3; i++, fasit++)
      EXPECT_EQ(a(i,j), fasit);
    EXPECT_EQ(a(4,j), 0);
  }

  a.expandRows(-2);
  std::cout <<"C:"<< a;
  fasit = 1;
  for (size_t j = 1; j <= a.cols(); j++, fasit++)
    for (size_t i = 1; i <= 2; i++, fasit++)
      EXPECT_EQ(a(i,j), fasit);
}


TEST(TestMatrix, Norm)
{
  utl::matrix<double> a(4,5);
  std::iota(a.begin(),a.end(),1.0);
  std::cout <<"A:"<< a;

  EXPECT_FLOAT_EQ(a.sum(),210.0);
  EXPECT_FLOAT_EQ(a.sum(5),34.0);
  EXPECT_FLOAT_EQ(a.asum(5),34.0);
  EXPECT_NEAR(a.norm2(5),sqrt(414.0),1.0e-15);
}


TEST(TestMatrix3D, GetColumn)
{
  utl::matrix3d<int> A(4,3,2);
  std::iota(A.begin(),A.end(),1);
  std::cout <<"A:"<< A;

  int value = 0;
  for (size_t c = 1; c <= A.dim(3); c++)
    for (size_t r = 1; r <= A.dim(2); r++)
    {
      utl::vector<int> column = A.getColumn(r,c);
      EXPECT_EQ(column.size(),A.dim(1));
      for (int v : column)
        EXPECT_EQ(++value, v);
    }

  EXPECT_EQ(value, (int)A.size());
}


TEST(TestMatrix3D, DumpRead)
{
  utl::matrix3d<double> A(2,3,4);
  double i = 0.0;
  for (auto p = A.begin(); p != A.end(); ++p, ++i)
    *p = 3.14159*i*i;

  const char* fname = "/tmp/testMatrix3D.dat";
  std::ofstream os(fname);
  os << std::setprecision(16) << A;
  std::ifstream is(fname,std::ios::in);
  utl::matrix3d<double> B(is);
  B -= A;
  ASSERT_NEAR(B.norm2(), 0.0, 1.0e-13);
}
