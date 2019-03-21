//==============================================================================
//!
//! \file TestMatrix.C
//!
//! \date Apr 11 2016
//!
//! \author Eivind Fonn / SINTEF
//!
//! \brief Unit tests for matrix and matrix3d.
//!
//==============================================================================

#include "MatVec.h"

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


TEST(TestMatrix, Multiply)
{
  utl::vector<double> u(14), v(9), x, y;
  utl::matrix<double> A(3,5);

  std::iota(u.begin(),u.end(),1.0);
  std::iota(v.begin(),v.end(),1.0);
  std::iota(A.begin(),A.end(),1.0);

  ASSERT_TRUE(A.multiply(u,x,1.0,0.0,false,3,4,1,2));
  ASSERT_TRUE(A.multiply(v,y,1.0,0.0,true,4,2));

  EXPECT_FLOAT_EQ(x(3),370.0);
  EXPECT_FLOAT_EQ(x(7),410.0);
  EXPECT_FLOAT_EQ(x(11),450.0);
  EXPECT_FLOAT_EQ(y(1),38.0);
  EXPECT_FLOAT_EQ(y(3),83.0);
  EXPECT_FLOAT_EQ(y(5),128.0);
  EXPECT_FLOAT_EQ(y(7),173.0);
  EXPECT_FLOAT_EQ(y(9),218.0);

  ASSERT_TRUE(A.multiply(u,x,1.0,-1.0,false,3,4,1,2));
  ASSERT_TRUE(A.multiply(v,y,1.0,-1.0,true,4,2));

  EXPECT_FLOAT_EQ(x.sum(),0.0);
  EXPECT_FLOAT_EQ(y.sum(),0.0);

  u.std::vector<double>::resize(5);
  v = 0.5*A*u;
  EXPECT_FLOAT_EQ(v(1),67.5);
  EXPECT_FLOAT_EQ(v(2),75.0);
  EXPECT_FLOAT_EQ(v(3),82.5);
}


TEST(TestMatrix, Norm)
{
  utl::matrix<double> a(4,5);
  std::iota(a.begin(),a.end(),1.0);
  std::cout <<"A:"<< a;

  EXPECT_FLOAT_EQ(a.sum(),210.0);
  EXPECT_FLOAT_EQ(a.sum(5),34.0);
  EXPECT_FLOAT_EQ(a.asum(5),34.0);
  EXPECT_FLOAT_EQ(a.trace(),34.0);
  EXPECT_NEAR(a.norm2(5),sqrt(414.0),1.0e-15);
}


TEST(TestMatrix3D, Trace)
{
  utl::matrix3d<double> a(4,3,3);
  std::iota(a.begin(),a.end(),1.0);
  std::cout <<"A:"<< a;

  for (size_t i = 1; i <= 4; i++)
    EXPECT_FLOAT_EQ(a.trace(i),3.0*i+48.0);
}


TEST(TestMatrix3D, GetColumn)
{
  utl::matrix3d<int> A(4,3,2);
  std::iota(A.begin(),A.end(),1);
  std::cout <<"A:"<< A;

  int value = 1;
  for (size_t c = 1; c <= A.dim(3); c++)
    for (size_t r = 1; r <= A.dim(2); r++)
    {
      utl::vector<int> column = A.getColumn(r,c);
      EXPECT_EQ(column.size(), A.dim(1));
      for (size_t i = 0; i < column.size(); i++, value++)
        EXPECT_EQ(value, column[i]);
    }

  EXPECT_EQ(value, 1+(int)A.size());
}


TEST(TestMatrix3D, DumpRead)
{
  int i = 0;
  utl::matrix3d<double> A(2,3,4);
  for (double& v : A) v = 3.14159*(++i);

  const char* fname = "/tmp/testMatrix3D.dat";
  std::ofstream os(fname);
  os << std::setprecision(16) << A;
  std::ifstream is(fname,std::ios::in);
  utl::matrix3d<double> B(is);
  B -= A;
  ASSERT_NEAR(B.norm2(), 0.0, 1.0e-13);
}


TEST(TestMatrix3D, Multiply)
{
  std::vector<double> a(10);
  utl::matrix<double> A(2,5);
  utl::matrix3d<double> B(5,4,3), C, D;

  std::iota(a.begin(),a.end(),1.0);
  std::iota(A.begin(),A.end(),1.0);
  std::iota(B.begin(),B.end(),1.0);

  C.multiply(A,B);
  ASSERT_TRUE(D.multiplyMat(a,B));

  std::vector<double>::const_iterator c = C.begin();
  for (double d : D)
    EXPECT_EQ(d,*(c++));
}
