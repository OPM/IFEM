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

#include "MatrixTests.h"

#include <numeric>
#include <iomanip>
#include <fstream>


TEST(TestVector, Add)
{
  vectorAddTest<double>();
  vectorAddTest<float>();
}


TEST(TestVector, Dot)
{
  vectorDotTest<double>();
  vectorDotTest<float>();
}


TEST(TestVector, Multiply)
{
  vectorMultiplyTest<double>();
  vectorMultiplyTest<float>();
}


TEST(TestVector, Norm)
{
  vectorMultiplyTest<double>();
  vectorMultiplyTest<float>();
}


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


TEST(TestMatrix, ExtractBlock)
{
  utl::matrix<int> a(3,3), b(2,2);

  std::iota(a.begin(), a.end(), 1);

  a.extractBlock(b,1,1);
  EXPECT_EQ(b(1,1), 1);
  EXPECT_EQ(b(2,1), 2);
  EXPECT_EQ(b(1,2), 4);
  EXPECT_EQ(b(2,2), 5);

  a.extractBlock(b,2,2,true);
  EXPECT_EQ(b(1,1), 1+5);
  EXPECT_EQ(b(2,1), 2+6);
  EXPECT_EQ(b(1,2), 4+8);
  EXPECT_EQ(b(2,2), 5+9);
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


TEST(TestMatrix, AugmentRows)
{
  utl::matrix<int> a(5,3), b(4,3), c(3,2);
  std::iota(a.begin(),a.end(),1);
  std::iota(b.begin(),b.end(),16);
  std::cout <<"A:"<< a;
  std::cout <<"B:"<< b;
  size_t nA = a.size();
  size_t na = a.rows();
  size_t nb = b.rows();
  ASSERT_TRUE(a.augmentRows(b));
  ASSERT_FALSE(a.augmentRows(c));
  std::cout <<"C:"<< a;
  for (size_t j = 1; j <= a.cols(); j++)
    for (size_t i = 1; i <= a.rows(); i++)
      if (i <= na)
        EXPECT_EQ(a(i,j), i+na*(j-1));
      else
        EXPECT_EQ(a(i,j), nA-na+i+nb*(j-1));
}


TEST(TestMatrix, AugmentCols)
{
  utl::matrix<int> a(3,5), b(3,4), c(2,3);
  std::iota(a.begin(),a.end(),1);
  std::iota(b.begin(),b.end(),16);
  std::cout <<"A:"<< a;
  std::cout <<"B:"<< b;
  ASSERT_TRUE(a.augmentCols(b));
  ASSERT_FALSE(a.augmentCols(c));
  std::cout <<"C:"<< a;
  int fasit = 1;
  for (size_t j = 1; j <= a.cols(); j++)
    for (size_t i = 1; i <= a.rows(); i++, fasit++)
      EXPECT_EQ(a(i,j), fasit);
}


TEST(TestMatrix, SumCols)
{
  utl::matrix<int> a(5,3);
  std::iota(a.begin(),a.end(),1);
  std::cout <<"A:"<< a;
  EXPECT_EQ(a.sum(-1), 15);
  EXPECT_EQ(a.sum(-2), 40);
  EXPECT_EQ(a.sum(-3), 65);
}


TEST(TestMatrix, Multiply)
{
  multiplyTest<double>();
  multiplyTest<float>();
}


TEST(TestMatrix, Norm)
{
  normTest<double>();
  normTest<float>();
}


TEST(TestMatrix, OuterProduct)
{
  outerProductTest<double>();
  outerProductTest<float>();
}


TEST(TestMatrix, Read)
{
  utl::vector<double> a(26);
  utl::matrix<double> A(2,3);
  std::iota(a.begin(),a.end(),1.0);
  std::iota(A.begin(),A.end(),1.0);
  std::cout <<"a:"<< a;
  std::cout <<"A:"<< A;

  auto&& checkVector = [&a](const char* fname)
  {
    std::ifstream is(fname,std::ios::in);
    utl::vector<double> b;
    is >> b;
    std::cout <<"b:"<< b;
    ASSERT_EQ(a.size(),b.size());
    for (size_t i = 1; i <= a.size(); i++)
      EXPECT_NEAR(a(i), b(i), 1.0e-13);
  };

  auto&& checkMatrix = [&A](const char* fname)
  {
    std::ifstream is(fname,std::ios::in);
    utl::matrix<double> B;
    is >> B;
    std::cout <<"B:"<< B;
    ASSERT_EQ(A.rows(),B.rows());
    ASSERT_EQ(A.cols(),B.cols());
    for (size_t i = 1; i <= A.rows(); i++)
      for (size_t j = 1; j <= A.cols(); j++)
        EXPECT_NEAR(A(i,j), B(i,j), 1.0e-13);
  };

  const char* fname0 = "/tmp/testVector.dat";
  const char* fname1 = "/tmp/testMatrix1.dat";
  const char* fname2 = "/tmp/testMatrix2.dat";
  const char* fname3 = "/tmp/testMatrix3.dat";
  const char* fname4 = "/tmp/testMatrix4.dat";

  std::ofstream os(fname0);
  os << a.size() << a;
  os.close();

  os.open(fname1);
  os << A.rows() <<' '<< A.cols() << A;
  os.close();

  os.open(fname2);
  os << A.rows() <<' '<< A.cols();
  for (size_t i = 1; i <= A.rows(); i++)
    for (size_t j = 1; j <= A.cols(); j++)
      os << (j == 1 ? '\n' : ' ') << A(i,j);
  os <<'\n';
  os.close();

  checkVector(fname0);
  checkMatrix(fname1);
  checkMatrix(fname2);

  double value = 0.0;
  A.resize(6,6);
  for (size_t i = 1; i <= A.rows(); i++)
    for (size_t j = i; j <= A.cols(); j++)
      A(i,j) = A(j,i) = ++value;
  std::cout <<"Symmetric A:"<< A;

  os.open(fname3);
  os <<"Symmetric: "<< A.rows() << A;
  os.close();

  checkMatrix(fname3);

  std::iota(A.begin(),A.end(),1.0);
  std::cout <<"Non-symmetric A:"<< A;
  os.open(fname4);
  os <<"Column-oriented: "<< A.rows() <<" "<< A.cols();
  for (double v : A) os <<"\n"<< v;
  os <<"\n";
  os.close();

  checkMatrix(fname4);
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
  matrix3DMultiplyTest<double>();
  matrix3DMultiplyTest<float>();
}
