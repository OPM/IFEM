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

#include <catch2/catch_template_test_macros.hpp>

#include <iomanip>
#include <numeric>
#include <fstream>


TEMPLATE_TEST_CASE("TestVector.Add", "", float, double)
{
  vectorAddTest<TestType>();
}


TEMPLATE_TEST_CASE("TestVector.Dot", "", float, double)
{
  vectorDotTest<TestType>();
}


TEMPLATE_TEST_CASE("TestVector.Multiply", "", float, double)
{
  vectorMultiplyTest<TestType>();
}


TEMPLATE_TEST_CASE("TestVector.Norm", "", float, double)
{
  vectorNormTest<TestType>();
}


TEST_CASE("TestMatrix.AddBlock")
{
  utl::matrix<int> a(3,3),b(2,2);
  std::iota(a.begin(),a.end(),1);
  std::iota(b.begin(),b.end(),1);

  a.addBlock(b, 2, 2, 2, false);
  REQUIRE(a(2,2) == 7);
  REQUIRE(a(3,2) == 10);
  REQUIRE(a(2,3) == 14);
  REQUIRE(a(3,3) == 17);

  a.addBlock(b, 1, 1, 1, true);
  REQUIRE(a(1,1) == 2);
  REQUIRE(a(2,1) == 5);
  REQUIRE(a(1,2) == 6);
  REQUIRE(a(2,2) == 11);
}


TEST_CASE("TestMatrix.ExtractBlock")
{
  utl::matrix<int> a(3,3), b(2,2);

  std::iota(a.begin(), a.end(), 1);

  a.extractBlock(b,1,1);
  REQUIRE(b(1,1) == 1);
  REQUIRE(b(2,1) == 2);
  REQUIRE(b(1,2) == 4);
  REQUIRE(b(2,2) == 5);

  a.extractBlock(b,2,2,true);
  REQUIRE(b(1,1) == 1+5);
  REQUIRE(b(2,1) == 2+6);
  REQUIRE(b(1,2) == 4+8);
  REQUIRE(b(2,2) == 5+9);
}


TEST_CASE("TestMatrix.AddRows")
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
      REQUIRE(a(i,j) == fasit);
    REQUIRE(a(4,j) == 0);
  }

  a.expandRows(-2);
  std::cout <<"C:"<< a;
  fasit = 1;
  for (size_t j = 1; j <= a.cols(); j++, fasit++)
    for (size_t i = 1; i <= 2; i++, fasit++)
      REQUIRE(a(i,j) == fasit);
}


TEST_CASE("TestMatrix.AugmentRows")
{
  utl::matrix<int> a(5,3), b(4,3), c(3,2);
  std::iota(a.begin(),a.end(),1);
  std::iota(b.begin(),b.end(),16);
  std::cout <<"A:"<< a;
  std::cout <<"B:"<< b;
  size_t nA = a.size();
  size_t na = a.rows();
  size_t nb = b.rows();
  REQUIRE(a.augmentRows(b));
  REQUIRE(!a.augmentRows(c));
  std::cout <<"C:"<< a;
  for (size_t j = 1; j <= a.cols(); j++)
    for (size_t i = 1; i <= a.rows(); i++)
      if (i <= na)
        REQUIRE(a(i,j) == static_cast<int>(i+na*(j-1)));
      else
        REQUIRE(a(i,j) == static_cast<int>(nA-na+i+nb*(j-1)));
}


TEST_CASE("TestMatrix.AugmentCols")
{
  utl::matrix<int> a(3,5), b(3,4), c(2,3);
  std::iota(a.begin(),a.end(),1);
  std::iota(b.begin(),b.end(),16);
  std::cout <<"A:"<< a;
  std::cout <<"B:"<< b;
  REQUIRE(a.augmentCols(b));
  REQUIRE(!a.augmentCols(c));
  std::cout <<"C:"<< a;
  int fasit = 1;
  for (size_t j = 1; j <= a.cols(); j++)
    for (size_t i = 1; i <= a.rows(); i++, fasit++)
      REQUIRE(a(i,j) == fasit);
}


TEST_CASE("TestMatrix.SumCols")
{
  utl::matrix<int> a(5,3);
  std::iota(a.begin(),a.end(),1);
  std::cout <<"A:"<< a;
  REQUIRE(a.sum(-1) == 15);
  REQUIRE(a.sum(-2) == 40);
  REQUIRE(a.sum(-3) == 65);
}


TEST_CASE("TestMatrix.Fill")
{
  utl::matrix<int> a;
  std::vector<int> v(16);
  std::iota(v.begin(),v.end(),1);
  a.fill(v,3,4);
  std::cout <<"a:"<< a;
  REQUIRE(a(3,1) == 3);
  a.fill(v,4,4);
  std::cout <<"a:"<< a;
  REQUIRE(a(4,1) == 4);
  a.fill(v,5,4);
  std::cout <<"a:"<< a;
  REQUIRE(a(5,1) == 0);
}


TEST_CASE("TestMatrix.Zero")
{
  utl::matrix<double> A(4,5);
  REQUIRE(A.zero());
  A(1,2) = 1.0e-8;
  REQUIRE(!A.zero());
  REQUIRE(A.zero(1.0e-6));
}


TEMPLATE_TEST_CASE("TestMatrix.Multiply", "", float, double)
{
  multiplyTest<TestType>();
}


TEMPLATE_TEST_CASE("TestMatrix.Norm", "", float, double)
{
  normTest<TestType>();
}


TEMPLATE_TEST_CASE("TestMatrix.OuterProduct", "", float, double)
{
  outerProductTest<TestType>();
}


TEST_CASE("TestMatrix.Read")
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
    REQUIRE(a.size() == b.size());
    for (size_t i = 1; i <= a.size(); i++)
      REQUIRE_THAT(a(i), WithinRel(b(i), 1.0e-13));
  };

  auto&& checkMatrix = [&A](const char* fname)
  {
    std::ifstream is(fname,std::ios::in);
    utl::matrix<double> B;
    is >> B;
    std::cout <<"B:"<< B;
    REQUIRE(A.rows() == B.rows());
    REQUIRE(A.cols() == B.cols());
    for (size_t i = 1; i <= A.rows(); i++)
      for (size_t j = 1; j <= A.cols(); j++)
        REQUIRE_THAT(A(i,j), WithinRel(B(i,j), 1.0e-13));
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


TEST_CASE("TestMatrix3D.Trace")
{
  utl::matrix3d<double> a(4,3,3);
  std::iota(a.begin(),a.end(),1.0);
  std::cout <<"A:"<< a;

  for (size_t i = 1; i <= 4; i++)
    REQUIRE_THAT(a.trace(i), WithinRel(3.0*i+48.0));
}


TEST_CASE("TestMatrix3D.GetColumn")
{
  utl::matrix3d<int> A(4,3,2);
  std::iota(A.begin(),A.end(),1);
  std::cout <<"A:"<< A;

  int value = 1;
  for (size_t c = 1; c <= A.dim(3); c++)
    for (size_t r = 1; r <= A.dim(2); r++)
    {
      utl::vector<int> column = A.getColumn(r,c);
      REQUIRE(column.size() == A.dim(1));
      for (size_t i = 0; i < column.size(); i++, value++)
        REQUIRE(value == column[i]);
    }

  REQUIRE(value == static_cast<int>(1+A.size()));
}


TEST_CASE("TestMatrix3D.DumpRead")
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
  REQUIRE_THAT(B.norm2(), WithinAbs(0.0, 1.0e-13));
}


TEMPLATE_TEST_CASE("TestMatrix3D.Multiply", "", float, double)
{
  matrix3DMultiplyTest<TestType>();
}
