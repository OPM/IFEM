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

#include "Catch2Support.h"

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
  CHECK(a(2,2) == 7);
  CHECK(a(3,2) == 10);
  CHECK(a(2,3) == 14);
  CHECK(a(3,3) == 17);

  a.addBlock(b, 1, 1, 1, true);
  CHECK(a(1,1) == 2);
  CHECK(a(2,1) == 5);
  CHECK(a(1,2) == 6);
  CHECK(a(2,2) == 11);
}


TEST_CASE("TestMatrix.ExtractBlock")
{
  utl::matrix<int> a(3,3), b(2,2);

  std::iota(a.begin(), a.end(), 1);

  a.extractBlock(b,1,1);
  CHECK(b(1,1) == 1);
  CHECK(b(2,1) == 2);
  CHECK(b(1,2) == 4);
  CHECK(b(2,2) == 5);

  a.extractBlock(b,2,2,true);
  CHECK(b(1,1) == 1+5);
  CHECK(b(2,1) == 2+6);
  CHECK(b(1,2) == 4+8);
  CHECK(b(2,2) == 5+9);
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
      CHECK(a(i,j) == fasit);
    CHECK(a(4,j) == 0);
  }

  a.expandRows(-2);
  std::cout <<"C:"<< a;
  fasit = 1;
  for (size_t j = 1; j <= a.cols(); j++, fasit++)
    for (size_t i = 1; i <= 2; i++, fasit++)
      CHECK(a(i,j) == fasit);

  a.expandRows(3,true);
  std::cout <<"D:"<< a;
  fasit = 1;
  for (size_t j = 1; j <= a.cols(); j++, fasit++)
  {
    for (size_t i = 1; i <= 2; i++, fasit++)
      CHECK(a(i,j) == fasit);
    CHECK(a(3,j) == 0);
  }

  a.expandRows(1,true);
  std::cout <<"E:"<< a;
  fasit = 1;
  for (size_t j = 1; j <= a.cols(); j++, fasit += 3)
    CHECK(a(1,j) == fasit);
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
        CHECK(a(i,j) == static_cast<int>(i+na*(j-1)));
      else
        CHECK(a(i,j) == static_cast<int>(nA-na+i+nb*(j-1)));
  REQUIRE(a.augmentRows(b,true));
  std::cout <<"C:"<< a;
  for (size_t j = 1; j <= a.cols(); j++)
    for (size_t i = 1; i <= a.rows(); i++)
      if (i <= nb)
        CHECK(a(i,j) == b(i,j));
      else if (i > nb+na)
        CHECK(a(i,j) == b(i-nb-na,j));
      else
        CHECK(a(i,j) == static_cast<int>(i-nb+na*(j-1)));
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
      CHECK(a(i,j) == fasit);
}


TEST_CASE("TestMatrix.SumCols")
{
  utl::matrix<int> a(5,3);
  std::iota(a.begin(),a.end(),1);
  std::cout <<"A:"<< a;
  CHECK(a.sum(-1) == 15);
  CHECK(a.sum(-2) == 40);
  CHECK(a.sum(-3) == 65);

  int fasit = 15;
  for (size_t i = 1; i <= a.cols(); i++, fasit += 25)
    CHECK(a.colsum(i) == fasit);

  fasit = 18;
  for (size_t i = 1; i <= a.rows(); i++, fasit += a.cols())
    CHECK(a.rowsum(i) == fasit);
}


TEST_CASE("TestMatrix.Fill")
{
  utl::matrix<int> a;
  std::vector<int> v(16);
  std::iota(v.begin(),v.end(),1);
  a.fill(v,3,4);
  std::cout <<"a:"<< a;
  CHECK(a(3,1) == 3);
  a.fill(v,4,4);
  std::cout <<"a:"<< a;
  CHECK(a(4,1) == 4);
  a.fill(v,5,4);
  std::cout <<"a:"<< a;
  CHECK(a(5,1) == 0);
}


TEST_CASE("TestMatrix.Zero")
{
  utl::matrix<double> A(4,5);
  CHECK(A.zero());
  A(1,2) = 1.0e-8;
  CHECK(!A.zero());
  CHECK(A.zero(1.0e-6));
}


TEST_CASE("TestMatrix.Transpose")
{
  utl::matrix<int> A(4,5);
  std::iota(A.begin(),A.end(),1);
  utl::matrix<int> B(A,true);
  std::cout <<"A:"<< A <<"B:"<< B;
  REQUIRE(A.rows() == B.cols());
  REQUIRE(A.cols() == B.rows());
  for (size_t i = 1; i <= A.rows(); i++)
    for (size_t j = 1; j <= A.cols(); j++)
      CHECK(A(i,j) == B(j,i));
}


TEST_CASE("TestMatrix.External")
{
  utl::vector<int> a(6);
  std::iota(a.begin(),a.end(),1);
  utl::matrix<int> A(a);
  A.resize(2,3);
  std::cout <<"A:"<< A;
  int value = 1;
  for (size_t j = 1; j <= A.cols(); j++)
    for (size_t i = 1; i <= A.rows(); i++, value++)
      CHECK(A(i,j) == value);
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
      CHECK_THAT(a(i), WithinRel(b(i), 1.0e-13));
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
        CHECK_THAT(A(i,j), WithinRel(B(i,j), 1.0e-13));
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
    CHECK_THAT(a.trace(i), WithinRel(3.0*i+48.0));
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
        CHECK(value == column[i]);
    }

  CHECK(value == static_cast<int>(1+A.size()));
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
  CHECK_THAT(B.norm2(), WithinAbs(0.0, 1.0e-13));
}


TEMPLATE_TEST_CASE("TestMatrix3D.Multiply", "", float, double)
{
  matrix3DMultiplyTest<TestType>();
}
