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
  ASSERT_EQ(a(2,2), 7);
  ASSERT_EQ(a(3,2), 10);
  ASSERT_EQ(a(2,3), 14);
  ASSERT_EQ(a(3,3), 17);

  a.addBlock(b, 1, 1, 1, true);
  ASSERT_EQ(a(1,1), 2);
  ASSERT_EQ(a(2,1), 5);
  ASSERT_EQ(a(1,2), 6);
  ASSERT_EQ(a(2,2), 11);
}


TEST(TestMatrix, AddRows)
{
  utl::matrix<int> a(3,5);
  std::iota(a.begin(),a.end(),1);
  std::cout <<"A:"<< a;

  a.expandRows(1);
  std::cout <<"B:"<< a;
  int fasit = 1;
  for (size_t j = 1; j <= 5; j++)
  {
    for (size_t i = 1; i <= 3; i++, fasit++)
      ASSERT_EQ(a(i,j), fasit);
    ASSERT_EQ(a(4,j), 0);
  }

  a.expandRows(-2);
  std::cout <<"C:"<< a;
  fasit = 1;
  for (size_t j = 1; j <= 5; j++, fasit++)
    for (size_t i = 1; i <= 2; i++, fasit++)
      ASSERT_EQ(a(i,j), fasit);
}


TEST(TestMatrix3D, DumpRead)
{
  size_t n0(2), n1(3), n2(4);
  utl::matrix3d<double> A(n0, n1, n2);
  Real* p = A.ptr();
  for (size_t i = 0; i<A.size(); ++i, ++p)
    *p = 3.14159*i*i;

  std::string fname = std::tmpnam(nullptr);
  std::ofstream os(fname);
  os << std::setprecision(16) << A;
  std::ifstream is(fname, std::ios::in);
  utl::matrix3d<double> B(is);
  B -= A;
  ASSERT_NEAR(B.norm2(), 0.0, 1.0e-13);
}
