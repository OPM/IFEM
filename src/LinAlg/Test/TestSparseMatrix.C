//==============================================================================
//!
//! \file TestSparseMatrix.C
//!
//! \date Apr 9 2019
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Unit tests for sparse matrices
//!
//==============================================================================

#include "SparseMatrix.h"

#include "gtest/gtest.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif


TEST(TestSparseMatrix, CalcCSR)
{
  SparseMatrix Mat1(2,3);
  Mat1(1,1) = 1.0;
  Mat1(2,1) = 2.0;
  Mat1(2,3) = 3.0;
  IntVec IA1, JA1;
  Mat1.calcCSR(IA1, JA1);

  EXPECT_EQ(IA1.size(), 3U);
  EXPECT_EQ(IA1[0], 0);
  EXPECT_EQ(IA1[1], 1);

  EXPECT_EQ(JA1.size(), 3U);
  EXPECT_EQ(JA1[0], 0);
  EXPECT_EQ(JA1[1], 0);
  EXPECT_EQ(JA1[2], 2);
}


TEST(TestSparseMatrix, split)
{
  SparseMatrix Amat(10,10);

  size_t i, j;
  for (i = 1; i <= 10; i++)
    Amat(i,i) = i*10;

  Amat(1,2) = Amat(2,1) = 1.0;
  Amat(1,3) = Amat(3,1) = 2.0;
  Amat(1,6) = Amat(6,1) = 3.0;
  Amat(2,4) = Amat(4,2) = 4.0;
  Amat(2,8) = Amat(8,2) = 5.0;
  Amat(3,5) = Amat(5,3) = 6.0;
  Amat(4,7) = Amat(7,4) = 7.0;
  Amat(5,7) = Amat(7,5) = 8.0;
  Amat(6,7) = Amat(7,6) = 9.0;

  Amat.dump(std::cout,LinAlg::MATRIX_MARKET,"A_mmarket");
  Amat.dump(std::cout,LinAlg::MATLAB,"A_matlab");
  Amat.dump(std::cout,LinAlg::FLAT,"A");

  std::array<SparseMatrix,4> Asub;
  ASSERT_TRUE(Amat.split(Asub,{2,5,7}));
  std::cout <<"A11: "<< Asub[0];
  std::cout <<"A21: "<< Asub[1];
  std::cout <<"A12: "<< Asub[2];
  std::cout <<"A22: "<< Asub[3];

  EXPECT_EQ(Asub[0](1,1),10.0);
  EXPECT_EQ(Asub[0](2,2),30.0);
  EXPECT_EQ(Asub[0](3,3),40.0);
  EXPECT_EQ(Asub[0](4,4),60.0);
  EXPECT_EQ(Asub[0](5,5),80.0);
  EXPECT_EQ(Asub[0](6,6),90.0);
  EXPECT_EQ(Asub[0](7,7),100.0);

  EXPECT_EQ(Asub[3](1,1),20.0);
  EXPECT_EQ(Asub[3](2,2),50.0);
  EXPECT_EQ(Asub[3](3,3),70.0);

  EXPECT_EQ(Asub[0](1,2),2.0);
  EXPECT_EQ(Asub[0](1,4),3.0);
  EXPECT_EQ(Asub[2](1,1),1.0);
  EXPECT_EQ(Asub[2](2,2),6.0);
  EXPECT_EQ(Asub[2](3,1),4.0);
  EXPECT_EQ(Asub[2](3,3),7.0);
  EXPECT_EQ(Asub[2](4,3),9.0);
  EXPECT_EQ(Asub[2](5,1),5.0);
  EXPECT_EQ(Asub[3](2,3),8.0);

  for (i = 1; i <= Asub[0].dim(1); i++)
    for (j = 1; j < i; j++)
      EXPECT_EQ(Asub[0](i,j),Asub[0](j,i));

  for (i = 1; i <= Asub[1].dim(1); i++)
    for (j = 1; j <= Asub[1].dim(1); j++)
      EXPECT_EQ(Asub[1](i,j),Asub[2](j,i));

  for (i = 1; i <= Asub[3].dim(1); i++)
    for (j = 1; j < i; j++)
      EXPECT_EQ(Asub[3](i,j),Asub[3](j,i));

  Vector aVec;
  ASSERT_TRUE(Asub[2].getColumn(3,aVec));
  EXPECT_EQ(aVec(1),0.0);
  EXPECT_EQ(aVec(2),0.0);
  EXPECT_EQ(aVec(3),7.0);
  EXPECT_EQ(aVec(4),9.0);
  EXPECT_EQ(aVec(5),0.0);
  EXPECT_EQ(aVec(6),0.0);
  EXPECT_EQ(aVec(7),0.0);
}

class TestSparseMatrix : public testing::Test,
                         public testing::WithParamInterface<SparseMatrix::SparseSolver>
{
};


TEST_P(TestSparseMatrix, MxV)
{
  constexpr size_t size = 10;
  SparseMatrix A(size, size, GetParam());
  StdVector x(size);
  StdVector y(size);

  for (size_t i = 1; i <= size; ++i) {
    if (i > 1)
      A(i, i-1) = 1.0;
    A(i, i) = -4.0;
    if (i < size)
      A(i, i+1) = 2.0;
    x(i) = i;
  }

  auto checkVec = [&y]()
  {
    for (size_t i = 1; i <= size; ++i)
      EXPECT_DOUBLE_EQ(y(i), (i > 1 ? i-1 : 0.0) -
                             (4.0 * i) +
                             (i < size ? 2.0*(i+1) : 0.0));
  };

  A.multiply(x, y);
  checkVec();

  A.compressPattern();

#ifdef USE_OPENMP
  omp_set_num_threads(1);
#endif
  A.multiply(x, y);
  checkVec();

#ifdef USE_OPENMP
  omp_set_num_threads(2);
  A.multiply(x, y);
  checkVec();
#endif
}

INSTANTIATE_TEST_SUITE_P(TestSparseMatrix, TestSparseMatrix,
                         testing::Values(SparseMatrix::SUPERLU, SparseMatrix::S_A_M_G));
