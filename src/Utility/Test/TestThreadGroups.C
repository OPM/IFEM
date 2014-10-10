//==============================================================================
//!
//! \file TestThreadGroups.C
//!
//! \date Oct 13 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for threading group partitioning.
//!
//==============================================================================

#include "ThreadGroups.h"
#ifdef USE_OPENMP
#include <omp.h>
#endif
#include <fstream>

#include "gtest/gtest.h"

typedef std::vector< std::vector<int> > IntMat;

static IntMat readIntMatrix(size_t r, const std::string& file)
{
  std::vector< std::vector<int> > result;
  result.resize(r);
  std::ifstream f(file);
  for (size_t i=0;i<r;++i) {
    size_t size;
    f >> size;
    result[i].resize(size);
    for (size_t j=0;j<size;++j)
      f >> result[i][j];
  }

  return result;
}

#define CHECK_INTMATRICES_EQUAL(A,path) \
  do { \
  IntMat B = readIntMatrix(A.size(), path); \
  for (size_t i=0;i<A.size();++i) \
    for (size_t j=0;j<A[i].size();++j) \
      ASSERT_EQ(A[i][j], B[i][j]); \
  } while(0);

TEST(TestThreadGroups, Groups2D)
{
  ThreadGroups groups;
#ifdef USE_OPENMP
  omp_set_num_threads(2);
#endif

  groups.calcGroups(8, 8, 2);

#ifdef USE_OPENMP
  ASSERT_EQ(groups.size(), 2U);
  CHECK_INTMATRICES_EQUAL(groups[0], "src/Utility/Test/refdata/ThreadGroups_2D_p1.ref");
  CHECK_INTMATRICES_EQUAL(groups[1], "src/Utility/Test/refdata/ThreadGroups_2D_p2.ref");
#else
  ASSERT_EQ(groups.size(), 1U);
  CHECK_INTMATRICES_EQUAL(groups[0], "src/Utility/Test/refdata/ThreadGroups_2D_1.ref");
#endif

#ifdef USE_OPENMP
  omp_set_num_threads(3);
#endif
  ThreadGroups groups2;
  std::vector<bool> b1(12, true);
  b1[3] = false;
  std::vector<bool> b2(12, true);
  b1[6] = false;
  groups2.calcGroups(b1, b2, 2, 2);
#ifdef USE_OPENMP
  ASSERT_EQ(groups.size(), 2U);
  CHECK_INTMATRICES_EQUAL(groups2[0], "src/Utility/Test/refdata/ThreadGroups_2D_p1_empty.ref");
  CHECK_INTMATRICES_EQUAL(groups2[1], "src/Utility/Test/refdata/ThreadGroups_2D_p2_empty.ref");
#else
  ASSERT_EQ(groups.size(), 1U);
  CHECK_INTMATRICES_EQUAL(groups2[0], "src/Utility/Test/refdata/ThreadGroups_2D_1_empty.ref");
#endif
}

TEST(TestThreadGroups, Groups3D)
{
  ThreadGroups groups;
#ifdef USE_OPENMP
  omp_set_num_threads(2);
#endif

  groups.calcGroups(8, 8, 8, 2);

#ifdef USE_OPENMP
  ASSERT_EQ(groups.size(), 2U);
  CHECK_INTMATRICES_EQUAL(groups[0], "src/Utility/Test/refdata/ThreadGroups_3D_p1.ref");
  CHECK_INTMATRICES_EQUAL(groups[1], "src/Utility/Test/refdata/ThreadGroups_3D_p2.ref");
#else
  ASSERT_EQ(groups.size(), 1U);
  CHECK_INTMATRICES_EQUAL(groups[0], "src/Utility/Test/refdata/ThreadGroups_3D_1.ref");
#endif

#ifdef USE_OPENMP
  omp_set_num_threads(3);
#endif
  ThreadGroups groups2;
  std::vector<bool> b1(14, true);
  b1[3] = false;
  std::vector<bool> b2(14, true);
  b2[6] = false;
  std::vector<bool> b3(14, true);
  b3[7] = false;
  groups2.calcGroups(b1, b2, b3, 2, 2, 2);

#ifdef USE_OPENMP
  ASSERT_EQ(groups.size(), 2U);
  CHECK_INTMATRICES_EQUAL(groups2[0], "src/Utility/Test/refdata/ThreadGroups_3D_p1_empty.ref");
  CHECK_INTMATRICES_EQUAL(groups2[1], "src/Utility/Test/refdata/ThreadGroups_3D_p2_empty.ref");
#else
  ASSERT_EQ(groups.size(), 1U);
  CHECK_INTMATRICES_EQUAL(groups2[0], "src/Utility/Test/refdata/ThreadGroups_3D_1_empty.ref");
#endif
}
