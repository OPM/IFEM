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


TEST(TestMatrix, AddBlock)
{
    utl::matrix<int> a(3,3);
    int data_a[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    a.fill(data_a, 9);

    utl::matrix<int> b(2,2);
    int data_b[9] = {1, 2, 3, 4};
    b.fill(data_b, 4);

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
