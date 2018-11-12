//==============================================================================
//!
//! \file TestImageReader.C
//!
//! \date Sept 2018
//!
//! \author Kjetil A. Johannessen / SINTEF
//!
//! \brief Tests for image reader.
//!
//==============================================================================

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "gtest/gtest.h"


TEST(TestImageReader, square)
{
  int width, height, nrChannels;
  unsigned char* image = stbi_load("src/Utility/Test/refdata/chess.png",
                                   &width, &height, &nrChannels, 0);
  ASSERT_EQ(width,      512);
  ASSERT_EQ(height,     512);
  ASSERT_EQ(nrChannels,   4);
  int i, j, k = 0;
  for (j = 0; j < height; j++)
    for (i = 0; i < width; i++) {
      bool is_white = ((i/64)%2 == (j/64)%2);
      EXPECT_EQ(image[k++], is_white * 255); // red
      EXPECT_EQ(image[k++], is_white * 255); // green
      EXPECT_EQ(image[k++], is_white * 255); // blue
      EXPECT_EQ(image[k++], 255);            // alpha
    }
}


TEST(TestImageReader, rectangular)
{
  int width, height, nrChannels;
  unsigned char* image = stbi_load("src/Utility/Test/refdata/rectangular_chess.png",
                                   &width, &height, &nrChannels, 0);
  ASSERT_EQ(width,      256);
  ASSERT_EQ(height,     128);
  ASSERT_EQ(nrChannels,   4);
  int i, j, k = 0;
  for (j = 0; j < height; j++)
    for (i = 0; i < width; i++) {
      bool is_white = ((i/64)%2 == (j/64)%2);
      EXPECT_EQ(image[k++], is_white * 255); // red
      EXPECT_EQ(image[k++], is_white * 255); // green
      EXPECT_EQ(image[k++], is_white * 255); // blue
      EXPECT_EQ(image[k++], 255);            // alpha
    }
}
