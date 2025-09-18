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

#include <catch2/catch_test_macros.hpp>


TEST_CASE("TestImageReader.Square")
{
  int width, height, nrChannels;
  unsigned char* image = stbi_load("src/Utility/Test/refdata/chess.png",
                                   &width, &height, &nrChannels, 0);
  REQUIRE(width == 512);
  REQUIRE(height == 512);
  REQUIRE(nrChannels ==  4);
  int i, j, k = 0;
  for (j = 0; j < height; j++)
    for (i = 0; i < width; i++) {
      bool is_white = ((i/64)%2 == (j/64)%2);
      REQUIRE(image[k++] == is_white * 255); // red
      REQUIRE(image[k++] == is_white * 255); // green
      REQUIRE(image[k++] == is_white * 255); // blue
      REQUIRE(image[k++] == 255);            // alpha
    }
}


TEST_CASE("TestImageReader.Rectangular")
{
  int width, height, nrChannels;
  unsigned char* image = stbi_load("src/Utility/Test/refdata/rectangular_chess.png",
                                   &width, &height, &nrChannels, 0);
  REQUIRE(width == 256);
  REQUIRE(height == 128);
  REQUIRE(nrChannels == 4);
  int i, j, k = 0;
  for (j = 0; j < height; j++)
    for (i = 0; i < width; i++) {
      bool is_white = ((i/64)%2 == (j/64)%2);
      REQUIRE(image[k++] == is_white * 255); // red
      REQUIRE(image[k++] == is_white * 255); // green
      REQUIRE(image[k++] == is_white * 255); // blue
      REQUIRE(image[k++] == 255);            // alpha
    }
}
