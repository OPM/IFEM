// $Id$
//==============================================================================
//!
//! \file StbImage.C
//!
//! \date Mar 7 2019
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Singleton wrapper for stb_image.
//!
//==============================================================================

#include "StbImage.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

unsigned char* stb::loadImage(const char* file, int& width,
                              int& height, int& nrChannels)
{
  return stbi_load(file,&width,&height,&nrChannels,0);
}
