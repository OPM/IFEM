// $Id$
//==============================================================================
//!
//! \file StbImage.h
//!
//! \date Mar 7 2019
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Singleton wrapper for stb_image.
//!
//==============================================================================

#ifndef STBIMAGE_H_
#define STBIMAGE_H_

//! \brief Namespace for STB routines.
namespace stb {
  //! \brief Load an image from file.
  //! \param[in] file The file to load
  //! \param[out] width The width of the image
  //! \param[out] height The height of the image
  //! \param[out] nrChannels The number of channels in the image
  unsigned char* loadImage(const char* file, int& width, int& height, int& nrChannels);
}

#endif
