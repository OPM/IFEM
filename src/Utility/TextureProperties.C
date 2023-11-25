// $Id$
//==============================================================================
//!
//! \file TextureProperties.C
//!
//! \date Jan 3 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Properties defined through a texture map.
//!
//==============================================================================

#include "TextureProperties.h"
#include "IFEM.h"
#include "HDF5Reader.h"
#include "ProcessAdm.h"
#include "Utilities.h"
#include "Vec3.h"

#include "tinyxml2.h"
#include "StbImage.h"


void TextureProperties::parse(const tinyxml2::XMLElement* elem)
{
  const tinyxml2::XMLElement* child = elem->FirstChildElement("property");
  for (; child; child = child->NextSiblingElement()) {
    std::string prop;
    utl::getAttribute(child,"name",prop);
    if (prop.empty()) {
      std::cerr << "No name for property, skipping.." << std::endl;
      continue;
    }

    std::string textureFile;
    utl::getAttribute(child, "file", textureFile);

    if (textureFile.find(".h5") != std::string::npos ||
        textureFile.find(".hdf5") != std::string::npos) {
      properties[prop].prescaled = true;
      ProcessAdm adm;
      HDF5Reader reader(textureFile, adm);
      reader.read3DArray(prop, properties[prop].textureData);
      properties[prop].min = *std::min_element(properties[prop].textureData.begin(),
                                               properties[prop].textureData.end());
      properties[prop].max = *std::max_element(properties[prop].textureData.begin(),
                                               properties[prop].textureData.end());
    } else {
      int width, height, nrChannels;
      unsigned char* image = stb::loadImage(textureFile.c_str(),
                                            width, height, nrChannels);
      if (!image) {
        std::cerr << "File not found: " << textureFile << std::endl;
        continue;
      }

      int nx, ny, nz;
      utl::getAttribute(child,"nx",nx);
      utl::getAttribute(child,"ny",ny);
      utl::getAttribute(child,"nz",nz);
      if (nx*ny*nz != width*height) {
        std::cerr << "Invalid dimensions specified " << nx << "x" << ny <<"x" << nz
                  << " image " << width << " " << height << std::endl;
        free(image);
        continue;
      }

      if (nrChannels != 1) {
        std::cerr << "Expect a grayscale image" << std::endl;
        free(image);
        continue;
      }

      properties[prop].textureData.resize(nx,ny,nz);
      const unsigned char* data = image;
      for (int i = 1; i <= nx; ++i)
        for (int j = 1; j <= ny; ++j)
          for (int k = 1; k <= nz; ++k)
            properties[prop].textureData(i,j,k) = double(*data++) / 255.0;

      free(image);

      utl::getAttribute(child,"min",properties[prop].min);
      utl::getAttribute(child,"max",properties[prop].max);
    }
  }
}

void TextureProperties::printLog() const
{
  for (const auto& prop : properties)
    IFEM::cout << "\n\t\tProperty with name " << prop.first
               << " (min = " << prop.second.min
               << ", max = " << prop.second.max << ")";
}


bool TextureProperties::getProperty(const std::string& name,
                                    const Vec3& X, double& val) const
{
  auto it = properties.find(name);
  if (it == properties.end())
    return false;

  const Property& prop = it->second;

  const Vec4* X4 = static_cast<const Vec4*>(&X);
  if (!X4)
    return false;

  int i = std::round(X4->u[0]*(prop.textureData.dim(1)-1));
  int j = std::round(X4->u[1]*(prop.textureData.dim(2)-1));
  int k = std::round(X4->u[2]*(prop.textureData.dim(3)-1));

  if (prop.prescaled)
    val = prop.textureData(i+1,j+1,k+1);
  else
    val = prop.min + (prop.max-prop.min) * prop.textureData(i+1,j+1,k+1);

  return true;
}


bool TextureProperties::hasProperty(const std::string& name) const
{
  return properties.find(name) != properties.end();
}
