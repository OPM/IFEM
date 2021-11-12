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
#include "Functions.h"
#include "HDF5Reader.h"
#include "ProcessAdm.h"
#include "Utilities.h"
#include "Vec3.h"

#include "tinyxml.h"
#include "StbImage.h"


void TextureProperties::parse(const TiXmlElement* elem)
{
  const TiXmlElement* child = elem->FirstChildElement("property");
  for (; child; child = child->NextSiblingElement()) {
    std::string prop;
    utl::getAttribute(child,"name",prop);
    if (prop.empty()) {
      std::cerr << "No name for property, skipping.." << std::endl;
      continue;
    }

    std::string textureFile, function;
    utl::getAttribute(child, "file", textureFile);
    utl::getAttribute(child, "function", function);
    int comp = 1;
    utl::getAttribute(child,"comp",comp);

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
      if (!function.empty()) {
        properties[prop].func_definition = function;
        FunctionBase* func;
        if (comp == 1)
          func = utl::parseRealFunc(function.c_str());
        else
          func = utl::parseVecFunc(function.c_str());

        properties[prop].function.reset(func);
      }

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
  for (const auto& prop : properties) {
    IFEM::cout << "\n\t\tProperty with name " << prop.first
               << " (min = " << prop.second.min
               << ", max = " << prop.second.max << ")";
    if (!prop.second.func_definition.empty())
      IFEM::cout << "\n\t\t\tfunction = " << prop.second.func_definition;
  }
}


bool TextureProperties::getProperty(const std::string& name,
                                    const Vec3& X, double& val) const
{
  auto it = properties.find(name);
  if (it == properties.end())
    return false;

  const Property& prop = it->second;
  val = this->getValue(prop, X);

  if (prop.function) {
    Vec3 f;
    f.x = val;
    val = prop.function->getValue(f).front();
  }

  return true;
}

bool TextureProperties::getProperty(const std::string& name,
                                    const Vec3& X, Vec3& val) const
{
  auto it = properties.find(name);
  if (it == properties.end())
    return false;

  const Property& prop = it->second;
  double value = this->getValue(prop, X);

  if (prop.function) {
    Vec3 f;
    f.x = value;
    val = prop.function->getValue(f);
  } else
    val = value;

  return true;
}


bool TextureProperties::hasProperty(const std::string& name) const
{
  return properties.find(name) != properties.end();
}


double TextureProperties::getValue(const Property& prop, const Vec3& X) const
{
  const Vec4* X4 = static_cast<const Vec4*>(&X);
  if (!X4)
    return false;

  int i = std::round(X4->u[0]*(prop.textureData.dim(1)-1));
  int j = std::round(X4->u[1]*(prop.textureData.dim(2)-1));
  int k = std::round(X4->u[2]*(prop.textureData.dim(3)-1));

  if (prop.prescaled)
    return prop.textureData(i+1,j+1,k+1);
  else
    return prop.min + (prop.max-prop.min) * prop.textureData(i+1,j+1,k+1);

}
