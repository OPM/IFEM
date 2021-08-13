// $Id$
//==============================================================================
//!
//! \file HDF5Reader.C
//!
//! \date Nov 1 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Read data from HDF5 file.
//!
//==============================================================================

#include "HDF5Reader.h"

#include <array>
#include <iostream>
#include <sstream>

#ifdef HAS_HDF5
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#endif


HDF5Reader::HDF5Reader (const std::string& name, const ProcessAdm& adm)
  : HDF5Base(((name.find(".hdf5") == std::string::npos &&
               name.find(".h5") == std::string::npos) ? name+".hdf5": name), adm)
{
}


/*!
  \brief Helper for reading data from HDF5 file
*/

#ifdef HAS_HDF5
template<class T>
static bool hdf5Read(hid_t file, hid_t type,
                     const std::string& path, std::vector<T>& vec)
{
  hid_t set = H5Dopen2(file,path.c_str(),H5P_DEFAULT);
  hsize_t size = H5Dget_storage_size(set) / sizeof(T);
  if (size != 0) {
    vec.resize(size);
    H5Dread(set,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,&vec[0]);
    H5Dclose(set);
    return true;
  }
  return false;
}
#endif


bool HDF5Reader::readVector (const std::string& path, std::vector<int>& vec)
{
#ifdef HAS_HDF5
  if (!this->openFile(H5F_ACC_RDONLY))
    return false;

  return hdf5Read(m_file, H5T_NATIVE_INT, path, vec);
#else
  std::cerr << "HDF5Reader: compiled without HDF5 support, no data read" << std::endl;
  return false;
#endif
}


bool HDF5Reader::readVector (const std::string& path, std::vector<double>& vec)
{
#ifdef HAS_HDF5
  if (!this->openFile(H5F_ACC_RDONLY))
    return false;

  return hdf5Read(m_file, H5T_NATIVE_DOUBLE, path, vec);
#else
  std::cerr << "HDF5Reader: compiled without HDF5 support, no data read" << std::endl;
  return false;
#endif
}


bool HDF5Reader::readDouble (const std::string& path, double& out)
{
#ifdef HAS_HDF5
  if (!this->openFile(H5F_ACC_RDONLY))
    return false;

  std::vector<double> vec;
  if (!hdf5Read(m_file, H5T_NATIVE_DOUBLE, path, vec))
    return false;

  out = vec.front();
  return true;
#else
  std::cerr << "HDF5Reader: compiled without HDF5 support, no data read" << std::endl;
  return false;
#endif
}


bool HDF5Reader::readString (const std::string& name, std::string& out)
{
#ifdef HAS_HDF5
  if (!this->openFile(H5F_ACC_RDONLY))
    return false;

  hid_t set = H5Dopen2(m_file,name.c_str(),H5P_DEFAULT);
  hid_t space = H5Dget_space(set);
  hsize_t siz = H5Sget_simple_extent_npoints(space);
  out.resize(siz);
  H5Dread(set,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,H5P_DEFAULT,&out[0]);
  H5Dclose(set);
  return true;
#else
  std::cerr << "HDF5Reader: compiled without HDF5 support, no data read" << std::endl;
  return false;
#endif
}


bool HDF5Reader::read3DArray (const std::string& name, Matrix3D& out)
{
#ifdef HAS_HDF5
  if (!this->openFile(H5F_ACC_RDONLY))
    return false;

  hid_t set = H5Dopen2(m_file,name.c_str(),H5P_DEFAULT);
  hid_t space = H5Dget_space(set);
  hsize_t dims = H5Sget_simple_extent_ndims(space);

  if (dims != 3) {
    std::cerr << "HDF5Reader: Wrong dimension for field " << name << std::endl;
    return false;
  }

  std::array<hsize_t, 3> ndim;

  H5Sget_simple_extent_dims(space, ndim.data(), nullptr);

  out.resize(ndim[0],ndim[1],ndim[2]);

  std::vector<double> test(ndim[0]*ndim[1]*ndim[2]);
  H5Dread(set,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,test.data());
  H5Dclose(set);

  auto it = test.begin();
  for (hsize_t i = 0; i < ndim[0]; ++i)
    for (hsize_t j = 0; j < ndim[1]; ++j)
      for (hsize_t k = 0; k < ndim[2]; ++k, ++it)
        out(i+1, j+1, k+1) = *it;

  return true;
#else
  std::cerr << "HDF5Reader: compiled without HDF5 support, no data read" << std::endl;
  return false;
#endif
}


int HDF5Reader::getFieldSize (const std::string& fieldPath)
{
#ifdef HAS_HDF5
  if (!this->openFile(H5F_ACC_RDONLY))
    return 0;

  int patches = 0;
  while (true) {
    ++patches;
    std::stringstream str;
    str << fieldPath << "/" << patches;
    if (!checkGroupExistence(m_file, str.str().c_str()))
      break;
  }

  return patches-1;
#else
  std::cerr << "HDF5Reader: compiled without HDF5 support, no data read" << std::endl;
  return 0;
#endif
}
