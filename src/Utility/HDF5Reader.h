// $Id$
//==============================================================================
//!
//! \file HDF5Reader.h
//!
//! \date Nov 1 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Read data from HDF5 file.
//!
//==============================================================================

#ifndef _HDF5_READER_H
#define _HDF5_READER_H

#include "HDF5Base.h"
#include "MatVec.h"

#include <string>
#include <vector>


/*!
  \brief Read data from a HDF5 file.
*/

class HDF5Reader : public HDF5Base
{
public:
  //! \brief The constructor opens a named HDF5-file.
  //! \param[in] name The name (without extension) of the data file
  //! \param[in] adm The process administrator
  HDF5Reader(const std::string& name, const ProcessAdm& adm);

  //! \brief Reads an integer vector.
  //! \param[in] path The path to the dataset to read
  //! \param[out] vec The vector to read data into
  bool readVector(const std::string& path, std::vector<int>& vec);

  //! \brief Reads an double precision vector.
  //! \param[in] path The path to the dataset to read
  //! \param[out] vec The vector to read data into
  bool readVector(const std::string& path, std::vector<double>& vec);

  //! \brief Reads a single double value.
  //! \param[in] name The name (path in HDF5 file) to the string
  //! \param[out] out The double to read data into
  bool readDouble(const std::string& name, double& out);

  //! \brief Reads a text string.
  //! \param[in] name The name (path in HDF5 file) to the string
  //! \param[out] out The string to read data into
  bool readString(const std::string& name, std::string& out);

  //! \brief Reads a 3D array.
  //! \param[in] name The name (path in HDF5 file) to the string
  //! \param[out] data The 3D array to read data into
  bool read3DArray(const std::string& name, Matrix3D& data);

  //! \brief Returns number of patches for a field.
  //! \param[in] fieldPath Path to field in hdf5 file
  int getFieldSize(const std::string& fieldPath);
};

#endif
