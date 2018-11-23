// $Id$
//==============================================================================
//!
//! \file HDF5Restart.h
//!
//! \date Nov 1 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Output of restart data to HDF5.
//!
//==============================================================================

#ifndef _HDF5_RESTART_H
#define _HDF5_RESTART_H

#include "HDF5Base.h"

#include <map>
#include <string>


class ProcessAdm;
class TimeStep;


/*!
  \brief Write and read restart data using a HDF5 file.

  \details The HDF5 restart hanlder writes and reads data using a HDF5 file.
  It supports parallel I/O, and can be used to add restart capability
  to applications.
*/

class HDF5Restart : public HDF5Base
{
public:
  typedef std::map<std::string, std::string> SerializeData; //!< Convenience typedef

  //! \brief The constructor opens a named HDF5-file.
  //! \param[in] name The name (without extension) of the data file
  //! \param[in] adm The process administrator
  //! \param[in] stride Restart data output stride
  HDF5Restart(const std::string& name, const ProcessAdm& adm, int stride);

  //! \brief Returns whether or not restart data should be output.
  //! \param tp Time stepping parameter
  bool dumpStep(const TimeStep& tp);

  //! \brief Writes restart data to file.
  //! \param[in] tp Time stepping information
  //! \param[in] data Data to write
  bool writeData(const TimeStep& tp, const SerializeData& data);

  //! \brief Reads restart data from file.
  //! \param[out] data The map to store data in
  //! \param[in] level Level to read (-1 to read last level in file)
  //! \param[in] basis If \e true, read basis instead of data
  //! \returns Negative value on error, else restart level loaded
  int readData(SerializeData& data, int level = -1, bool basis = false);

protected:
  //! \brief Internal helper function reading into an array of chars.
  //! \param[in] group The HDF5 dataset to read data from
  //! \param[in] name The name of the array
  //! \param[out] len The length of the data read
  //! \param[out] data The array to read data into
  void readArray(hid_t group, const std::string& name, int& len, char*& data);

  //! \brief Internal helper function writing a data array to file.
  //! \param[in] group The HDF5 group to write data into
  //! \param[in] name The name of the array
  //! \param[in] len The length of the array
  //! \param[in] data The array to write
  //! \param[in] type The HDF5 type for the data (see H5T)
  void writeArray(hid_t group, const std::string& name,
                  int len, const void* data, hid_t type);

  //! \brief Static helper used for reading data.
  static herr_t read_restart_data(hid_t group_id, const char* member_name, void* data);

private:
  int m_stride;        //!< Stride between outputs
};

#endif
