// $Id$
//==============================================================================
//!
//! \file HDF5Base.h
//!
//! \date Nov 1 2018
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Base class for interacting with HDF5 files.
//!
//==============================================================================

#ifndef _HDF5_BASE_H
#define _HDF5_BASE_H

#include <string>
#ifdef HAS_HDF5
#include <hdf5.h>
#endif

class ProcessAdm;


/*!
  \brief Base class for interacting with HDF5 files.
*/

class HDF5Base
{
public:
  //! \brief The constructor opens a named HDF5-file.
  //! \param[in] name The name (without extension) of the data file
  //! \param[in] adm The process administrator
  HDF5Base(const std::string& name, const ProcessAdm& adm);
  //! \brief The destructor closes the file.
  virtual ~HDF5Base() { this->closeFile(); }

protected:
  //! \brief Opens the HDF5 file.
  //! \param[in] flag Mode to open file using
  //! \param[in] latest If true, open using latest version
  bool openFile(unsigned int flag, bool latest = false);
  //! \brief Closes the HDF5 file.
  void closeFile();

#ifdef HAS_HDF5
  //! \brief Internal helper function checking if a group exists in the file.
  //! \param[in] parent The HDF5 group of the parent
  //! \param[in] path Path of dataset
  //! \return \e true if group exists, otherwise \e false
  static bool checkGroupExistence(hid_t parent, const char* path);

  hid_t        m_file;      //!< The HDF5 handle for our file
#endif
  std::string  m_hdf5_name; //!< The file name of the HDF5 file
#ifdef HAVE_MPI
  const ProcessAdm& m_adm;  //!< Pointer to process administrator in use
#endif
};

#endif
