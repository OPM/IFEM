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
#include <vector>

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
#ifndef HAS_HDF5
  using hid_t = int; //!< Type alias providing hid_t without the hdf5 library
  using herr_t = int; //!< Type alias providing herr_t without the hdf5 library
#endif
  //! \brief The constructor opens a named HDF5-file.
  //! \param[in] name The name (without extension) of the data file
  //! \param[in] adm The process administrator
  HDF5Base(const std::string& name, const ProcessAdm& adm);

  //! \brief The destructor closes the file.
  virtual ~HDF5Base();

protected:
  //! \brief Open the HDF5 file.
  //!\ param flag Mode to open file using
  bool openFile (unsigned flag);

  //! \brief Close the HDF5 file.
  void closeFile();

  //! \brief Internal helper function checking if a group exists in the file.
  //! \param[in] parent The HDF5 group of the parent
  //! \param[in] path Path of dataset
  //! \return \e true if group exists, otherwise \e false
  bool checkGroupExistence (hid_t parent, const char* path);

  hid_t        m_file; //!< The HDF5 handle for our file
  std::string  m_hdf5_name; //!< The file name of the HDF5 file
#ifdef HAVE_MPI
  const ProcessAdm& m_adm;   //!< Pointer to process adm in use
#endif
};

#endif
