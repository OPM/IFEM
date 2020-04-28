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
  typedef std::map<std::string,std::string> SerializeData; //!< Convenience type

  //! \brief The constructor opens a named HDF5-file.
  //! \param[in] name The name (without extension) of the data file
  //! \param[in] adm The process administrator
  //! \param[in] stride Restart data output stride
  HDF5Restart(const std::string& name, const ProcessAdm& adm, int stride);

  //! \brief Returns whether or not restart data should be output.
  //! \param[in] tp Time stepping information
  bool dumpStep(const TimeStep& tp);

  //! \brief Writes restart data to file.
  //! \param[in] tp Time stepping information
  //! \param[in] data Data to write
  bool writeData(const TimeStep& tp, const SerializeData& data);

  //! \brief Reads restart data from file.
  //! \param[out] data The map to store data in
  //! \param[in] level Level to read (-1 to read last level in file)
  //! \returns Negative value on error, else restart level loaded
  int readData(SerializeData& data, int level = -1);

private:
  int m_stride; //!< Stride between outputs
};

#endif
