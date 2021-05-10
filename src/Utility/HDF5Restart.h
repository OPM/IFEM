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

  \details The HDF5 restart handler writes and reads data using a HDF5 file.
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
  //! \param[in] level0 Time level to start output at (in case of new restart)
  HDF5Restart(const std::string& name, const ProcessAdm& adm,
              int stride = 1, int level0 = 0);

  //! \brief Returns whether or not restart data should be output.
  //! \param[in] tp Time stepping information
  bool dumpStep(const TimeStep& tp);

  //! \brief Writes restart data to file.
  //! \param[in] data Data to write
  bool writeData(const SerializeData& data);

  //! \brief Reads restart data from file.
  //! \param[out] data The map to store data in
  //! \param[in] level Level to read (-1 to read last level in file)
  //! \param[in] basis If \e true, read basis instead of data
  //! \returns Negative value on error, else restart level loaded
  int readData(SerializeData& data, int level = -1, bool basis = false);

  //! \brief Returns current time level.
  int getTimeLevel() const { return m_level; }

private:
  int m_stride; //!< Time step stride between outputs
  int m_level;  //!< Current time level for restart output
  int m_last;   //!< Last time step that was written
};

#endif
