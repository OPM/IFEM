// $Id$
//==============================================================================
//!
//! \file HDF5Writer.h
//!
//! \date Jul 7 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Output of model and results to HDF5 file.
//!
//==============================================================================

#ifndef _HDF5_WRITER_H
#define _HDF5_WRITER_H

#include "DataExporter.h"
#include "HDF5Base.h"

class SIMbase;


/*!
  \brief Write data to a HDF5 file.

  \details The HDF5 writer writes data to a HDF5 file. It supports parallel I/O.
*/

class HDF5Writer : public DataWriter, public HDF5Base
{
public:
  //! \brief The constructor opens a named HDF5-file.
  //! \param[in] name The name (without extension) of the data file
  //! \param[in] adm The process administrator
  //! \param[in] append Whether to append to or overwrite an existing file
  HDF5Writer(const std::string& name, const ProcessAdm& adm,
             bool append = false);

  //! \brief Empty destructor.
  virtual ~HDF5Writer() {}

  //! \brief Returns the last time level stored in the HDF5 file.
  virtual int getLastTimeLevel();

  //! \brief Opens the file at a given time level.
  //! \param[in] level The requested time level
  virtual void openFile(int level);

  //! \brief Closes the file.
  //! \param[in] level Level we just wrote to the file
  virtual void closeFile(int level);

  //! \brief Writes a vector to file.
  //! \param[in] level The time level to write the vector at
  //! \param[in] entry The DataEntry describing the vector
  virtual void writeVector(int level, const DataEntry& entry);

  //! \brief Writes data from a SIM to file.
  //! \param[in] level The time level to write the data at
  //! \param[in] entry The DataEntry describing the data to write
  //! \param[in] geometryUpdated Whether or not geometries should be written
  //! \param[in] prefix Field name prefix
  //!
  //! \details If prefix is non-empty and we are asked to dump secondary
  //! solutions, we assume they come from different projections
  //! \sa SIMbase::project
  virtual void writeSIM(int level, const DataEntry& entry,
                        bool geometryUpdated, const std::string& prefix);

  //! \brief Writes nodal forces to file.
  //! \param[in] level The time level to write the data at
  //! \param[in] entry The DataEntry describing the data to write
  virtual void writeNodalForces(int level, const DataEntry& entry);

  //! \brief Writes knot span field to file.
  //! \param[in] level The time level to write the data at
  //! \param[in] entry The DataEntry describing the field
  //! \param[in] prefix Prefix for field
  virtual void writeKnotspan(int level, const DataEntry& entry,
                             const std::string& prefix);

  //! \brief Writes a basis to file.
  //! \param[in] level The time level to write the basis at
  //! \param[in] entry The DataEntry describing the basis
  //! \param[in] prefix Prefix for basis
  virtual void writeBasis(int level, const DataEntry& entry,
                          const std::string& prefix);

  //! \brief Writes time stepping info to file.
  //! \param[in] level The time level to write the info at
  //! \param[in] interval The number of time steps between each data dump
  //! \param[in] tp The current time stepping info
  virtual bool writeTimeInfo(int level, int interval, const TimeStep& tp);

  //! \brief Write a log to output file.
  //! \param data Text to write
  //! \param name Name of log
  virtual bool writeLog(const std::string& data, const std::string& name);

  //! \brief Write analytical solutions to file.
  //! \param aSol Map of functions
  //! \param name Name of simulator
  bool writeAnaSol(const std::map<std::string, std::string>& aSols,
                   const std::string& name);

#ifdef HAS_HDF5
protected:
  //! \brief Internal helper function writing a data array to file.
  //! \param[in] group The HDF5 group to write data into
  //! \param[in] name The name of the array
  //! \param[in] patch Patch number of the array
  //! \param[in] len The length of the array
  //! \param[in] data The array to write
  //! \param[in] type The HDF5 type for the data (see H5T)
  void writeArray(hid_t group, const std::string& name,
                  int patch, int len, const void* data, hid_t type);

  //! \brief Internal helper function writing a SIM's basis (geometry) to file.
  //! \param[in] SIM The SIM we want to write basis for
  //! \param[in] name The name of the basis
  //! \param[in] basis 1/2 Write primary or secondary basis from SIM
  //! \param[in] level The time level to write the basis at
  //! \param[in] redundant Whether or not basis is redundant across processes
  //! \param[in] l2g True to write local-to-global node numbers
  void writeBasis(const SIMbase* SIM, const std::string& name,
                  int basis, int level, bool redundant = false,
                  bool l2g = false);

private:
  unsigned int m_flag; //!< The file flags to open HDF5 file with
#endif
};

#endif
