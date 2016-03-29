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

class SIMbase;


/*!
  \brief Write data to a HDF5 file.

  \details The HDF5 writer writes data to a HDF5 file.
  It supports parallel I/O, and can be used to add restart capability
  to applications.
*/

class HDF5Writer : public DataWriter
{
public:
  //! \brief The constructor opens a named HDF5-file.
  //! \param[in] name The name (filename without extension) of data file
  //! \param[in] append Whether to append to or overwrite an existing file
  //! \param[in] keepopen Whether to always keep the HDF5 open
  HDF5Writer(const std::string& name, const ProcessAdm& adm, bool append = false,
             bool keepopen = false);

  //! \brief Empty destructor.
  virtual ~HDF5Writer() {}

  //! \brief Returns the last time level stored in the HDF5 file.
  virtual int getLastTimeLevel();

  //! \brief Opens the file at a given time level.
  //! \param[in] level The requested time level
  virtual void openFile(int level);

  //! \brief Closes the file.
  //! \param[in] level Level we just wrote to the file
  //! \param[in] force If \e true, close even if \a keepopen was \e true
  virtual void closeFile(int level, bool force = false);

  //! \brief Writes a vector to file.
  //! \param[in] level The time level to write the vector at
  //! \param[in] entry The DataEntry describing the vector
  virtual void writeVector(int level, const DataEntry& entry);

  //! \brief Reads a vector from file.
  //! \param[in] level The time level to read the vector at
  //! \param[in] entry The DataEntry describing the vector
  virtual bool readVector(int level, const DataEntry& entry);

  //! \brief Writes data from a SIM to file.
  //! \param[in] level The time level to write the data at
  //! \param[in] entry The DataEntry describing the vector
  //! \param[in] geometryUpdated Whether or not geometries should be written
  //! \param[in] prefix Field name prefix
  //!
  //! \details If prefix is non-empty and we are asked to dump secondary
  //! solutions, we assume they come from different projections
  //! \sa SIMbase::project
  virtual void writeSIM(int level, const DataEntry& entry,
                        bool geometryUpdated, const std::string& prefix);

  //! \brief Reads data from a file into a SIM.
  //! \param[in] level The time level to read the data at
  //! \param[in] entry The DataEntry describing the SIM
  virtual bool readSIM(int level, const DataEntry& entry);

  //! \brief Writes nodal forces to file.
  //! \param[in] level The time level to write the data at
  //! \param[in] entry The DataEntry describing the vector
  virtual void writeNodalForces(int level, const DataEntry& entry);

  //! \brief Writes knot span field to file.
  //! \param[in] level The time level to write the data at
  //! \param[in] entry The DataEntry describing the field
  //! \param[in] prefix Prefix for field
  virtual void writeKnotspan(int level, const DataEntry& entry,
                             const std::string& prefix);

  //! \brief Write a basis to file
  //! \param[in] level The time level to write the basis at
  //! \param[in] entry The DataEntry describing the basis
  //! \param[in] prefix Prefix for basis
  virtual void writeBasis(int level, const DataEntry& entry,
                          const std::string& prefix);

  //! \brief Writes time stepping info to file.
  //! \param[in] level The time level to write the info at
  //! \param[in] order The temporal order
  //! \param[in] interval The number of time steps between each data dump
  //! \param[in] tp The current time stepping info
  virtual bool writeTimeInfo(int level, int order, int interval,
                             const TimeStep& tp);

  //! \brief Reads a text string.
  //! \param[in] name The name (path in HDF5 file) to the string
  //! \param[out] out The string to read data into
  //! \param[in] close If \e false, keep the HDF5-file open after reading
  void readString(const std::string& name, std::string& out, bool close = true);

  //! \brief Reads a double vector.
  //! \param[in] level The time level to read at
  //! \param[in] name The name (path in HDF5 file) to the string
  //! \param[in] patch The patch to read
  //! \param[out] vec The vector to read data into
  bool readVector(int level, const std::string& name,
                  int patch, std::vector<double>& vec);

  //! \brief Reads an integer vector.
  //! \param[in] level The time level to read at
  //! \param[in] name The name (path in HDF5 file) to the string
  //! \param[in] patch The patch to read
  //! \param[out] vec The vector to read data into
  bool readVector(int level, const std::string& name,
                  int patch, std::vector<int>& vec);

  //! \brief Reads a double value.
  //! \param[in] level The time level to read at
  //! \param[in] name The name (path in HDF5 file) to the array
  //! \param[in] group The HDF5 group to read from
  //! \param[out] data The variable to read data into
  bool readDouble(int level, const std::string& group,
                  const std::string& name, double& data);

  //! \brief Checks if updated geometries exist in file at given time level.
  //! \param[in] level The time level to check
  //! \param[in] basisName Check for a particular basis
  bool hasGeometries(int level, const std::string& basisName = "");

protected:
  //! \brief Internal helper function. Writes a data array to HDF5 file.
  //! \param[in] group The HDF5 group to write data into
  //! \param[in] name The name of the array
  //! \param[in] len The length of the array
  //! \param[in] data The array to write
  //! \param[in] type The HDF5 type for the data (see H5T)
  void writeArray(int group, const std::string& name,
                  int len, const void* data, int type);

  //! \brief Internal helper function. Writes a SIM's basis (geometry) to file.
  //! \param[in] SIM The SIM we want to write basis for
  //! \param[in] name The name of the basis
  //! \param[in] basis 1/2 Write primary or secondary basis from SIM
  //! \param[in] level The time level to write the basis at
  //! \param[in] redundant Whether or not basis is redundant across processes
  void writeBasis(SIMbase* SIM, const std::string& name, int basis,
                  int level, bool redundant=false);

  //! \brief Internal helper function. Reads an array into a array of doubles.
  //! \param[in] group The HDF5 group to read data from
  //! \param[in] name The name of the array
  //! \param[in] len The length of the data to read
  //! \param[out] data The array to read data into
  void readArray(int group, const std::string& name, int& len, double*& data);

  //! \brief Internal helper function. Reads an array into a array of integers.
  //! \param[in] group The HDF5 group to read data from
  //! \param[in] name The name of the array
  //! \param[in] len The length of the data to read
  //! \param[out] data The array to read data into
  void readArray(int group, const std::string& name, int& len, int*& data);

  //! \brief Internal helper function. Checks if a group exists in the file.
  //! \param[in] parent The HDF5 group of the parent
  //! \param[in] group The name of the group to check for
  //! \return \e true if group exists, otherwise \e false
  bool checkGroupExistence(int parent, const char* group);

private:
  int          m_file; //!< The HDF5 handle for our file
  unsigned int m_flag; //!< The file flags to open HDF5 file with
  bool     m_keepOpen; //!< If \e true, we always keep the file open
#ifdef HAVE_MPI
  const ProcessAdm& m_adm;   //!< Pointer to process adm in use
#endif
};

#endif
