// $Id$
//==============================================================================
//!
//! \file DataExporter.h
//!
//! \date May 7 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Administer and write data using DataWriters.
//!
//==============================================================================

#pragma once

#include <map>
#include <string>
#include <vector>

class DataWriter;
class TimeStep;


/*!
  \brief Administer and write data using DataWriters.

  \details This class holds a list of data writers,
  and the SIM classes or vectors to write.
*/

class DataExporter
{
 public:
   //! \brief Supported field types
  enum FieldType {
    VECTOR,
    SIM
  };

  //! \brief An enum used to describe the results to write
  enum Results {
    PRIMARY=0,
    SECONDARY=1,
    NORMS=2
  };

  //! \brief A structure holding information about registered fields
  struct FileEntry {
    std::string description; //!< The description of the field
    FieldType   field;       //!< The type of the field
    int         results;     //!< Which results to store
    const void* data;        //!< Pointer to the primary data (e.g. a SIM class)
    const void* data2;       //!< Pointer to the secondary data (e.g. a vector)
  };

  //! \brief Default constructor.
  //! \param[in] dynWriters If \e true, delete the writers on destruction.
  //! \param[in] ndump Interval between dumps
  //! \param[in] order The temporal order of simulations
  //! (always dumps order solutions in a row)
  DataExporter(bool dynWriters = false, int ndump=1, int order=1) :
    m_delete(dynWriters), m_level(-1), m_ndump(ndump), m_order(order) {}

  //! \brief The destructor deletes the writers if \a dynWriters was \e true.
  ~DataExporter();

  //! \brief Registers an entry for storage.
  //! \param[in] name Name of entry
  //! \param[in] description Description of entry
  //! \param[in] field Type of entry
  //! \param[in] results Which results to store
  bool registerField(const std::string& name,
                     const std::string& description,
                     FieldType field, int results = PRIMARY);

  //! \brief Register a data writer
  //! \param[in] writer A pointer to the datawriter we want registered
  bool registerWriter(DataWriter* writer);

  //! \brief Set the data values for a registered field
  //! \param[in] name Name the field is registered with
  //! \param[in] data The value to set the field to
  //! \param[in] data2 (optional) The secondary data of the field
  bool setFieldValue(const std::string& name,
                     const void* data,
                     const void* data2 = NULL);

  //! \brief This dumps all registered fields using all registered writers
  //! \param[in] tp Current time stepping info
  //! \param[in] geometryUpdated Whether or not geometries are updated
  bool dumpTimeLevel(const TimeStep* tp=NULL, bool geometryUpdated=false);

  //! \brief Loads last time level with first registered writer by default.
  //! \param[in] level Time level to load, defaults to last time level
  //! \param[in] info DataWriter to read the info from (e.g. the XML writer)
  //! \param[in] input DataWriter to read the data from (e.g. the HDF5 writer)
  bool loadTimeLevel(int level=-1, DataWriter* info=NULL, DataWriter* input=NULL);

  //! \brief Return the current time level of the exporter
  int getTimeLevel();

  //! \brief Set the prefixes used for norm output
  //! \param[in] prefix The prefixes
  void setNormPrefixes(const char** prefix);
protected:
  //! \brief Internal helper function
  int getWritersTimeLevel() const;

  //! \brief A map of field names -> field info structures
  std::map<std::string,FileEntry> m_entry;
  //! \brief A vector of registered data writers
  std::vector<DataWriter*> m_writers;
  //! \brief If true, we are in charge of freeing up datawriters
  bool m_delete;
  //! \brief Current time level
  int m_level;
  //! \brief A stride for dumping. We dump at every m_ndump'th time level
  int m_ndump;
  //! \brief The temporal order used. We need this to facilitate restart of > first order simulations.
  int m_order;
};

//! \brief Convenience type
typedef std::pair<std::string,DataExporter::FileEntry> DataEntry;


/*!
 \brief Stores and reads data from a file

 \details A DataWriter is a backend for the DataExporter,
 they abstract different file formats.
*/

class DataWriter
{
protected:
  //! \brief Protected constructor as this is a purely virtual class
  DataWriter(const std::string& name);

public:
  //! \brief Empty destructor
  virtual ~DataWriter() {}

  //! \brief Return the last time level stored in file
  virtual int getLastTimeLevel() = 0;

  //! \brief Open the file at a given time level
  //! \param[in] level The requested time level
  virtual void openFile(int level) = 0;

  //! \brief Close the file
  //! \param[in] level Level we just wrote to the file
  //! \param[in] force If true, we always close the actual file,
  //                   else it's up to the individual writers
  virtual void closeFile(int level, bool force=false) = 0;

  //! \brief Write a vector to file
  //! \param[in] level The time level to write the vector at
  //! \param[in] entry The DataEntry describing the vector
  virtual void writeVector(int level, const DataEntry& entry) = 0;

  //! \brief Read a vector from file
  //! \param[in] level The time level to read the vector at
  //! \param[in] entry The DataEntry describing the vector
  virtual bool readVector(int level, const DataEntry& entry) = 0;

  //! \brief Write data from a SIM to file
  //! \param[in] level The time level to write the data at
  //! \param[in] entry The DataEntry describing the vector
  //! \param[in] geometryUpdated Whether or not geometries should be written
  virtual void writeSIM(int level, const DataEntry& entry,
                        bool geometryUpdated) = 0;

  //! \brief Read data from a file into SIM
  //! \param[in] level The time level to read the data at
  //! \param[in] entry The DataEntry describing the SIM
  virtual bool readSIM(int level, const DataEntry& entry) = 0;

  //! \brief Write time stepping info to file (currently a dummy)
  //! \param[in] level The time level to write the info at
  //! \param[in] order The temporal order
  //! \param[in] interval The number of time steps between each data dump
  //! \param[in] tp The current time stepping info
  virtual bool writeTimeInfo(int level, int order, int interval,
                             const TimeStep& tp) = 0;

  //! \brief Set the prefixes used for norm output
  //! \param[in] prefix The prefixes
  void setNormPrefixes(const char** prefix)
  {
    m_prefix = prefix;
  }

protected:
  std::string m_name; //!< File name

  const char** m_prefix; //!< The norm prefixes

  int m_size; //!< Number of MPI nodes (processors)
  int m_rank; //!< MPI rank (processor ID)
};
