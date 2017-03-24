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

#ifndef _DATA_EXPORTER_H
#define _DATA_EXPORTER_H

#include "ControlFIFO.h"
#include <map>
#include <string>
#include <vector>

class DataWriter;
class ProcessAdm;
class TimeStep;


/*!
  \brief Administer and write data using DataWriters.

  \details This class holds a list of data writers,
  and the SIM classes or vectors to write.
*/

class DataExporter : public ControlCallback
{
public:
  typedef std::map<std::string, std::string> SerializeData; //!< Convenience typedef

   //! \brief Supported field types
  enum FieldType {
    VECTOR,
    INTVECTOR,
    KNOTSPAN,
    SIM,
    NODALFORCES,
    BASIS
  };

  //! \brief An enum used to describe the results to write from a SIM
  enum Results {
    PRIMARY      = 1,   //!< Storage of primary solutions
    DISPLACEMENT = 2,   //!< Storage of vector fields as displacements
    SECONDARY    = 4,   //!< Storage of secondary field
    NORMS        = 8,   //!< Storage of norms
    EIGENMODES   = 16,  //!< Storage of eigenmodes
    ONCE         = 32,  //!< Only write field once
    GRID         = 128, //!< Always store an updated grid
    REDUNDANT    = 256  //!< Field is redundantly calculated on all processes
  };

  //! \brief A structure holding information about registered fields
  struct FileEntry {
    std::string description; //!< The description of the field
    FieldType   field;       //!< The type of the field
    int         results;     //!< \brief Which results to store.
                             //! \details A negative value indicates that we
                             //! want to use the description as name for the
                             //! primary vector, not the name of the Integrand.
    const void* data;        //!< Pointer to the primary data (e.g. a SIM class)
    const void* data2;       //!< Pointer to the secondary data (e.g. a vector)
    std::string prefix;      //!< Field name prefix
    bool enabled;            //!< Whether or not field is enabled
    int  ncmps;              //!< Number of components. Use to override SIM info
  };

  //! \brief Default constructor.
  //! \param[in] dynWriters If \e true, delete the writers on destruction
  //! \param[in] ndump Interval between dumps
  //! \param[in] nrestart Restart stride. 0 to disable
  DataExporter(bool dynWriters = false, int ndump=1, int nrestart=0) :
    m_delete(dynWriters), m_level(-1), m_ndump(ndump),
    m_last_step(-1), m_nrestart(nrestart), m_infoReader(0), m_dataReader(0) {}

  //! \brief The destructor deletes the writers if \a dynWriters was \e true.
  virtual ~DataExporter();

  //! \brief Registers an entry for storage.
  //! \param[in] name Name of entry
  //! \param[in] description Description of entry
  //! \param[in] field Type of entry
  //! \param[in] results Which results to store
  //! \param[in] prefix Field name prefix
  //! \param[in] ncmps Number of field components
  bool registerField(const std::string& name,
                     const std::string& description,
                     FieldType field, int results = PRIMARY,
                     const std::string& prefix = "", int ncmps = 0);

  //! \brief Registers a data writer.
  //! \param[in] writer A pointer to the data writer we want registered
  //! \param info If \e true, set as default info reader
  //! \param reader If \e true, set as default data reader
  bool registerWriter(DataWriter* writer, bool info=false, bool reader=false);

  //! \brief Sets the data values for a registered field.
  //! \param[in] name Name the field is registered with
  //! \param[in] data The value to set the field to
  //! \param[in] data2 (optional) The secondary data of the field
  bool setFieldValue(const std::string& name,
                     const void* data,
                     const void* data2=nullptr);

  //! \brief Dumps all registered fields using the registered writers.
  //! \param[in] tp Current time stepping info
  //! \param[in] geometryUpdated Whether or not geometries are updated
  //! \param[in] serializeData Serialized data from simulators for restart files
  bool dumpTimeLevel(const TimeStep* tp=nullptr, bool geometryUpdated=false,
                     SerializeData* serializeData = nullptr);

  //! \brief Loads last time level with first registered writer by default.
  //! \param[in] level Time level to load, defaults to last time level
  //! \param[in] info DataWriter to read the info from (e.g. the XML writer)
  //! \param[in] input DataWriter to read the data from (e.g. the HDF5 writer)
  bool loadTimeLevel(int level=-1,
                     DataWriter* info=nullptr, DataWriter* input=nullptr);

  //! \brief Returns the current time level of the exporter.
  int getTimeLevel();

  //! \brief Calculates the real time level taking ndump into account.
  int realTimeLevel(int filelevel) const;

  //! \brief Sets the prefices used for norm output.
  void setNormPrefixes(const char** prefix);

  //! \brief Callback on receiving a XML control block from external controller.
  virtual void OnControl(const TiXmlElement* context);
  //! \brief Returns context name for callback for external controller.
  virtual std::string GetContext() const { return "datawriter"; }

  //! \brief Return name from data writer
  std::string getName() const;

  //! \brief Returns visualization data stride
  int getStride() const { return m_ndump; }

  //! \brief Returns whether current step should be saved for restart or not.
  bool dumpForRestart(const TimeStep* tp) const;

protected:
  //! \brief Internal helper function.
  int getWritersTimeLevel() const;

  //! A map of field names -> field info structures
  std::map<std::string,FileEntry> m_entry;
  //! A vector of registered data writers
  std::vector<DataWriter*>        m_writers;

  bool m_delete;    //!< If true, we are in charge of freeing up datawriters
  int  m_level;     //!< Current time level
  int  m_ndump;     //!< Time level stride for dumping
  int  m_last_step; //!< Last time step we dumped for
  int  m_nrestart;  //!< Stride for restart data dumping

  DataWriter* m_infoReader; //!< DataWriter to read data information from
  DataWriter* m_dataReader; //!< DataWriter to read numerical data from
};

//! \brief Convenience type
typedef std::pair<std::string,DataExporter::FileEntry> DataEntry;


/*!
  \brief Stores and reads data from a file.

  \details A DataWriter is a backend for the DataExporter,
  they abstract different file formats.
*/

class DataWriter
{
protected:
  //! \brief Protected constructor as this is a purely virtual class.
  DataWriter(const std::string& name, const ProcessAdm& adm,
             const char* defaultExt = nullptr);

public:
  //! \brief Empty destructor.
  virtual ~DataWriter() {}

  //! \brief Returns the last time level stored in file.
  virtual int getLastTimeLevel() = 0;

  //! \brief Opens the file at a given time level.
  //! \param[in] level The requested time level
  virtual void openFile(int level) = 0;

  //! \brief Closes the file.
  //! \param[in] level Level we just wrote to the file
  //! \param[in] force If true, we always close the actual file,
  //! otherwise it's up to the individual writers
  virtual void closeFile(int level, bool force = false) = 0;

  //! \brief Writes a vector to file.
  //! \param[in] level The time level to write the vector at
  //! \param[in] entry The DataEntry describing the vector
  virtual void writeVector(int level, const DataEntry& entry) = 0;

  //! \brief Reads a vector from file.
  //! \param[in] level The time level to read the vector at
  //! \param[in] entry The DataEntry describing the vector
  virtual bool readVector(int level, const DataEntry& entry) = 0;

  //! \brief Writes data from a SIM object to file.
  //! \param[in] level The time level to write the data at
  //! \param[in] entry The DataEntry describing the vector
  //! \param[in] geometryUpdated Whether or not geometries should be written
  //! \param[in] prefix Field name prefix
  virtual void writeSIM(int level, const DataEntry& entry,
                        bool geometryUpdated, const std::string& prefix) = 0;

  //! \brief Writes nodal forces to file.
  //! \param[in] level The time level to write the data at
  //! \param[in] entry The DataEntry describing the vector
  virtual void writeNodalForces(int level, const DataEntry& entry) = 0;

  //! \brief Writes a knotspan field to file.
  //! \param[in] level The time level to write the data at
  //! \param[in] entry The DataEntry describing the field
  //! \param[in] prefix Prefix for field
  virtual void writeKnotspan(int level, const DataEntry& entry,
                             const std::string& prefix) = 0;

  //! \brief Write a basis to file
  //! \param[in] level The time level to write the basis at
  //! \param[in] entry The DataEntry describing the basis
  //! \param[in] prefix Prefix for basis
  virtual void writeBasis(int level, const DataEntry& entry,
                          const std::string& prefix) = 0;

  //! \brief Writes time stepping info to file.
  //! \param[in] level The time level to write the info at
  //! \param[in] interval The number of time steps between each data dump
  //! \param[in] tp The current time stepping info
  virtual bool writeTimeInfo(int level, int interval,
                             const TimeStep& tp) = 0;

  //! \brief Sets the prefices used for norm output.
  void setNormPrefixes(const char** prefix) { m_prefix = prefix; }

  //! \brief Returns the name of the file
  const std::string& getName() const { return m_name; }

  //! \brief Write restart data.
  //! \param[in] step Level to write restart data at
  //! \param[in] data Data to write
  virtual bool writeRestartData(int level, const DataExporter::SerializeData& data) = 0;

protected:
  std::string  m_name;   //!< File name
  const char** m_prefix; //!< The norm prefixes

  int m_size; //!< Number of MPI nodes (processors)
  int m_rank; //!< MPI rank (processor ID)
};

#endif
