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
namespace utl {
class LogStream;
}


/*!
  \brief Administer and write data using DataWriters.

  \details This class holds a list of data writers,
  and the SIM classes or vectors to write.
*/

class DataExporter : public ControlCallback
{
public:
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
    PRIMARY      = 1, //!< Storage of primary solutions
    DISPLACEMENT = 2, //!< Storage of vector fields as displacements
    SECONDARY    = 4, //!< Storage of secondary field
    NORMS        = 8, //!< Storage of norms
    EIGENMODES   = 16, //!< Storage of eigenmodes
    ONCE         = 32, //!< Only write field once
    GRID         = 128, //!< Always store an updated grid
    REDUNDANT    = 256, //!< Field is redundantly calculated on all processes
    L2G_NODE     = 512 //!< Store local-to-global node mapping
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
    std::vector<const void*> data2; //!< Pointers to the secondary data (e.g. a vector)
    std::string prefix;      //!< Field name prefix
    bool enabled;            //!< Whether or not field is enabled
    int  ncmps;              //!< Number of components. Use to override SIM info
  };

  //! \brief Default constructor.
  //! \param[in] dynWriters If \e true, delete the writers on destruction
  //! \param[in] ndump Interval between dumps
  //! \param[in] level0 Time level to start dumping from (in case of restart)
  DataExporter(bool dynWriters = false, int ndump = 1, int level0 = 0) :
    m_delete(dynWriters), m_ndump(ndump), m_level(level0-1), m_last_step(-1) {}

  //! \brief The destructor deletes the writers if \a dynWriters was \e true.
  virtual ~DataExporter();

  //! \brief Adds the data \a writer to the list of registered writers.
  void registerWriter(DataWriter* writer) { m_writers.push_back(writer); }

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

  //! \brief Sets the data values for a registered field.
  //! \param[in] name Name the field is registered with
  //! \param[in] data The value to set the field to
  //! \param[in] data2 (optional) The secondary data of the field
  //! \param[in] data3 (optional) The third data of the field
  //! \param[in] data4 (optional) The fourth data of the field
  bool setFieldValue(const std::string& name,
                     const void* data,
                     const void* data2=nullptr,
                     const void* data3=nullptr,
                     const void* data4=nullptr);

  //! \brief Dumps all registered fields using the registered writers.
  //! \param[in] tp Current time stepping info
  //! \param[in] geoUpd If \e true, write updated geometry as well
  //! \param[in] doLog If \e true, write a log message when dumping data
  bool dumpTimeLevel(const TimeStep* tp = nullptr,
                     bool geoUpd = false, bool doLog = false);

  //! \brief Callback on receiving a XML control block from external controller.
  virtual void OnControl(const tinyxml2::XMLElement* context);
  //! \brief Returns context name for callback for external controller.
  virtual std::string GetContext() const { return "datawriter"; }

  //! \brief Returns name from (first) data writer.
  std::string getName() const;

  //! \brief Returns current time level of the exporter.
  int getTimeLevel() const { return m_level; }
  //! \brief Returns the data dump stride.
  int getStride() const { return m_ndump; }

protected:
  //! A map of field names -> field info structures
  std::map<std::string,FileEntry> m_entry;
  //! A vector of registered data writers
  std::vector<DataWriter*>        m_writers;

  bool m_delete;    //!< If true, we are in charge of freeing up datawriters
  int  m_ndump;     //!< Time level stride for dumping
  int  m_level;     //!< Current time level
  int  m_last_step; //!< Last time step we dumped for
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
  virtual void closeFile(int level) = 0;

  //! \brief Writes a vector to file.
  //! \param[in] level The time level to write the vector at
  //! \param[in] entry The DataEntry describing the vector
  virtual void writeVector(int level, const DataEntry& entry) = 0;

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

  //! \brief Write a log to output file.
  //! \param data Text to write
  //! \param name Name of log
  virtual bool writeLog(const std::string& data, const std::string& name) = 0;

  //! \brief Returns the name of the file
  const std::string& getName() const { return m_name; }

protected:
  std::string m_name; //!< File name

  int m_size; //!< Number of MPI nodes (processors)
  int m_rank; //!< MPI rank (processor ID)
};

#endif
