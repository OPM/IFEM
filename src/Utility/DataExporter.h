// $Id$

#pragma once

#include <map>
#include <string>
#include <vector>


class DataWriter;

class DataExporter
{
 public:
  enum FieldType {
    VECTOR,
    SIM
  };

  struct FileEntry {
    std::string description;
    FieldType field;
    int size;
    void* data;
    void* data2;
  };

  //! \brief Default constructor.
  //! \param[in] dynWriters If \e true, delete the writers on destruction.
  DataExporter(bool dynWriters = false) : m_delete(dynWriters), m_level(-1) {}
  //! \brief The desctructor deletes the writers if \a dynWriters was \e true.
  ~DataExporter();

  //! \brief Registers an entry for storage.
  //! param[in] name Name of entry
  //! param[in] description Description of entry
  //! param[in] field Type of entry
  //! param[in] size Number of entries in an array,
  //! the time level to use for SIM
  bool registerField(const std::string& name,
		     const std::string& description,
		     FieldType field, size_t size=0);

  bool registerWriter(DataWriter* writer);

  bool setFieldValue(const std::string& name, void* data, void* data2=NULL);

  bool dumpTimeLevel();

  //! \brief Loads last time level with first registered writer by default.
  bool loadTimeLevel(int level=-1, DataWriter* input=NULL);
  int getTimeLevel();

protected:
  int getWritersTimeLevel() const;

  std::map<std::string,FileEntry> m_entry;
  std::vector<DataWriter*> m_writers;
  bool m_delete;
  int m_level;
};


typedef std::pair<std::string,DataExporter::FileEntry> DataEntry;

class DataWriter
{
protected:
  DataWriter(const std::string& name);

public:
  virtual ~DataWriter() {}

  virtual int getLastTimeLevel() = 0;

  virtual void openFile(int level) = 0;
  virtual void closeFile(int level) = 0;

  virtual void writeVector(int level, const DataEntry& entry) = 0;
  virtual bool readVector(int level, const DataEntry& entry) = 0;
  virtual void writeSIM(int level, const DataEntry& entry) = 0;
  virtual bool readSIM(int level, const DataEntry& entry) = 0;

protected:
  std::string m_name; //!< File name

  int m_size; //!< Number of MPI nodes (processors)
  int m_rank; //!< MPI rank (processor ID)
};
