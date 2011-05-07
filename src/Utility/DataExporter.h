// $Id$

#pragma once

#include <map>
#include <string>
#include <vector>


class DataWriter;

class DataExporter {
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

    DataExporter(bool dynamicWriters = false);
    virtual ~DataExporter();

    //! \brief Registers an entry for storage.
    //! param[in] name Name of entry
    //! param[in] description Description of entry
    //! param[in] Type of entry
    //! param[in] size  set to number of entries in an array,
    //                 the time level to use for SIM
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
    int getWritersTimeLevel();

    std::map<std::string,FileEntry> m_entry;
    std::vector<DataWriter*> m_writers;
    bool deleteWriters;
    int m_level;
};


typedef std::pair<std::string,DataExporter::FileEntry> DataEntry;

class DataWriter
{
protected:
  DataWriter() {}

public:
  virtual ~DataWriter() {}

  virtual int getLastTimeLevel() = 0;

  virtual void openFile(int level) = 0;
  virtual void closeFile(int level) = 0;

  virtual void writeVector(int level, const DataEntry& entry) = 0;
  virtual void readVector(int level, const DataEntry& entry) = 0;
  virtual void writeSIM(int level, const DataEntry& entry) = 0;
  virtual void readSIM(int level, const DataEntry& entry) = 0;
};
