// $Id$

#pragma once

#include "DataExporter.h"

class TiXmlDocument;
class TiXmlNode;
class SIMparameters;


class XMLWriter : public DataWriter
{
public:
  struct Entry {
    std::string name;
    std::string description;
    std::string basis;
    int patches;
    int components;
    double timestep;
    int order;
    int interval;
    std::string type;
  };

  XMLWriter(const std::string& name);
  virtual ~XMLWriter() {}

  virtual int getLastTimeLevel();

  //! \brief Calculate the real time level, taking order and ndump into account
  int realTimeLevel(int filelevel) const;

  //! \brief Calculate the real time level, taking order and ndump into account
  int realTimeLevel(int filelevel, int order, int interval) const;

  void readInfo();
  const std::vector<Entry>& getEntries() const { return m_entry; }

  virtual void openFile(int level);
  virtual void closeFile(int level, bool force=false);

  virtual void writeVector(int level, const DataEntry& entry);
  virtual bool readVector(int level, const DataEntry& entry);
  virtual void writeSIM(int level, const DataEntry& entry,
                        bool geometryUpdated);
  virtual bool readSIM(int level, const DataEntry& entry);
  virtual bool writeTimeInfo(int level, int order, int interval, 
                             SIMparameters& tp);

protected:
  void addField(const std::string& name, const std::string& description,
                const std::string& geometry, int components, int patches,
                const std::string& type="field");

  std::vector<Entry> m_entry;

  TiXmlDocument* m_doc;
  TiXmlNode* m_node;
  double m_dt;
  int m_order;
  int m_interval;
};
