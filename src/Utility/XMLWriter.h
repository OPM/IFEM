// $Id$

#pragma once

#include "DataExporter.h"

class TiXmlDocument;
class TiXmlNode;


class XMLWriter : public DataWriter
{
public:
  struct Entry {
    std::string name;
    std::string description;
    std::string basis;
    int patches;
    int components;
    std::string type;
  };

  XMLWriter(const std::string& name);
  virtual ~XMLWriter() {}

  virtual int getLastTimeLevel();

  void readInfo();
  const std::vector<Entry>& getEntries() const { return m_entry; }

  virtual void openFile(int level);
  virtual void closeFile(int level);

  virtual void writeVector(int level, const DataEntry& entry);
  virtual bool readVector(int level, const DataEntry& entry);
  virtual void writeSIM(int level, const DataEntry& entry);
  virtual bool readSIM(int level, const DataEntry& entry);

protected:
  void addField(const std::string& name, const std::string& description,
                const std::string& geometry, int components, int patches,
                const std::string& type="field");

  std::vector<Entry> m_entry;

  TiXmlDocument* m_doc;
  TiXmlNode* m_node;
};
