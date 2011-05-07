// $Id$

#pragma once

#include "DataExporter.h"
#include "tinyxml.h"


class XMLWriter : public DataWriter {
  public:
    XMLWriter(const std::string& name);
    virtual ~XMLWriter();

    int getLastTimeLevel();

    void openFile(int level);
    void closeFile(int level);

    void writeVector(int level, const DataEntry& entry);
    void readVector(int level, const DataEntry& entry);
    void writeSIM(int level, const DataEntry& entry);
    void readSIM(int level, const DataEntry& entry);

  protected:
    void addField(const std::string& name, const std::string& description,
                  const std::string& type, int patches);
    std::string m_xml;

    int m_rank; // MPI rank
    int m_size; // number of MPI nodes

    TiXmlDocument* m_doc;
    TiXmlNode* m_node;
};
