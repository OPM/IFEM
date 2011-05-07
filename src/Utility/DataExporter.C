// $Id$

#include "DataExporter.h"
#include <iostream>
#include <algorithm>


DataExporter::DataExporter(bool dynWrts) : deleteWriters(dynWrts), m_level(-1)
{
}

DataExporter::~DataExporter()
{
  if (deleteWriters)
    for (size_t i = 0; i < m_writers.size(); i++)
      delete m_writers[i];
}


bool DataExporter::registerField(const std::string& name,
                                 const std::string& description,
                                 FieldType field, size_t size)
{
  if (m_entry.find(name) != m_entry.end())
    return false;

  FileEntry entry;
  entry.description = description;
  entry.field = field;
  entry.size = size;
  entry.data = NULL;
  m_entry.insert(make_pair(name,entry));

  return true;
}

bool DataExporter::registerWriter(DataWriter* writer)
{
  m_writers.push_back(writer);
  return true;
}

bool DataExporter::setFieldValue(const std::string& name,
				 void* data, void* data2)
{
  std::map<std::string,FileEntry>::iterator it = m_entry.find(name);
  if (it == m_entry.end())
    return false;

  it->second.data = data;
  it->second.data2 = data2;
  return true;
}

bool DataExporter::dumpTimeLevel()
{
  if (m_level == -1)
    m_level = getWritersTimeLevel()+1;

  std::map<std::string,FileEntry>::iterator it;
  std::vector<DataWriter*>::iterator it2;
  for (it2 = m_writers.begin(); it2 != m_writers.end(); ++it2) {
    (*it2)->openFile(m_level);
    for (it = m_entry.begin(); it != m_entry.end(); ++it) {
      if (!it->second.data)
        return false;
      switch (it->second.field) {
        case VECTOR:
          (*it2)->writeVector(m_level,*it);
          break;
        case SIM:
          (*it2)->writeSIM(m_level,*it);
          break;
        default:
	  std::cout <<"DataExporter: Invalid field type registered, skipping"
		    << std::endl;
          break;
      }
    }
    (*it2)->closeFile(m_level);
  }
  m_level++;

  return true;
}

bool DataExporter::loadTimeLevel(int level, DataWriter* input)
{
  if (!input && m_writers.empty())
    return false;
  if (!input)
    input = m_writers.front();
  if (level == -1)
    level = m_level = input->getLastTimeLevel();
  if (level == -1)
    return false;

  input->openFile(level);
  std::map<std::string,FileEntry>::iterator it;
  for (it = m_entry.begin(); it != m_entry.end(); ++it)
  {
    if (!it->second.data)
      return false;
    switch (it->second.field)
    {
      case VECTOR:
        input->readVector(level,*it);
        break;
      case SIM:
        input->readSIM(level,*it);
        break;
      default:
	std::cout << "DataExporter::loadTimeLevel: Invalid field type "
		  << "registered, skipping" << std::endl;
        break;
    }
  }
  input->closeFile(level);
  m_level++;

  return true;
}

int DataExporter::getTimeLevel()
{
  if (m_level == -1)
    m_level = getWritersTimeLevel();

  return m_level;
}

int DataExporter::getWritersTimeLevel()
{
  std::vector<int> levels;
  std::vector<DataWriter*>::const_iterator it2;
  for (it2 = m_writers.begin(); it2 != m_writers.end(); ++it2)
    levels.push_back((*it2)->getLastTimeLevel());
  return *min_element(levels.begin(),levels.end());
}
