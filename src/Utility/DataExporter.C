// $Id$

#include "DataExporter.h"
#include <iostream>
#include <algorithm>
#ifdef PARALLEL_PETSC
#include <mpi.h>
#endif
#include "SIMparameters.h"


DataWriter::DataWriter (const std::string& name) : m_name(name)
{
#ifdef PARALLEL_PETSC
  MPI_Comm_size(MPI_COMM_WORLD,&m_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&m_rank);
#else
  m_size = 1;
  m_rank = 0;
#endif
}


DataExporter::~DataExporter ()
{
  if (m_delete)
    for (size_t i = 0; i < m_writers.size(); i++)
      delete m_writers[i];
}


bool DataExporter::registerField(const std::string& name,
                                 const std::string& description,
                                 FieldType field, int results)
{
  if (m_entry.find(name) != m_entry.end())
    return false;

  FileEntry entry;
  entry.description = description;
  entry.field = field;
  entry.results = results;
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


bool DataExporter::dumpTimeLevel(SIMparameters* tp)
{
  if (tp && tp->step % m_ndump 
         && tp->step % m_ndump > m_order-1)
    return true;

  if (m_level == -1)
    m_level = this->getWritersTimeLevel()+1;

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
	  std::cerr <<"DataExporter: Invalid field type registered, skipping"
		    << std::endl;
          break;
      }
    }
    if (tp)
      (*it2)->writeTimeInfo(m_level,m_order,m_ndump,*tp);
    (*it2)->closeFile(m_level);
  }
  m_level++;

  return true;
}


bool DataExporter::loadTimeLevel (int level, DataWriter* info,
                                  DataWriter* input)
{
  if (!input)
    if (m_writers.empty())
      return false;
    else
      input = m_writers.front();

  if (!info)
    if (m_writers.empty())
      return false;
    else
      info = m_writers.front();

  int level2=level;
  if (level == -1)
    if ((m_level = info->getLastTimeLevel()) < 0)
      return false;
    else
      level2 = m_level;

  bool ok = true;
  input->openFile(level2);
  std::map<std::string,FileEntry>::iterator it;
  for (it = m_entry.begin(); it != m_entry.end() && ok; ++it)
    if (!it->second.data)
      ok = false;
    else switch (it->second.field)
    {
      case VECTOR:
        ok = input->readVector(level2,*it);
        break;
      case SIM:
        ok = input->readSIM(level2,*it);
        break;
      default:
        ok = false;
	std::cerr <<" *** DataExporter: Invalid field type registered "
		  << it->second.field << std::endl;
        break;
    }
  input->closeFile(level2);
  // if we load the last time level, we want to advance
  // if we load a specified time level, we do not want to advance
  if (level == -1)
    m_level++;

  return ok;
}


int DataExporter::getTimeLevel ()
{
  if (m_level == -1)
    m_level = this->getWritersTimeLevel();

  return m_level;
}


int DataExporter::getWritersTimeLevel () const
{
  std::vector<int> levels;
  std::vector<DataWriter*>::const_iterator it2;
  for (it2 = m_writers.begin(); it2 != m_writers.end(); ++it2)
    levels.push_back((*it2)->getLastTimeLevel());
  return *min_element(levels.begin(),levels.end());
}
