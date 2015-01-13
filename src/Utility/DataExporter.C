// $Id$
//==============================================================================
//!
//! \file DataExporter.C
//!
//! \date May 7 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Administer and write data using DataWriters.
//!
//==============================================================================

#include "DataExporter.h"
#include "Utilities.h"
#include "ProcessAdm.h"
#include "TimeStep.h"
#include "tinyxml.h"
#include <iostream>
#include <algorithm>


DataWriter::DataWriter (const std::string& name,
                        const ProcessAdm& adm,
                        const char* defaultExt)
{
  if (defaultExt && name.find(defaultExt) == std::string::npos)
    m_name = name + defaultExt;
  else
    m_name = name;

  m_prefix = NULL;

  m_size = adm.getNoProcs();
  m_rank = adm.getProcId();
}


DataExporter::~DataExporter ()
{
  if (m_delete)
    for (size_t i = 0; i < m_writers.size(); i++)
      delete m_writers[i];
}


bool DataExporter::registerField (const std::string& name,
				  const std::string& description,
				  FieldType field, int results,
                                  const std::string& prefix,
                                  int ncmps)
{
  if (m_entry.find(name) != m_entry.end())
    return false;

  FileEntry entry;
  entry.description = description;
  entry.field = field;
  entry.results = results;
  entry.data = entry.data2 = NULL;
  entry.prefix = prefix;
  if (!prefix.empty())
    entry.prefix += ' ';
  entry.enabled = true;
  entry.ncmps = ncmps;
  m_entry.insert(std::make_pair(name,entry));

  return true;
}


bool DataExporter::registerWriter (DataWriter* writer, bool info, bool data)
{
  m_writers.push_back(writer);

  if (info)
    m_infoReader = writer;
  if (data)
    m_dataReader = writer;

  return true;
}


bool DataExporter::setFieldValue (const std::string& name,
				  const void* data, const void* data2)
{
  std::map<std::string,FileEntry>::iterator it = m_entry.find(name);
  if (it == m_entry.end())
    return false;

  it->second.data = data;
  it->second.data2 = data2;
  return true;
}


bool DataExporter::dumpTimeLevel (const TimeStep* tp, bool geometryUpdated)
{
  if (tp && tp->step % m_ndump && tp->step % m_ndump > m_order)
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
        case INTVECTOR:
        case VECTOR:
          (*it2)->writeVector(m_level,*it);
          break;
        case SIM:
          (*it2)->writeSIM(m_level,*it,geometryUpdated,it->second.prefix);
          break;
        case NODALFORCES:
          (*it2)->writeNodalForces(m_level,*it);
          break;
        case KNOTSPAN:
          (*it2)->writeKnotspan(m_level,*it,it->second.prefix);
          break;
        case BASIS:
          (*it2)->writeBasis(m_level,*it,it->second.prefix);
          break;
        default:
          std::cerr <<"  ** DataExporter: Invalid field type registered "
                    << it->second.field <<", skipping"<< std::endl;
          break;
      }
    }
    if (tp)
      (*it2)->writeTimeInfo(m_level,m_order,m_ndump,*tp);

    (*it2)->closeFile(m_level);
  }
  m_level++;

  // disable fields marked as once
  for (it = m_entry.begin(); it != m_entry.end(); ++it)
    if (abs(it->second.results) & ONCE)
      it->second.enabled = false;

  return true;
}


bool DataExporter::loadTimeLevel (int level, DataWriter* info,
                                  DataWriter* input)
{
  if (!input)
    if (m_writers.empty())
      return false;
    else if (m_dataReader)
      input = m_dataReader;
    else
      input = m_writers.front();

  if (!info)
    if (m_infoReader)
      info = m_infoReader;
    else if (m_writers.size() < 2)
      return false;
    else
      info = m_writers[1];

  int level2=level;
  if (level == -1)
    if ((m_level = info->getLastTimeLevel()) < 0)
      return false;
    else
      level2 = m_level;

  bool ok = true;
  input->openFile(level2);
  std::map<std::string,FileEntry>::iterator it;
  for (it = m_entry.begin(); it != m_entry.end() && ok; ++it) {
    if (!it->second.data)
      ok = false;
    else switch (it->second.field)
    {
      case SIM:
        ok = input->readSIM(level2,*it);
        break;
      default:
        break;
    }
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


void DataExporter::setNormPrefixes(const char** prefix)
{
  for (std::vector<DataWriter*>::iterator it  = m_writers.begin();
                                          it != m_writers.end(); ++it)
    (*it)->setNormPrefixes(prefix);
}


int DataExporter::realTimeLevel(int filelevel) const
{
  return realTimeLevel(filelevel,m_order,m_ndump);
}


int DataExporter::realTimeLevel(int filelevel, int order, int interval) const
{
  return filelevel/order*interval;
}


void DataExporter::OnControl(const TiXmlElement* context)
{
  const TiXmlElement* child = context->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (strcasecmp(child->Value(),"enable_field") == 0) {
      std::string name;
      bool enable = true;
      utl::getAttribute(child,"name",name);
      utl::getAttribute(child,"enable",enable);
      std::map<std::string,FileEntry>::iterator it = m_entry.find(name);
      if (it != m_entry.end()) {
        std::cout <<"  * DataWriter: "<< (enable ? "Enabled " : "Disabled ")
                  << name << std::endl;
        it->second.enabled = enable;
      }
    }
    else if (strcasecmp(child->Value(),"set_stride") == 0)
      if (utl::getAttribute(child,"value",m_ndump))
        std::cout <<"  * DataWriter: set stride "<< m_ndump << std::endl;
}


std::string DataExporter::getName() const
{
  if (!m_writers.empty())
    return m_writers.front()->getName().substr(0,m_writers.front()->getName().rfind('.'));

  return "";
}
