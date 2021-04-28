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
#include "IFEM.h"
#include "Utilities.h"
#include "Profiler.h"
#include "ProcessAdm.h"
#include "TimeStep.h"
#include "tinyxml.h"


DataWriter::DataWriter (const std::string& name,
                        const ProcessAdm& adm,
                        const char* defaultExt)
{
  if (defaultExt && name.find(defaultExt) == std::string::npos)
    m_name = name + defaultExt;
  else
    m_name = name;

  m_size = adm.getNoProcs();
  m_rank = adm.getProcId();
}


DataExporter::~DataExporter ()
{
  if (IFEM::memoryLog)
    for (DataWriter* writer : m_writers) {
      writer->openFile(0);
      writer->writeLog(IFEM::memoryLog->str(), "cout");
      writer->closeFile(0);
    }

  if (m_delete)
    for (DataWriter* writer : m_writers)
      delete writer;
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
  entry.prefix = prefix;
  if (!prefix.empty())
    entry.prefix += ' ';
  entry.enabled = true;
  entry.ncmps = ncmps;
  m_entry.insert(std::make_pair(name,entry));

  return true;
}


bool DataExporter::setFieldValue (const std::string& name,
                                  const void* data,
                                  const void* data2,
                                  const void* data3,
                                  const void* data4)
{
  std::map<std::string,FileEntry>::iterator it = m_entry.find(name);
  if (it == m_entry.end())
    return false;

  it->second.data = data;
  it->second.data2.clear();
  it->second.data2.push_back(data2);
  it->second.data2.push_back(data3);
  it->second.data2.push_back(data4);
  return true;
}


bool DataExporter::dumpTimeLevel (const TimeStep* tp, bool geometryUpdated)
{
  if (tp) {
    if (tp->step == m_last_step)
      return true; // ignore multiple calls for the same time step
    else if (tp->step % m_ndump > 0)
      return true; // write only every m_ndump step

    m_last_step = tp->step;
  }

  PROFILE1("DataExporter::dumpTimeLevel");

  if (m_level == -1)
    m_level = this->getWritersTimeLevel()+1;

  for (DataWriter* writer : m_writers) {
    writer->openFile(m_level);
    for (const DataEntry& it : m_entry) {
      if (!it.second.data)
        return false;
      switch (it.second.field) {
        case INTVECTOR:
        case VECTOR:
          writer->writeVector(m_level,it);
          break;
        case SIM:
          writer->writeSIM(m_level,it,geometryUpdated,it.second.prefix);
          break;
        case NODALFORCES:
          writer->writeNodalForces(m_level,it);
          break;
        case KNOTSPAN:
          writer->writeKnotspan(m_level,it,it.second.prefix);
          break;
        case BASIS:
          writer->writeBasis(m_level,it,it.second.prefix);
          break;
        default:
          std::cerr <<"  ** DataExporter: Invalid field type registered "
                    << it.second.field <<", skipping"<< std::endl;
          break;
      }
    }
    if (tp && tp->multiSteps())
      writer->writeTimeInfo(m_level,m_ndump,*tp);

    writer->closeFile(m_level);
  }
  m_level++;

  // disable fields marked as once
  for (std::pair<const std::string,FileEntry>& it : m_entry)
    if (abs(it.second.results) & ONCE)
      it.second.enabled = false;

  return true;
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
  for (DataWriter* writer : m_writers)
    levels.push_back(writer->getLastTimeLevel());
  return *min_element(levels.begin(),levels.end());
}


void DataExporter::setNormPrefixes(const std::vector<std::string>& prefix)
{
  for (DataWriter* writer : m_writers)
    writer->setNormPrefixes(prefix);
}


int DataExporter::realTimeLevel(int filelevel) const
{
  return filelevel*m_ndump;
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
  if (m_writers.empty())
    return "";

  const std::string& name = m_writers.front()->getName();
  return name.substr(0,name.rfind('.'));
}
