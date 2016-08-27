// $Id$
//==============================================================================
//!
//! \file XMLWriter.C
//!
//! \date May 7 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Output of metadata associated with HDF5 to XML.
//!
//==============================================================================

#include "XMLWriter.h"
#include "GlbForceVec.h"
#include "SIMbase.h"
#include "IntegrandBase.h"
#include "TimeStep.h"
#include "StringUtils.h"
#include "tinyxml.h"
#include <fstream>
#include <cstdio>


XMLWriter::XMLWriter (const std::string& name, const ProcessAdm& adm) :
  DataWriter(name,adm,".xml")
{
  m_doc = nullptr;
  m_node = nullptr;
  m_dt = 0;
  m_order = m_interval = 1;
}


int XMLWriter::getLastTimeLevel()
{
  TiXmlDocument doc(m_name.c_str());
  doc.LoadFile();
  TiXmlHandle handle(&doc);
  TiXmlElement* levels = handle.FirstChild("info").
                                            FirstChild("levels").ToElement();
  if (levels && levels->FirstChild())
    return atoi(levels->FirstChild()->Value());

  return -1;
}


void XMLWriter::openFile(int level)
{
  if (m_rank != 0)
    return;

  m_doc = new TiXmlDocument;
  TiXmlElement element("info");
  m_node = m_doc->InsertEndChild(element);
}


void XMLWriter::closeFile(int level, bool force)
{
  if (!m_doc || m_rank != 0)
    return;

  TiXmlElement element2("levels");
  TiXmlNode *pNewNode = m_node->InsertEndChild(element2);
  char temp[32];
  sprintf(temp,"%i",level);
  TiXmlText value(temp);
  pNewNode->InsertEndChild(value);

  // TODO: support variable time steps
  TiXmlElement element3("timestep");
  sprintf(temp,"%f",m_dt);
  element3.SetAttribute("constant","1");
  element3.SetAttribute("order",m_order);
  element3.SetAttribute("interval",m_interval);
  pNewNode = m_node->InsertEndChild(element3);
  TiXmlText value2(temp);
  pNewNode->InsertEndChild(value2);

  m_doc->SaveFile(m_name);
  delete m_doc;
  m_doc = nullptr;
}


void XMLWriter::readInfo()
{
  TiXmlDocument doc(m_name);
  doc.LoadFile();
  TiXmlHandle handle(&doc);
  TiXmlElement* elem = handle.FirstChild("info").
                                         FirstChild("entry").ToElement();
  TiXmlElement* timestep = handle.FirstChild("info").FirstChild("timestep").ToElement();
  while (elem) {
    if (strcasecmp(elem->Attribute("type"),"field") == 0 ||
        strcasecmp(elem->Attribute("type"),"knotspan") == 0 ||
        strcasecmp(elem->Attribute("type"),"displacement") == 0 ||
        strcasecmp(elem->Attribute("type"), "eigenmodes") == 0 ||
        strcasecmp(elem->Attribute("type"), "nodalforces") == 0) {
      Entry entry;
      entry.name = elem->Attribute("name");
      entry.description = elem->Attribute("description");
      if (elem->Attribute("patches"))
        entry.patches = atoi(elem->Attribute("patches"));
      else
        entry.patches = 0;
      if (elem->Attribute("components"))
        entry.components = atoi(elem->Attribute("components"));
      else
        entry.components = 0;
      entry.type = elem->Attribute("type");
      if (elem->Attribute("once") && strcasecmp(elem->Attribute("once"),"true") == 0)
        entry.once = true;
      else
        entry.once = false;
      if (timestep)
        entry.timestep = atof(timestep->FirstChild()->Value());
      else
        entry.timestep = 0;
      if (elem->Attribute("basis"))
        entry.basis = elem->Attribute("basis");
      m_entry.push_back(entry);
    }
    elem = elem->NextSiblingElement("entry");
  }
}


bool XMLWriter::readVector(int level, const DataEntry& entry)
{
  return true;
}


void XMLWriter::writeVector(int level, const DataEntry& entry)
{
  if (m_rank != 0)
    return;

  TiXmlElement element("entry");
  element.SetAttribute("name",entry.first.c_str());
  element.SetAttribute("description",entry.second.description.c_str());
  if (entry.second.field == DataExporter::INTVECTOR)
    element.SetAttribute("type","intvector");
  else
    element.SetAttribute("type","vector");
  Vector* vec = (Vector*)entry.second.data;
  element.SetAttribute("size",vec->size());
  m_node->InsertEndChild(element);
}


void XMLWriter::writeNodalForces(int level, const DataEntry& entry)
{
  if (m_rank != 0)
    return;

  TiXmlElement element("entry");
  element.SetAttribute("name",entry.first.c_str());
  element.SetAttribute("description",entry.second.description.c_str());
  element.SetAttribute("type","nodalforces");
  GlbForceVec* vec = (GlbForceVec*)entry.second.data;
  element.SetAttribute("size",vec->size());
  m_node->InsertEndChild(element);
}


void XMLWriter::writeKnotspan(int level, const DataEntry& entry,
                              const std::string& prefix)
{
  if (m_rank != 0)
    return;

  const SIMbase* sim = static_cast<const SIMbase*>(entry.second.data);
  std::string basisname;
  if (prefix.empty())
    basisname = sim->getName()+"-1";
  else
    basisname = prefix+sim->getName()+"-1";

  addField(prefix+entry.first,entry.second.description,
           basisname,1,sim->getNoPatches(),"knotspan");
}


bool XMLWriter::readSIM (int level, const DataEntry& entry)
{
  return true;
}


void XMLWriter::writeSIM (int level, const DataEntry& entry, bool,
                          const std::string& prefix)
{
  if (m_rank != 0)
    return;

  const SIMbase* sim = static_cast<const SIMbase*>(entry.second.data);
  const IntegrandBase* prob = sim->getProblem();

  int results = entry.second.results;
  bool usedescription=false;
  if (results < 0) {
    results = -results;
    usedescription = true;
  }

  int cmps = entry.second.ncmps>0?entry.second.ncmps:prob->getNoFields(1);

  std::string basisname;
  if (prefix.empty())
    basisname = sim->getName()+"-1";
  else
    basisname = prefix+sim->getName()+"-1";

  // restart vector
  if (results & DataExporter::RESTART) {
    addField(prefix+"restart",entry.second.description,basisname,
             cmps,sim->getNoPatches(),"restart");
  }

  if (results & DataExporter::PRIMARY) {
    if (entry.second.results < 0)
      addField(entry.second.description, entry.second.description,
               basisname, cmps, sim->getNoPatches(), "field");
    else if (sim->mixedProblem())
    {
      // primary solution vector
      addField(prefix+entry.first,entry.second.description,basisname,
               0,sim->getNoPatches(),"restart");

      // primary solution fields
      for (int b=1; b <= sim->getNoBasis(); ++b) {
        std::stringstream str;
        str << sim->getName() << "-" << b;
        addField(prefix+prob->getField1Name(10+b),"primary",str.str(),
                 sim->getNoFields(b),sim->getNoPatches(),
                 "field",results & DataExporter::ONCE?true:false);
      }
    }
    else
    {
      // primary solution
      addField(usedescription ? entry.second.description:
                                prefix+prob->getField1Name(11),
               entry.second.description,basisname,
               cmps,sim->getNoPatches(),"field",results & DataExporter::ONCE?true:false);
    }
  }

  if (results & DataExporter::DISPLACEMENT) {
    if (sim->mixedProblem())
    {
      // primary solution vector
      addField(prefix+entry.first,entry.second.description,basisname,
               prob->getNoFields(1),sim->getNoPatches(),"displacement");
    }
    else
    {
      // primary solution
      addField(usedescription ? entry.second.description:
                                prefix+prob->getField1Name(11),
               entry.second.description,basisname,
               prob->getNoFields(1),sim->getNoPatches(), "displacement");
    }
  }

  // secondary solution fields
  size_t i, j;
  if (results & DataExporter::SECONDARY)
    for (j = 0; j < prob->getNoFields(2); j++)
      addField(prefix+prob->getField2Name(j),"secondary", basisname,
               1,sim->getNoPatches());

  // norms
  if (results & DataExporter::NORMS) {
    // since the norm data isn't available, we have to instance the object
    NormBase* norm = sim->getNormIntegrand();
    if (norm) {
      for (i = 1; i <= norm->getNoFields(0); ++i) {
        for (j = 1; j <= norm->getNoFields(i); ++j) {
          if (norm->hasElementContributions(i,j))
            addField(prefix+norm->getName(i,j,(i>1&&m_prefix)?m_prefix[i-2]:0),"knotspan wise norm",
                     basisname,1,sim->getNoPatches(),"knotspan");
        }
      }
      delete norm;
    }
  }

  if (results & DataExporter::EIGENMODES) {
    const std::vector<Mode>* vec = static_cast<const std::vector<Mode>* >(entry.second.data2);
    addField(prefix+"eigenmode", "eigenmode", sim->getName()+"-1", vec->size(), sim->getNoPatches(), "eigenmodes");
  }
}


bool XMLWriter::writeTimeInfo (int level, int order, int interval,
                               const TimeStep& tp)
{
  m_dt = tp.time.dt;
  m_order = order;
  m_interval = interval;
  return true;
}


void XMLWriter::addField (const std::string& name,
                          const std::string& description,
                          const std::string& basis, int components, int patches,
                          const std::string& type, bool once)
{
  TiXmlElement element("entry");
  element.SetAttribute("name",name.c_str());
  element.SetAttribute("description",description.c_str());
  element.SetAttribute("type",type);
  if (!basis.empty())
    element.SetAttribute("basis",basis.c_str());
  element.SetAttribute("patches",patches);
  element.SetAttribute("components",components);
  if (once)
    element.SetAttribute("once","true");

  m_node->InsertEndChild(element);
}
