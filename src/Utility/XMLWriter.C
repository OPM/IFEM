// $Id$

#include "XMLWriter.h"
#include "SIMbase.h"
#include "IntegrandBase.h"
#include "StringUtils.h"
#include "tinyxml.h"
#include <fstream>
#include <cstdio>


XMLWriter::XMLWriter(const std::string& name) : DataWriter(name+".xml")
{
  m_doc = NULL;
  m_node = NULL;
}


int XMLWriter::getLastTimeLevel()
{
  int result=-1;
  TiXmlDocument doc(m_name.c_str());
  doc.LoadFile();
  TiXmlHandle handle(&doc);
  TiXmlElement* levels = handle.FirstChild("info").
                                            FirstChild("levels").ToElement();
  if (levels && levels->FirstChild())
    result = atoi(levels->FirstChild()->Value());

  return result;
}


void XMLWriter::openFile(int level)
{
  if (m_rank != 0)
    return;

  m_doc = new TiXmlDocument;
  TiXmlElement element("info");
  m_node = m_doc->InsertEndChild(element);
}

void XMLWriter::closeFile(int level)
{
  if (!m_doc || m_rank != 0)
    return;

  TiXmlElement element2("levels");
  TiXmlNode *pNewNode = m_node->InsertEndChild(element2);
  char temp[8];
  sprintf(temp,"%i",level);
  TiXmlText value(temp);
  pNewNode->InsertEndChild(value);

  m_doc->SaveFile(m_name);
  delete m_doc;
  m_doc = NULL;
}


void XMLWriter::readInfo()
{
  TiXmlDocument doc(m_name);
  doc.LoadFile();
  TiXmlHandle handle(&doc);
  TiXmlElement* elem = handle.FirstChild("info").
                                         FirstChild("entry").ToElement();
  while (elem) {
    if (strcasecmp(elem->Attribute("type"),"field") == 0) {
      Entry entry;
      entry.name = elem->Attribute("name");
      entry.description = elem->Attribute("description");
      entry.patches = atoi(elem->Attribute("patches"));
      entry.components = atoi(elem->Attribute("components"));
      if (elem->Attribute("patchfile"))
	entry.patchfile = elem->Attribute("patchfile");
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
  element.SetAttribute("type","vector");
  element.SetAttribute("size",entry.second.size);
  m_node->InsertEndChild(element);
}


bool XMLWriter::readSIM (int level, const DataEntry& entry)
{
  return true;
}

void XMLWriter::writeSIM (int level, const DataEntry& entry)
{
  if (m_rank != 0)
    return;

  const SIMbase*   sim  = static_cast<const SIMbase*>(entry.second.data);
  const Integrand* prob = sim->getProblem();

  std::string g2file;
  if (prob->mixedFormulation())
  {
    // primary solution vector
    addField(entry.first,entry.second.description,sim->getName()+"-0",
	     prob->getNoFields(1),sim->getNoPatches());

    // Assuming that basis2 is used for secondary variables
    // primary solution fields
    addField(prob->getField1Name(11),"primary",sim->getName()+"-0",
	     sim->getNoFields(1),sim->getNoPatches());
    addField(prob->getField1Name(12),"primary",sim->getName()+"-1",
	     sim->getNoFields(2),sim->getNoPatches());
  }
  else
  {
    // primary solution
    addField(prob->getField1Name(11),entry.second.description,sim->getName()+"-0",
	     prob->getNoFields(1),sim->getNoPatches());
  }

  // secondary solution fields
  if (entry.second.size == -1)
    for (size_t j = 0; j < prob->getNoFields(2); j++)
      addField(prob->getField2Name(j),"secondary",sim->getName()+"-1",1,sim->getNoPatches());
}


void XMLWriter::addField (const std::string& name,
                          const std::string& description,
                          const std::string& patchfile,
                          int components, int patches)
{
  TiXmlElement element("entry");
  element.SetAttribute("name",name.c_str());
  element.SetAttribute("description",description.c_str());
  element.SetAttribute("type","field");
  if (!patchfile.empty())
    element.SetAttribute("patchfile",patchfile.c_str());
  element.SetAttribute("patches",patches);
  element.SetAttribute("components",components);
  m_node->InsertEndChild(element);
}
