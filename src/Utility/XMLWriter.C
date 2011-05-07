// $Id$

#include "XMLWriter.h"
#include "SIMbase.h"
#include "IntegrandBase.h"

#ifdef PARALLEL_PETSC
#include <mpi.h>
#endif


XMLWriter::XMLWriter(const std::string& name) : m_xml(name+".xml")
{
  m_doc = NULL;
#ifdef PARALLEL_PETSC
  MPI_Comm_rank(MPI_COMM_WORLD,&m_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&m_size);
#else
  m_rank = 0;
#endif
}

XMLWriter::~XMLWriter()
{
}

int XMLWriter::getLastTimeLevel()
{
  int result=-1;
  TiXmlDocument doc(m_xml.c_str());
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

  m_doc->SaveFile(m_xml);
  delete m_doc;
  m_doc = NULL;
}

void XMLWriter::readVector(int level, const DataEntry& entry)
{
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

void XMLWriter::readSIM(int level, const DataEntry& entry)
{
}

void XMLWriter::writeSIM(int level, const DataEntry& entry)
{
  if (m_rank != 0)
    return;

  int nPatch = static_cast<SIMbase*>(entry.second.data)->getNoPatches();
  const Integrand* prb = static_cast<SIMbase*>(entry.second.data)->getProblem();

  // primary solution
  addField(entry.first,entry.second.description,"field",nPatch);

  // secondary solutions
  if (entry.second.size == -1)
    for (size_t j = 0; j < prb->getNoFields(2); j++)
      addField(prb->getField2Name(j),"secondary","field",nPatch);
}

void XMLWriter::addField(const std::string& name, const std::string& description,
                         const std::string& type, int patches)
{
  TiXmlElement element("entry");
  element.SetAttribute("name",name.c_str());
  element.SetAttribute("description",description.c_str());
  element.SetAttribute("type",type.c_str());
  element.SetAttribute("patches",patches);
  m_node->InsertEndChild(element);
}
