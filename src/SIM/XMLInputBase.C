// $Id$
//==============================================================================
//!
//! \file XMLInputBase.C
//!
//! \date Jul 16 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Base class for XML input parsing functionality.
//!
//==============================================================================

#include "XMLInputBase.h"
#include "IFEM.h"
#include "tinyxml2.h"
#include <algorithm>
#include <cstring>

namespace {

/*!
  \brief Helper method to load an XML file and print error message if failure.
*/

const tinyxml2::XMLElement* loadXMLfile (tinyxml2::XMLDocument& doc,
                                         const char* fileName)
{
  tinyxml2::XMLError err = doc.LoadFile(fileName);
  if (err == tinyxml2::XML_SUCCESS)
    return doc.RootElement();

  std::cerr <<" *** Failed to load XML-file \""<< fileName
            <<"\".\n" << doc.ErrorStr() << std::endl;
  return nullptr;
}

/*!
  \brief Helper to handle includes.
*/

struct IncludeInjector : public tinyxml2::XMLVisitor
{
  //! \brief The constructor creates the new document.
  IncludeInjector() : new_doc(true, tinyxml2::COLLAPSE_WHITESPACE) {}

  //! \brief If elem is include, replace with contents of include file, else copy element.
  bool VisitEnter(const tinyxml2::XMLElement& elem,
                  const tinyxml2::XMLAttribute* attribute) override
  {
      tinyxml2::XMLElement* e;
      if (std::string(elem.Value()) == "include") {
        tinyxml2::XMLDocument inc_doc(true, tinyxml2::COLLAPSE_WHITESPACE);
          if (inc_doc.LoadFile(elem.GetText()) != tinyxml2::XML_SUCCESS) {
            std::cerr << "** Error parsing include file " << elem.GetText() << std::endl;
            return false;
        }
        currElem->InsertEndChild(inc_doc.RootElement()->DeepClone(&new_doc));
        include_found = true;
      } else {
        e = new_doc.NewElement(elem.Name());
        if (elem.GetText())
          e->SetText(elem.GetText());
        if (!currElem)
          new_doc.InsertEndChild(e);
        else
          currElem->InsertEndChild(e);
        while (attribute)
        {
          e->SetAttribute(attribute->Name(), attribute->Value());
          attribute = attribute->Next();
        }
        currElem = e;
      }
      return true;
  }

  //! brief Set element to insert into to the parent.
  bool VisitExit(const tinyxml2::XMLElement& elem) override
  {
      if (currElem && currElem->Parent())
        currElem = const_cast<tinyxml2::XMLElement*>(currElem->Parent()->ToElement());
      return true;
  }

  tinyxml2::XMLDocument new_doc; //!< New document
  tinyxml2::XMLElement* currElem = nullptr; //!< Current root element
  bool include_found = false; //!< True if an include tag was found
};

}


const tinyxml2::XMLElement*
XMLInputBase::loadFile (tinyxml2::XMLDocument& doc,
                        const char* fileName, bool verbose)
{
  const tinyxml2::XMLElement* tag = loadXMLfile(doc,fileName);

  if (tag && strcmp(tag->Value(),"simulation")) {
    std::cerr <<" *** Malformatted XML-file \""<< fileName <<"\"."<< std::endl;
    verbose = false;
    tag = nullptr;
  }

  if (verbose)
    IFEM::cout <<"\nParsing input file "<< fileName << std::endl;

  if (tag) {
    for (size_t i = 0; i < 10; ++i) { // Maximum 10 levels of include files
      IncludeInjector v;
      doc.RootElement()->Accept(&v);
      if (v.include_found)
        v.new_doc.DeepCopy(&doc);
      else
        break;
    }
    tag = doc.RootElement();
  }

#ifdef SP_DEBUG
  if (verbose) {
    std::cout <<"\nHere is the input-file content:"<< std::endl;
    doc.Print();
  }
#endif

  return tag;
}


bool XMLInputBase::readXML (const char* fileName, bool verbose)
{
  tinyxml2::XMLDocument doc(true, tinyxml2::COLLAPSE_WHITESPACE);
  const tinyxml2::XMLElement* tag = this->loadFile(doc,fileName,verbose);
  if (!tag) return false;

  std::vector<const tinyxml2::XMLElement*> parsed;
  if (!this->handlePriorityTags(tag,parsed,verbose))
    return false;

  for (tag = tag->FirstChildElement(); tag; tag = tag->NextSiblingElement())
    if (std::find(parsed.begin(),parsed.end(),tag) == parsed.end()) {
      if (verbose)
        IFEM::cout <<"\nParsing <"<< tag->Value() <<">"<< std::endl;
      if (!this->parse(tag)) {
        std::cerr <<" *** SIMadmin::read: Failure occured while parsing \""
                  << tag->Value() <<"\"."<< std::endl;
        return false;
      }
    }

  if (verbose)
    IFEM::cout <<"\nParsing input file succeeded."<< std::endl;

  return true;
}


bool XMLInputBase::loadXML (const char* xml)
{
  tinyxml2::XMLDocument doc;
  doc.Parse(xml);
  const tinyxml2::XMLElement* tag = doc.RootElement();
  return tag ? this->parse(tag) : false;
}


bool XMLInputBase::handlePriorityTags (const tinyxml2::XMLElement* base,
                                       std::vector<const tinyxml2::XMLElement*>& parsed,
                                       bool verbose)
{
  const char** q = this->getPrioritizedTags();
  if (!q) return true; // No prioritized tags defined

  while (*q) {
    const tinyxml2::XMLElement* elem = base->FirstChildElement(*(q++));
    if (elem) {
      if (verbose)
        IFEM::cout <<"\nParsing <"<< elem->Value() <<">"<< std::endl;
      if (!this->parse(elem)) {
        std::cerr <<" *** SIMadmin::read: Failure occured while parsing \""
                  << elem->Value() <<"\"."<< std::endl;
        return false;
      }
      parsed.push_back(elem);
    }
  }

  return true;
}
