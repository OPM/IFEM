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

#include "tinyxml.h"
#include <algorithm>


void XMLInputBase::injectIncludeFiles (TiXmlElement* tag) const
{
  static int nLevels = 0; ++nLevels;
  bool foundIncludes = false;
  TiXmlElement* elem = tag->FirstChildElement();
  for (; elem; elem = elem->NextSiblingElement())
    if (strcasecmp(elem->Value(),"include"))
      this->injectIncludeFiles(elem);
    else if (elem->FirstChild() && elem->FirstChild()->Value()) {
      TiXmlDocument doc;
      if (doc.LoadFile(elem->FirstChild()->Value())) {
        for (int i = 1; i < nLevels; i++) IFEM::cout <<"  ";
        IFEM::cout <<"Loaded included file "<< elem->FirstChild()->Value()
                   << std::endl;
        elem = tag->ReplaceChild(elem,*doc.RootElement())->ToElement();
        TiXmlElement* elem2 = doc.RootElement()->NextSiblingElement();
        for (; elem2; elem2 = elem2->NextSiblingElement())
          tag->LinkEndChild(new TiXmlElement(*elem2));
        foundIncludes = true;
      }
      else
	std::cerr << __PRETTY_FUNCTION__ <<": Failed to load "
		  << elem->FirstChild()->Value()
		  <<"\n\tError at line "<< doc.ErrorRow() <<": "
		  << doc.ErrorDesc() << std::endl;
    }

  if (foundIncludes)
    this->injectIncludeFiles(tag);

  --nLevels;
}


bool XMLInputBase::readXML (const char* fileName)
{
  TiXmlDocument doc;
  if (!doc.LoadFile(fileName)) {
    std::cerr << __PRETTY_FUNCTION__ <<": Failed to load "<< fileName
	      <<"\n\tError at line "<< doc.ErrorRow() <<": "
	      << doc.ErrorDesc() << std::endl;
    return false;
  }

  const TiXmlElement* tag = doc.RootElement();
  if (!tag || strcmp(tag->Value(),"simulation")) {
    std::cerr << __PRETTY_FUNCTION__ <<": Malformatted input file "<< fileName
	      << std::endl;
    return false;
  }

  IFEM::cout <<"\nParsing input file "<< fileName << std::endl;

  this->injectIncludeFiles(const_cast<TiXmlElement*>(tag));

  std::vector<const TiXmlElement*> parsed;
  if (!handlePriorityTags(doc.RootElement(),parsed))
    return false;

  for (tag = tag->FirstChildElement(); tag; tag = tag->NextSiblingElement())
    if (std::find(parsed.begin(),parsed.end(),tag) == parsed.end()) {
      IFEM::cout <<"\nParsing <"<< tag->Value() <<">"<< std::endl;
      if (!this->parse(tag)) {
        std::cerr <<" *** SIMinput::read: Failure occured while parsing \""
                  << tag->Value() <<"\""<< std::endl;
        return false;
      }
    }

  IFEM::cout <<"\nParsing input file succeeded."<< std::endl;

  return true;
}


bool XMLInputBase::handlePriorityTags (const TiXmlElement* base,
				       std::vector<const TiXmlElement*>& parsed)
{
  const char** q = this->getPrioritizedTags();
  if (!q) return true; // No prioritized tags defined

  for (const TiXmlElement* elem = 0; *q; q++)
    if ((elem = base->FirstChildElement(*q))) {
      IFEM::cout <<"\nParsing <"<< elem->Value() <<">"<< std::endl;
      if (!this->parse(elem)) {
        std::cerr <<" *** SIMinput::read: Failure occured while parsing \""
                  << elem->Value() <<"\""<< std::endl;
        return false;
      }
      parsed.push_back(elem);
    }

  return true;
}
