// $Id$
//==============================================================================
//!
//! \file SIMinput.C
//!
//! \date Jun 1 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for simulators with input parsing functionality.
//!
//==============================================================================

#include "SIMinput.h"
#include "Utilities.h"
#include "tinyxml.h"
#ifdef PARALLEL_PETSC
#include "petscsys.h"
#endif
#include <fstream>


int SIMinput::msgLevel = 2;


SIMinput::SIMinput ()
{
#ifdef PARALLEL_PETSC
  MPI_Comm_rank(PETSC_COMM_WORLD,&myPid);
  MPI_Comm_size(PETSC_COMM_WORLD,&nProc);
#else
  myPid = 0;
  nProc = 1;
#endif
}


bool SIMinput::read (const char* fileName)
{
  if (strcasestr(fileName,".xinp"))
    return this->readXML(fileName);
  else
    return this->readFlat(fileName);
}


bool SIMinput::readFlat (const char* fileName)
{
  std::ifstream is(fileName);
  if (is)
    std::cout <<"\nReading input file "<< fileName << std::endl;
  else
  {
    std::cerr <<"\n *** SIMinput::read: Failure opening input file "
	      << fileName << std::endl;
    return false;
  }

  char* keyWord = 0;
  while (is.good() && (keyWord = utl::readLine(is)))
    if (!this->parse(keyWord,is))
    {
      std::cerr <<" *** SIMinput::read: Failure occured while parsing \""
		<< keyWord <<"\""<< std::endl;
      return false;
    }

  std::cout <<"\nReading input file succeeded."<< std::endl;
  return true;
}


/*!
  \brief Recursive helper method for processing the \a include XML-tags.
*/

static void injectIncludeFiles (TiXmlElement* tag)
{
  TiXmlElement* elem = tag->FirstChildElement();
  while (elem) {
    if (strcasecmp(elem->Value(),"include"))
      injectIncludeFiles(elem);
    else if (elem->FirstChild() && elem->FirstChild()->Value()) {
      TiXmlDocument doc;
      if (doc.LoadFile(elem->FirstChild()->Value()))
	tag->ReplaceChild(elem,*doc.RootElement());
      else
	std::cerr << __PRETTY_FUNCTION__ <<": Failed to load "
		  << elem->FirstChild()->Value()
		  <<"\n\tError at line "<< doc.ErrorRow() <<": "
		  << doc.ErrorDesc() << std::endl;
    }
    elem = elem->NextSiblingElement();
  }
}


bool SIMinput::readXML (const char* fileName)
{
  TiXmlDocument doc;
  if (!doc.LoadFile(fileName)) {
    std::cerr << __PRETTY_FUNCTION__ <<": Failed to load "<< fileName
	      <<"\n\tError at line "<< doc.ErrorRow() <<": "
	      << doc.ErrorDesc() << std::endl;
    return false;
  }

  if (!doc.RootElement() || strcmp(doc.RootElement()->Value(),"simulation")) {
    std::cerr << __PRETTY_FUNCTION__ <<": Malformatted input file "<< fileName
	      << std::endl;
    return false;
  }

  injectIncludeFiles(doc.RootElement());

  std::vector<const TiXmlElement*> parsed;
  if (!handlePriorityTags(doc.RootElement(),parsed))
    return false;

  TiXmlElement* elem = doc.RootElement()->FirstChildElement();
  while (elem) {
    if (std::find(parsed.begin(),parsed.end(),elem) == parsed.end()) {
      if (!this->parse(elem)) {
        std::cerr <<" *** SIMinput::read: Failure occured while parsing \""
                  << elem->Value() <<"\""<< std::endl;
        return false;
      }
    }
    elem = elem->NextSiblingElement();
  }
  return true;
}


bool SIMinput::handlePriorityTags (const TiXmlElement* base,
				   std::vector<const TiXmlElement*>& parsed)
{
  // Add particular keywords to this list.
  // They will be parsed first, in the order specified in the list.
  const char* special[] = {"discretization","geometry","boundaryconditions",0};

  const TiXmlElement* elem = 0;
  for (const char** q = special; *q; q++)
    if ((elem = base->FirstChildElement(*q))) {
      if (!this->parse(elem)) {
        std::cerr <<" *** SIMinput::read: Failure occured while parsing \""
                  << elem->Value() <<"\""<< std::endl;
	return false;
      }
      parsed.push_back(elem);
    }

  return true;
}
