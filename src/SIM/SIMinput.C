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
#include "IFEM.h"
#include "Utilities.h"
#include "tinyxml.h"
#include "SystemMatrix.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#include <fstream>

#if defined(__MINGW32__) || defined(__MINGW64__)
#define strcasestr strstr
#endif

int SIMinput::msgLevel = 2;


SIMinput::SIMinput (const char* heading) : opt(myOpts)
#ifdef HAS_PETSC
  , adm(PETSC_COMM_WORLD)
#elif defined(HAVE_MPI)
  , adm(true)
#endif
{
  myPid = adm.getProcId();
  nProc = adm.getNoProcs();

  myOpts = IFEM::getOptions(); // Initialize options from command-line arguments

  if (heading) myHeading = heading;
}


SIMinput::SIMinput (SIMinput& anotherSIM) : opt(anotherSIM.myOpts)
{
  adm = anotherSIM.adm;
  myPid = anotherSIM.myPid;
  nProc = anotherSIM.nProc;
}


void SIMinput::printHeading (int& subStep) const
{
  if (myHeading.empty())
    return;

  size_t n = myHeading.find_last_of('\n');
  if (n+1 < myHeading.size()) n = myHeading.size()-n;
  IFEM::cout <<"\n"<< ++subStep <<". "<< myHeading <<"\n";
  for (size_t i = 0; i < 3+n && i < 3+myHeading.size(); i++) IFEM::cout <<'=';
  if (subStep > 9) IFEM::cout <<'=';
  IFEM::cout << std::endl;
}


bool SIMinput::read (const char* fileName)
{
  static int substep = 0;
  this->printHeading(substep);

  bool result;
  if (strcasestr(fileName,".xinp"))
    result = this->readXML(fileName);
  else
    result = this->readFlat(fileName);

  // Let command-line options override settings on the input file
  IFEM::applyCommandLineOptions(opt);

  return result;
}


bool SIMinput::readFlat (const char* fileName)
{
  std::ifstream is(fileName);
  if (!is)
  {
    std::cerr <<"\n *** SIMinput::read: Failure opening input file "
              << fileName << std::endl;
    return false;
  }

  IFEM::cout <<"\nReading input file "<< fileName << std::endl;

  char* keyWord = 0;
  while (is.good() && (keyWord = utl::readLine(is)))
    if (!this->parse(keyWord,is))
    {
      std::cerr <<" *** SIMinput::read: Failure occured while parsing \""
                << keyWord <<"\""<< std::endl;
      return false;
    }

  IFEM::cout <<"\nReading input file succeeded."<< std::endl;

  return true;
}


bool SIMinput::parse (char*, std::istream&)
{
  std::cerr <<" *** SIMinput::parse(char*,std::istream&):"
            <<" The flat file format is depreciated.\n"
            <<"     Use the XML format instead."<< std::endl;
  return false;
}


void SIMinput::injectIncludeFiles (TiXmlElement* tag) const
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


bool SIMinput::readXML (const char* fileName)
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


bool SIMinput::handlePriorityTags (const TiXmlElement* base,
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
