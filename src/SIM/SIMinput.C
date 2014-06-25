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
#include "ASMbase.h"
#ifdef PARALLEL_PETSC
#include "petscsys.h"
#endif
#include <fstream>

#if defined(__MINGW32__) || defined(__MINGW64__)
#define strcasestr strstr
#endif

int SIMinput::msgLevel = 2;


SIMinput::SIMinput (const char* heading) : opt(myOpts)
#ifdef PARALLEL_PETSC
  , adm(PETSC_COMM_WORLD)
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
  if (myHeading.empty() || myPid > 0)
    return;

  size_t n = myHeading.find_last_of('\n');
  if (n+1 < myHeading.size()) n = myHeading.size()-n;
  std::cout <<"\n"<< ++subStep <<". "<< myHeading <<"\n";
  for (size_t i = 0; i < 3+n && i < 3+myHeading.size(); i++) std::cout <<'=';
  if (subStep > 9) std::cout <<'=';
  std::cout << std::endl;
}


bool SIMinput::read (const char* fileName)
{
#ifdef HAS_PETSC
  // In parallel simulations, we need to retain all DOFs in the equation system.
  // The fixed DOFs (if any) will receive a homogeneous constraint instead.
  if (opt.solver == SystemMatrix::PETSC)
    ASMbase::fixHomogeneousDirichlet = false;
#endif

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

  if (myPid == 0)
    std::cout <<"\nReading input file "<< fileName << std::endl;

  char* keyWord = 0;
  while (is.good() && (keyWord = utl::readLine(is)))
    if (!this->parse(keyWord,is))
    {
      std::cerr <<" *** SIMinput::read: Failure occured while parsing \""
                << keyWord <<"\""<< std::endl;
      return false;
    }

  if (myPid == 0)
    std::cout <<"\nReading input file succeeded."<< std::endl;

  return true;
}


bool SIMinput::parse (char*, std::istream&)
{
  std::cerr <<" *** SIMinput::parse(char*,std::istream&):"
            <<" The flat file format is depreciated.\n"
            <<"     Use the XML format instead."<< std::endl;
  return false;
}


/*!
  \brief Recursive helper method for processing the \a include XML-tags.
*/

static bool injectIncludeFiles (TiXmlElement* tag)
{
  bool result=false;
  TiXmlElement* elem = tag->FirstChildElement();
  for (; elem; elem = elem->NextSiblingElement())
    if (strcasecmp(elem->Value(),"include"))
      result |= injectIncludeFiles(elem);
    else if (elem->FirstChild() && elem->FirstChild()->Value()) {
      TiXmlDocument doc;
      if (doc.LoadFile(elem->FirstChild()->Value())) {
        elem = tag->ReplaceChild(elem,*doc.RootElement())->ToElement();
        TiXmlElement* elem2 = doc.RootElement()->NextSiblingElement();
        while (elem2) {
          tag->LinkEndChild(new TiXmlElement(*elem2));
          elem2 = elem2->NextSiblingElement();
        }
        result = true;
      } else
	std::cerr << __PRETTY_FUNCTION__ <<": Failed to load "
		  << elem->FirstChild()->Value()
		  <<"\n\tError at line "<< doc.ErrorRow() <<": "
		  << doc.ErrorDesc() << std::endl;
    }
  return result;
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

  if (myPid == 0)
    std::cout <<"\nParsing input file "<< fileName << std::endl;

  while (injectIncludeFiles(const_cast<TiXmlElement*>(tag)));

  std::vector<const TiXmlElement*> parsed;
  if (!handlePriorityTags(doc.RootElement(),parsed))
    return false;

  for (tag = tag->FirstChildElement(); tag; tag = tag->NextSiblingElement())
    if (std::find(parsed.begin(),parsed.end(),tag) == parsed.end()) {
      if (myPid == 0)
        std::cout <<"\nParsing <"<< tag->Value() <<">"<< std::endl;
      if (!this->parse(tag)) {
        std::cerr <<" *** SIMinput::read: Failure occured while parsing \""
                  << tag->Value() <<"\""<< std::endl;
        return false;
      }
    }

  if (myPid == 0)
    std::cout <<"\nParsing input file succeeded."<< std::endl;

  return true;
}


bool SIMinput::handlePriorityTags (const TiXmlElement* base,
				   std::vector<const TiXmlElement*>& parsed)
{
  const char** q = this->getPrioritizedTags();
  if (!q) return true; // No prioritized tags defined

  for (const TiXmlElement* elem = 0; *q; q++)
    if ((elem = base->FirstChildElement(*q))) {
      if (myPid == 0)
        std::cout <<"\nParsing <"<< elem->Value() <<">"<< std::endl;
      if (!this->parse(elem)) {
        std::cerr <<" *** SIMinput::read: Failure occured while parsing \""
                  << elem->Value() <<"\""<< std::endl;
        return false;
      }
      parsed.push_back(elem);
    }

  return true;
}
