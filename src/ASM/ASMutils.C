// $Id$
//==============================================================================
//!
//! \file ASMutils.C
//!
//! \date Aug 12 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Various utilities for assembly scope.
//!
//==============================================================================

#include "ASMutils.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "Vec3.h"
#include "tinyxml.h"
#include <sstream>
#include <cstring>
#include <cctype>


/*!
  This function reads a matlab file with arrays defining a standard FE mesh.
  The format corresponds to that of the output method SIMoutput::dumpMatlabGrid.
*/

bool ASM::readMatlab (std::istream& is, IntMat& MNPC, std::vector<Vec3>& nodes,
                      std::vector<NodeSet>& nodeSets)
{
  // Lambda function reading one data line ignoring all whitespace characters.
  auto&& readLine = [&is](std::string& s, char delim = '\n')
  {
    size_t i, n = 0;
    while (n == 0 && std::getline(is,s,delim))
      for (i = 0; i < s.length(); i++)
        if (!isspace(s[i])) s[n++] = s[i];

    s = s.substr(0,n);
    for (i = 0; i < n; i++)
      if (s[i] == ',') s[i] = ' ';

    return n > 0;
  };

  std::string cline;
  if (!readLine(cline,'=') ||
      cline.find("function") == std::string::npos ||
      cline.find("Node")     == std::string::npos ||
      cline.find("Element")  == std::string::npos)
  {
    std::cerr <<" *** ASM::readMatlab: Not a matlab file."<< std::endl;
    return false;
  }

  std::getline(is,cline);
  while (readLine(cline,'['))
    if (cline.find("Node=") == cline.length()-5)
    {
      while (readLine(cline))
      {
        size_t id; Vec3 X;
        std::istringstream(cline) >> id >> X;
#if SP_DEBUG > 1
        std::cout << id <<": "<< X << std::endl;
#endif
        if (id > nodes.size()) nodes.resize(id);
        nodes[id-1] = X;

        if (cline.back() == ';') break;
      }
#ifdef SP_DEBUG
      std::cout <<"Read "<< nodes.size() <<" nodes."<< std::endl;
#endif
    }
    else if (cline.find("Element=") == cline.length()-8)
    {
      while (readLine(cline))
      {
        size_t id; int node;
        IntVec mnpc; mnpc.reserve(4);
        std::istringstream selem(cline);
        selem >> id >> node;
        while (selem)
        {
          mnpc.push_back(node-1);
          selem >> node;
        }
#if SP_DEBUG > 1
        std::cout << id <<":";
        for (int n : mnpc) std::cout <<" "<< n;
        std::cout << std::endl;
#endif
        if (mnpc.size() == 4)
          std::swap(mnpc[2],mnpc[3]);
        else if (mnpc.size() != 2)
        {
          std::cerr <<"  ** ASM::readMatlab: "<< mnpc.size()
                    <<"-noded elements not supported (ignored).\n";
          continue;
        }
        if (id > MNPC.size()) MNPC.resize(id);
        MNPC[id-1] = mnpc;
        if (cline.back() == ';') break;
      }
#ifdef SP_DEBUG
      std::cout <<"Read "<< MNPC.size() <<" elements."<< std::endl;
#endif
    }
    else if (cline.back() == '=')
    {
      std::string setname = cline.substr(0,cline.length()-1);
      if (readLine(cline,']'))
      {
        size_t i = 0; // Remove the '...' continuation markers
        while ((i = cline.find_first_of('.',i)) != std::string::npos)
          cline.erase(i,1);

        int node;
        IntVec nodes;
        std::istringstream selem(cline);
        selem >> node;
        while (selem)
        {
          nodes.push_back(node);
          selem >> node;
        }
#if SP_DEBUG > 1
        std::cout <<"Node set \""<< setname <<"\":";
        for (int n : nodes) std::cout <<" "<< n;
        std::cout << std::endl;
#endif
        std::getline(is,cline);
        nodeSets.push_back(std::make_pair(setname,nodes));
      }
    }

  return true;
}


bool ASM::readXML (std::istream& is, IntMat& MNPC, std::vector<Vec3>& nodes,
                   std::vector<NodeSet>& nodeSets,
                   std::vector<NodeSet>* elemSets)
{
  char cline[128];
  if (!is.getline(cline,128))
    return false;
  else if (!strstr(cline,"<patch>"))
  {
    std::cerr <<" *** ASM::readXML: Failed to read patch geometry."<< std::endl;
    return false;
  }

  std::string data("<patch>\n");
  while (is.getline(cline,128))
  {
    data.append(cline);
    data.append("\n");
    if (strstr(cline,"</patch>"))
      break;
  }

  TiXmlDocument doc;
  doc.Parse(data.c_str(),nullptr,TIXML_ENCODING_UTF8);
  const TiXmlElement* tag = doc.RootElement();
  if (!tag)
  {
    std::cerr <<" *** ASM::readXML: Malformatted XML input."<< std::endl;
    return false;
  }

  // Lambda function for parsing nodal points from a string.
  auto&& parseNodes = [&nodes](const char* data)
  {
    Vec3 X;
    std::istringstream iss(data);
#if SP_DEBUG > 1
    for (size_t inod = 0; iss; inod++)
#else
    while (iss)
#endif
    {
      iss >> X;
      if (iss)
      {
#if SP_DEBUG > 1
        std::cout << inod <<": "<< X << std::endl;
#endif
        nodes.push_back(X);
      }
    }
  };

  // Lambda function for parsing element connectivities from a string.
  auto&& parseElements = [&MNPC](const char* data, size_t nenod)
  {
    IntVec mnpc(nenod);
    std::istringstream iss(data);
#if SP_DEBUG > 1
    for (size_t iel = 0; iss; iel++)
#else
    while (iss)
#endif
    {
      for (int& n : mnpc) iss >> n;
      if (iss)
      {
#if SP_DEBUG > 1
        std::cout << iel <<":";
        for (int n : mnpc) std::cout <<" "<< n;
        std::cout << std::endl;
#endif
        MNPC.push_back(mnpc);
      }
    }
  };

  for (tag = tag->FirstChildElement(); tag; tag = tag->NextSiblingElement())
    if (tag->Value() && tag->FirstChild())
    {
      std::vector<ASM::NodeSet>* nset = nullptr;
      if (!strcasecmp(tag->Value(),"nodes"))
        parseNodes(tag->FirstChild()->Value());
      else if (!strcasecmp(tag->Value(),"elements"))
      {
        size_t nenod = 2;
        utl::getAttribute(tag,"nenod",nenod);
        parseElements(tag->FirstChild()->Value(),nenod);
      }
      else if (!strcasecmp(tag->Value(),"nodeset"))
        nset = &nodeSets;
      else if (!strcasecmp(tag->Value(),"elementset"))
        nset = elemSets;
      if (nset)
      {
        int node;
        IntVec nodes;
        std::string name;
        utl::getAttribute(tag,"name",name);
        std::istringstream iss(tag->FirstChild()->Value());
        iss >> node;
        while (iss)
        {
          nodes.push_back(1+node);
          iss >> node;
        }
#if SP_DEBUG > 1
        if (nset == &nodeSets)
          std::cout <<"Node ";
        else
          std::cout <<"Element ";
        std::cout <<"set \""<< name <<"\":";
        for (int n : nodes) std::cout <<" "<< n;
        std::cout << std::endl;
#endif
        nset->push_back(std::make_pair(name,nodes));
      }
    }

  return true;
}
