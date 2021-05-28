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
#include "Vec3Oper.h"
#include "Vec3.h"
#include <sstream>
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
