// $Id$
//==============================================================================
//!
//! \file ASMu2DLag.C
//!
//! \date Aug 12 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of unstructured 2D Lagrange FE models.
//!
//==============================================================================

#include "ASMu2DLag.h"
#include "ElementBlock.h"
#include "Vec3Oper.h"
#include <numeric>
#include <sstream>
#include <cctype>


ASMu2DLag::ASMu2DLag (unsigned char n_s,
                      unsigned char n_f, char fType) : ASMs2DLag(n_s,n_f)
{
  fileType = fType;
}


ASMu2DLag::ASMu2DLag (const ASMu2DLag& p, unsigned char n_f) : ASMs2DLag(p,n_f)
{
  fileType = 0;
  nodeSets = p.nodeSets;
}


ASMu2DLag::ASMu2DLag (const ASMu2DLag& p) : ASMs2DLag(p)
{
  fileType = 0;
  nodeSets = p.nodeSets;
}


/*!
  \brief Helper method reading one data line ignoring all whitespace characters.
  \param is Input stream to read data from
  \param[out] s The read data is placed in this string
  \param[in] delim Read until this character is reached
*/

static bool readLine (std::istream& is, std::string& s, char delim = '\n')
{
  size_t i, n = 0;
  while (n == 0 && std::getline(is,s,delim))
    for (i = 0; i < s.length(); i++)
      if (!isspace(s[i])) s[n++] = s[i];

  s = s.substr(0,n);
  for (i = 0; i < n; i++)
    if (s[i] == ',') s[i] = ' ';

  return n > 0;
}


/*!
  This method reads a matlab file with arrays defining a standard FE mesh.
  The format corresponds to that of the output method SIMoutput::dumpMatlabGrid.
*/

bool ASMu2DLag::readMatlab (std::istream& is)
{
  std::string cline;
  if (!readLine(is,cline,'=') ||
      cline.find("function") == std::string::npos ||
      cline.find("Node")     == std::string::npos ||
      cline.find("Element")  == std::string::npos)
  {
    std::cerr <<" *** ASMu2DLag::readMatlab: Not a matlab file."<< std::endl;
    return false;
  }

  std::getline(is,cline);
  while (readLine(is,cline,'['))
    if (cline.find("Node=") == cline.length()-5)
    {
      while (readLine(is,cline))
      {
        size_t id; Vec3 X;
        std::stringstream(cline) >> id >> X;
#if SP_DEBUG > 1
        std::cout << id <<": "<< X << std::endl;
#endif
        this->setCoord(id,X);
        if (cline.back() == ';') break;
      }
#ifdef SP_DEBUG
      std::cout <<"Read "<< nnod <<" nodes."<< std::endl;
#endif
    }
    else if (cline.find("Element=") == cline.length()-8)
    {
      while (readLine(is,cline))
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
        else
        {
          std::cerr <<"  ** ASMu2DLag::readMatlab: "<< mnpc.size()
                    <<"-noded elements not supported (ignored).\n";
          continue;
        }
        if (id > myMNPC.size()) myMNPC.resize(id);
        myMNPC[id-1] = mnpc;
        if (cline.back() == ';') break;
      }
      nel = myMNPC.size();
#ifdef SP_DEBUG
      std::cout <<"Read "<< nel <<" elements."<< std::endl;
#endif
    }
    else if (cline.back() == '=')
    {
      std::string setname = cline.substr(0,cline.length()-1);
      if (readLine(is,cline,']'))
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


bool ASMu2DLag::read (std::istream& is)
{
  switch (fileType) {
  case 'm':
  case 'M':
    return this->readMatlab(is);
  default:
    std::cerr <<" *** ASMu2DLag::read: Undefined file format."<< std::endl;
    return false;
  }
}


bool ASMu2DLag::generateFEMTopology ()
{
  p1 = p2 = 2; // So far only linear elements supported

  myMLGN.resize(nnod);
  myMLGE.resize(nel);

  std::iota(myMLGN.begin(),myMLGN.end(),gNod+1);
  std::iota(myMLGE.begin(),myMLGE.end(),gEl+1);

  gNod += nnod;
  gEl  += nel;

  return true;
}


int ASMu2DLag::getNodeSetIdx (const std::string& setName) const
{
  int idx = 1;
  for (const auto& it : nodeSets)
    if (it.first == setName)
      return idx;
    else
      ++idx;

  return 0;
}


const IntVec& ASMu2DLag::getNodeSet (int idx) const
{
  int count = 0;
  for (const auto& it : nodeSets)
    if (++count == idx)
      return it.second;

  return this->ASMbase::getNodeSet(idx);
}


IntVec& ASMu2DLag::getNodeSet (const std::string& setName)
{
  for (NodeSet& it : nodeSets)
    if (it.first == setName)
      return it.second;

  nodeSets.push_back(std::make_pair(setName,IntVec()));
  return nodeSets.back().second;
}


void ASMu2DLag::generateThreadGroups (const Integrand&, bool, bool)
{
  // TODO: Add some coloring scheme later
  threadGroups.oneGroup(nel);
}


bool ASMu2DLag::tesselate (ElementBlock& grid, const int*) const
{
  grid.unStructResize(nel,nnod);

  size_t i, j, k;
  for (i = 0; i < nnod; i++)
    grid.setCoor(i,this->getCoord(1+i));

  for (i = k = 0; i < nel; i++)
    for (j = 0; j < MNPC[i].size(); j++)
      if (j > 1 && MNPC[i].size() == 4)
        grid.setNode(k++,MNPC[i][5-j]);
      else
        grid.setNode(k++,MNPC[i][j]);

  return true;
}
