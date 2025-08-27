// $Id$
//==============================================================================
//!
//! \file SIMmultiCpl.C
//!
//! \date Feb 12 2021
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Monolithic coupling of multiple simulators.
//!
//==============================================================================

#include "SIMmultiCpl.h"
#include "SIMoutput.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "IFEM.h"
#include "tinyxml2.h"


SIMmultiCpl::SIMmultiCpl (const std::vector<SIMoutput*>& sims)
  : SIMadmin(*sims.front()), mySims(sims)
{
  if (mySims.size() > 1)
    for (SIMoutput* sim : mySims)
      if (sim == mySims.front())
        sim->setMDflag(1);
      else if (sim == mySims.back())
        sim->setMDflag(2);
      else
        sim->setMDflag(3);
}


SIMmultiCpl::~SIMmultiCpl ()
{
  for (SIMoutput*& sim : mySims)
    delete sim;
}


bool SIMmultiCpl::parse (const tinyxml2::XMLElement* elem)
{
  bool results = true;
  if (!strcasecmp(elem->Value(),"coupling"))
  {
    const tinyxml2::XMLElement* child = elem->FirstChildElement("connection");
    for (; child; child = child->NextSiblingElement("connection"))
      results &= this->parseConnection(child);
  }
  else for (SIMbase* sim : mySims)
    if (!sim->parse(elem))
      return false;

  return results;
}


bool SIMmultiCpl::parseConnection (const tinyxml2::XMLElement* elem)
{
  IFEM::cout <<"  Parsing <"<< elem->Value() <<">"<< std::endl;

  std::string master, slave;
  utl::getAttribute(elem,"master",master);
  utl::getAttribute(elem,"slave",slave);
  SIMinput* mstSim = nullptr;
  SIMinput* slvSim = nullptr;
  for (SIMoutput* sim : mySims)
    if (sim->getEntity(master).size() == 1)
      mstSim = sim;
    else if (sim->getEntity(slave).size() == 1)
      slvSim = sim;

  if (!mstSim || !slvSim)
    return false;

  myCpl.push_back({mstSim,slvSim,
                   mstSim->getEntity(master).begin(),
                   slvSim->getEntity(slave).begin()});
  IFEM::cout <<"\tMaster point: \""<< master <<"\""<< *myCpl.back().master
             <<"\n\tSlave point:  \""<< slave <<"\""<< *myCpl.back().slave
             << std::endl;

  return true;
}


bool SIMmultiCpl::preprocessC (const IntVec& ignored, bool fixDup, double time0)
{
  // Preprocess the FE model of each sub-simulator
  IntVec empty;
  int nOffset = 0, pOffset = 0;
  std::map<SIMinput*,int> nSubNodes, nSubPatch;
  for (SIMoutput* sim : mySims)
    if (sim->preprocessC(nOffset == 0 ? ignored : empty, fixDup, time0))
    {
      nSubNodes[sim] = nOffset;
      nSubPatch[sim] = pOffset;
      nOffset += sim->getNoNodes();
      pOffset += sim->getNoPatches();
    }
    else
      return false;

  int substep = 10 + mySims.size();
  this->printHeading(substep);

  // Process the inter-sim couplings
  int misMatch = 0;
  std::map<int,int> cplNodes;
  for (const SIMcoupling& cpl : myCpl)
  {
    IntVec mNodes, sNodes;
    cpl.mstSim->getTopItemNodes(*cpl.master,mNodes);
    cpl.slvSim->getTopItemNodes(*cpl.slave,sNodes);
    for (int& n : mNodes) n += nSubNodes[cpl.mstSim];
    for (int& n : sNodes) n += nSubNodes[cpl.slvSim];
    if (sNodes.size() != mNodes.size())
    {
      std::cerr <<" *** SIMmultiCpl::preprocess: Mismatching topological items"
                <<" for inter-sim coupling "<< sNodes.size() <<" "
                << mNodes.size() << std::endl;
      return false;
    }
    else for (size_t i = 0; i < sNodes.size(); i++)
    {
      Vec3 Xm = cpl.mstSim->getNodeCoord(mNodes[i]-nSubNodes[cpl.mstSim]);
      Vec3 Xs = cpl.slvSim->getNodeCoord(sNodes[i]-nSubNodes[cpl.slvSim]);
      if (Xm.equal(Xs,1.0e-4))
        cplNodes[sNodes[i]] = mNodes[i];
      else
      {
        misMatch++;
        std::cerr <<" *** SIMmultiCpl::preprocess: Master node "<< mNodes[i]
                  <<" ("<< Xm <<") does not match slave node "<< sNodes[i]
                  <<" ("<< Xs <<")."<< std::endl;
      }
    }
  }
  if (misMatch > 0)
    return false;

  IFEM::cout <<"\nCoupling node mapping:";
  for (const std::pair<const int,int>& cp : cplNodes)
    IFEM::cout <<"\n\t"<< cp.first <<" -> "<< cp.second;
  IFEM::cout << std::endl;

  // Merge the equation systems into one monolithic system
  for (SIMoutput* sim : mySims)
    if (!mySims.front()->merge(sim,&cplNodes,nSubPatch[sim]))
      return false;

  return true;
}
