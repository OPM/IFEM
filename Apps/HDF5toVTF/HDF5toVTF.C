// $Id$
//==============================================================================
//!
//! \file HDF5toVTF.C
//!
//! \date Apr 08 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Convert a HDF5 results database to VTF for visualisation.
//!
//==============================================================================

#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "HDF5Writer.h"
#include "XMLWriter.h"
#include "StringUtils.h"
#include <sstream>
#include <stdlib.h>
#include <string.h>

typedef std::map< std::string,std::vector<XMLWriter::Entry> > ProcessList;


int main (int argc, char** argv)
{
  int format = 0;
  int n[3] = { 5, 5, 5 };
  int dims = 3;
  char* infile = 0;
  char* vtffile = 0;

  for (int i = 1; i < argc; i++)
    if (!strcmp(argv[i],"-vtf") && i < argc-1)
      format = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-nviz") && i < argc-1)
      n[0] = n[1] = n[2] = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-1D"))
      dims = 1;
    else if (!strcmp(argv[i],"-2D"))
      dims = 2;
    else if (!infile)
      infile = argv[i];
    else if (!vtffile)
      vtffile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile) {
    std::cout <<"usage: "<< argv[0]
              <<" <inputfile> [<vtffile>] [-nviz <nviz>] [-1D|-2D]"<< std::endl;
    return 0;
  }
  else if (!vtffile)
    vtffile = infile;

  std::cout <<"\n >>> IFEM HDF5 to VTF converter <<<"
            <<"\n ==================================\n"
            <<"\nInput file: "<< infile
            <<"\nOutput file: "<< vtffile
            <<"\nNumber of visualization points: "
            << n[0] <<" "<< n[1] << " " << n[2] << std::endl;

  HDF5Writer hdf(strtok(infile,"."),true);
  XMLWriter xml(infile);
  xml.readInfo();

  int levels = xml.getLastTimeLevel();
  if (levels > 0) SIMinput::msgLevel = 1;
  std::cout <<"Reading "<< infile <<": Time levels = "<< levels << std::endl;

  const std::vector<XMLWriter::Entry>& entry = xml.getEntries();
  std::vector<XMLWriter::Entry>::const_iterator it;

  ProcessList processlist;
  for (it = entry.begin(); it != entry.end(); ++it)
    processlist[it->patchfile].push_back(*it);

  ProcessList::const_iterator pit = processlist.begin();
  for (int j = 1; pit != processlist.end(); ++pit, ++j)
  {
    SIMbase* sim;
    if (dims == 1)
      sim = new SIM1D;
    else if (dims == 2)
      sim = new SIM2D;
    else
      sim = new SIM3D;

    std::stringstream dummy;
    dummy << "PATCHFILE " << pit->first;
    sim->parse(const_cast<char*>(dummy.str().c_str()),dummy);
    std::stringstream dummy2;
    std::string foo(pit->first);
    dummy2 << "NODEFILE " << replaceAll(foo,".g2",".gno");
    sim->parse(const_cast<char*>(dummy2.str().c_str()),dummy2);

    if (!sim->preprocess(std::vector<int>(),false))
      return 1;

    bool ok = true;
    if (processlist.size() > 1) {
      std::string temp(strtok(vtffile,"."));
      std::stringstream str;
      str <<"-"<< j;
      temp.append(str.str());
      ok = sim->writeGlv(temp.c_str(),n,format);
    }
    else
      ok = sim->writeGlv(vtffile,n,format);

    int block = 0;
    for (int i = 0; i <= levels && ok; ++i) {
      int k = 20;
      if (levels > 0) std::cout <<"\nTime level "<< i << std::endl;
      for (it = pit->second.begin(); it != pit->second.end() && ok; ++it) {
        Vector vec;
        std::cout <<"Reading \""<< it->name <<"\""<< std::endl;
        ok = hdf.readField(i,it->name,vec,sim,it->components);
        if (it->description == "displacement")
          ok &= sim->writeGlvS(vec,n,i+1,block,i,1,NULL,10,it->components);
        else if (it->components < 2)
          ok &= sim->writeGlvS(vec,n,i+1,block,i,2,it->name.c_str(),k++,1);
        else
          ok &= sim->writeGlvV(vec,it->name.c_str(),n,i+1,block);
      }
      if (ok)
        ok = sim->writeGlvStep(i+1,i);
    }
    delete sim;

    if (!ok) return 2;
  }

  return 0;
}
