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
    if (!it->basis.empty())
    {
      processlist["PATCHFILE " + it->basis].push_back(*it);
      std::cout << it->name <<"\t"<< it->description <<"\tnc="<< it->components
		<<"\t"<< it->basis << std::endl;
    }

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

    if (!sim->parse(const_cast<char*>(pit->first.c_str()),std::cin))
      return 1;
    else if (!sim->preprocess(std::vector<int>(),false))
      return 2;

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
        if (it->name == "displacement") // displacement field
          ok &= sim->writeGlvS(vec,n,i+1,block,i,1,NULL,10,it->components);
        else if (it->components < 2) // scalar field
          ok &= sim->writeGlvS(vec,n,i+1,block,i,2,it->name.c_str(),k++,1);
        else if (it->name.find('+') < it->name.size())
	{
	  // Temporary hack to split a vector into scalar fields.
	  // The big assumption here is that the individual scalar names
	  // are separated by '+'-characters in the vector field name
	  Matrix tmp(it->components,vec.size()/it->components);
	  tmp.fill(vec.ptr());
	  size_t pos = 0;
	  for (size_t r = 1; r <= tmp.rows() && pos < it->name.size(); r++)
	  {
	    size_t end = it->name.find('+',pos);
	    ok &= sim->writeGlvS(tmp.getRow(r),n,i+1,block,i,2,
				 it->name.substr(pos,end).c_str(),k++,1);
	    pos = end+1;
	  }
	}
	else // just write as vector field
          ok &= sim->writeGlvV(vec,it->name.c_str(),n,i+1,block);
      }
      if (ok)
        ok = sim->writeGlvStep(i+1,i);
    }
    delete sim;

    if (!ok) return 3;
  }

  return 0;
}
