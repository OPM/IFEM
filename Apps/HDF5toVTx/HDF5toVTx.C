// $Id$
//==============================================================================
//!
//! \file HDF5toVTF.C
//!
//! \date Apr 08 2011
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Convert a HDF5 results database to VTF/VTU for visualisation.
//!
//==============================================================================

#include "HDF5Writer.h"
#include "XMLWriter.h"
#include "StringUtils.h"
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include "ASMs1D.h"
#include "ASMs2D.h"
#include "ASMs3D.h"
#include "ElementBlock.h"
#include "VTF.h"
#include "VTU.h"

typedef std::map< std::string,std::vector<XMLWriter::Entry> > ProcessList;
typedef std::map< std::string, std::vector<int> > VTFList;

std::vector<ASMbase*> readBasis(const std::string& name, 
                                int patches, HDF5Writer& hdf, int dim)
{
  std::vector<ASMbase*> result;
  for (int i=0;i<patches;++i) {
    std::stringstream geom;
    geom << "/0/basis/";
    geom << name;
    geom << "/";
    geom << i+1;
    std::string out;
    hdf.readString(geom.str(),out);
    std::stringstream basis;
    basis << out;
    if (dim == 1)
      result.push_back(new ASMs1D(basis,1,1));
    if (dim == 2)
      result.push_back(new ASMs2D(basis,2,1));
    if (dim == 3)
      result.push_back(new ASMs3D(basis,false,1));
    result.back()->generateFEMTopology();
  }

  return result;
}


bool writeFieldPatch(const Vector& locvec, int components,
                     ASMbase& patch, RealArray* model, int geomID, int& nBlock,
                     const std::string& name, VTFList& vlist, VTFList& slist,
                     VTF& myVtf)
{
  Matrix field;
  if (!patch.evalSolution(field,locvec,model))
    return false;

  if (components > 1) {
    if (!myVtf.writeVres(field,++nBlock,geomID,components))
      return false;
    else
      vlist[name].push_back(nBlock);
  }

  for (size_t j = 0; j < field.rows(); j++) {
    std::string nam = name;
    if (field.rows() > 1) {
      nam += "_";
      nam += (char)('x'+j);
    }
    if (!myVtf.writeNres(field.getRow(1+j),++nBlock,geomID))
      return false;
    else 
      slist[nam].push_back(nBlock);
  }

  return true;
}


void writeFieldBlocks(VTFList& vlist, VTFList& slist, VTF& myvtf,
                      int iStep)
{
  int idBlock = 20;
  for (VTFList::iterator it = vlist.begin(); it != vlist.end(); ++it) {
    myvtf.writeVblk(it->second,it->first.c_str(),
                    it->first=="displacement"?10:idBlock++,iStep);
  }
  for (VTFList::iterator it = slist.begin(); it != slist.end(); ++it) {
    myvtf.writeSblk(it->second,
                    it->first.c_str(),idBlock++,iStep);
  }
}


void writePatchGeometry(ASMbase* patch, int id, VTF& myVtf, int* nViz)
{
  std::stringstream str;
  str << "Patch " << id;

  size_t nd = patch->getNoParamDim();
  ElementBlock* lvb = new ElementBlock(nd == 3 ? 8 : (nd == 2 ? 4 : 2));
  patch->tesselate(*lvb,nViz);
  myVtf.writeGrid(lvb,str.str().c_str());
}


std::vector<RealArray*> generateFEModel(std::vector<ASMbase*> patches,
                                        int dims, int* n)
{
  std::vector<RealArray*> result;
  for (size_t i=0;i<patches.size();++i) {
    RealArray* gpar = new RealArray[dims];
    for (int k=0;k<dims;++k) {
      if (dims == 2) {
        ASMs2D* patch = (ASMs2D*)patches[i];
        patch->getGridParameters(gpar[k],k,n[k]-1);
      }
      if (dims == 3) {
        ASMs3D* patch = (ASMs3D*)patches[i];
        patch->getGridParameters(gpar[k],k,n[k]-1);
      }
    }
    result.push_back(gpar);
  }

  return result;
}


int main (int argc, char** argv)
{
  int format = 0;
  int n[3] = { 5, 5, 5 };
  int dims = 3;
  int skip=1;
  int start=0;
  int end=-1;
  bool last=false;
  char* infile = 0;
  char* vtffile = 0;
  char* basis = 0;
  float starttime = -1, endtime = -1;

  for (int i = 1; i < argc; i++)
    if (!strcmp(argv[i],"-vtf") && i < argc-1)
      format = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-nviz") && i < argc-1)
      n[0] = n[1] = n[2] = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-1D"))
      dims = 1;
    else if (!strcmp(argv[i],"-2D"))
      dims = 2;
    else if (!strcmp(argv[i],"-last"))
      last = true;
    else if (!strcmp(argv[i],"-start") && i < argc-1)
      start = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-starttime") && i < argc-1)
      starttime = atof(argv[++i]);
    else if (!strcmp(argv[i],"-end") && i < argc-1)
      end = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-endtime") && i < argc-1)
      endtime = atof(argv[++i]);
    else if (!strcmp(argv[i],"-ndump") && i < argc-1)
      skip = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-basis") && i < argc-1)
      basis = argv[i++];
    else if (!infile)
      infile = argv[i];
    else if (!vtffile)
      vtffile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile) {
    std::cout <<"usage: "<< argv[0]
              <<" <inputfile> [<vtffile>] [<vtufile>] [-nviz <nviz>]" << std::endl
              << "[-ndump <ndump>] [-last] [-start <level] [-end <level>]" << std::endl
              << "[-starttime <time>] [-endtime <time>] [-basis <basis>] [-1D|-2D]"<< std::endl;
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

  HDF5Writer hdf(strtok(infile,"."),true,true);
  XMLWriter xml(infile);
  xml.readInfo();

  int levels = xml.getLastTimeLevel();
  std::cout <<"Reading "<< infile <<": Time levels = "<< levels << std::endl;

  const std::vector<XMLWriter::Entry>& entry = xml.getEntries();
  std::vector<XMLWriter::Entry>::const_iterator it;

  ProcessList processlist;
  std::map<std::string, std::vector<ASMbase*> > patches;
  std::vector<RealArray*> FEmodel;
  VTF* myVtf;
  if (strstr(vtffile,".vtf"))
    myVtf = new VTF(vtffile,1);
  else
    myVtf = new VTU(vtffile,last?1:0);

  for (it = entry.begin(); it != entry.end(); ++it) {
    if (!it->basis.empty() && it->type != "restart") {
      processlist[it->basis].push_back(*it);
      std::cout << it->name <<"\t"<< it->description <<"\tnc="<< it->components
		<<"\t"<< it->basis << std::endl;
      if (patches[it->basis].empty())
        patches[it->basis] = readBasis(it->basis,it->patches,hdf,dims);
    }
  }

  ProcessList::const_iterator pit = processlist.begin();

  // This is broken with time dependent geometries.
  // Luckily it's not fundamentally broken so we can remedy when it's needed
  std::vector<ASMbase*> gpatches;
  if (basis) 
    gpatches = patches[basis];
  else
    gpatches = patches.begin()->second;
  for (int i=0;i<pit->second[0].patches;++i)
    writePatchGeometry(gpatches[i],i+1,*myVtf,n);
  FEmodel = generateFEModel(gpatches,dims,n);

  // setup step boundaries and initial time
  if (starttime > 0)
    start = (int)(floor(starttime/pit->second.begin()->timestep));
  if (endtime > 0)
    end = int(endtime/pit->second.begin()->timestep+0.5f);
  if (end == -1)
    end = levels;
  double time=last?end  *pit->second.begin()->timestep:
                   start*pit->second.begin()->timestep;

  bool ok = true;
  int block = 0;
  int k=1;
  for (int i = last?end:start; i <= end && ok; i += skip) {
    if (levels > 0) std::cout <<"\nTime level "<< i << " (t=" << time << ")" << std::endl;
    VTFList vlist, slist;
    for (pit = processlist.begin(); pit != processlist.end(); ++pit) {
      for (it = pit->second.begin(); it != pit->second.end() && ok; ++it) {
        std::cout <<"Reading \""<< it->name <<"\""<< std::endl;
        for( int j=0;j<pit->second[0].patches;++j) {
          Vector vec;
          ok = hdf.readVector(i,it->name,j+1,vec);

          if (it->name.find('+') != std::string::npos) {
            // Temporary hack to split a vector into scalar fields.
            // The big assumption here is that the individual scalar names
            // are separated by '+'-characters in the vector field name
            Matrix tmp(it->components,vec.size()/it->components);
            tmp.fill(vec.ptr());
            size_t pos = 0;
            for (size_t r = 1; r <= tmp.rows() && pos < it->name.size(); r++) {
              size_t end = it->name.find('+',pos);

              ok &= writeFieldPatch(tmp.getRow(r),1,*patches[pit->first][j],
                                    FEmodel[j],j+1,
                                    block,it->name.substr(pos,end),vlist,
                                    slist,*myVtf);
              pos = end+1;
            }
          }
          else {
            ok &= writeFieldPatch(vec,it->components,
                                  *patches[pit->first][j],
                                  FEmodel[j],j+1,
                                  block,it->name,vlist,slist,*myVtf);
          }
        }
      }
    }
    writeFieldBlocks(vlist,slist,*myVtf,k);

    if (ok)
      myVtf->writeState(k++,"Time %g",time,1);
    else
      return 3;
    pit = processlist.begin();
    time += pit->second.begin()->timestep*skip;
  }
  hdf.closeFile(levels,true);
  delete myVtf;

  return 0;
}
