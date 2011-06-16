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
  char* infile = 0;
  char* vtffile = 0;
  char* basis = 0;

  for (int i = 1; i < argc; i++)
    if (!strcmp(argv[i],"-vtf") && i < argc-1)
      format = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-nviz") && i < argc-1)
      n[0] = n[1] = n[2] = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-1D"))
      dims = 1;
    else if (!strcmp(argv[i],"-2D"))
      dims = 2;
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
              <<" <inputfile> [<vtffile>] [-nviz <nviz>] [-ndump <ndump>]"
              <<" [-basis <basis>] [-1D|-2D]"<< std::endl;
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
  VTF myVtf(vtffile,1);

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
    writePatchGeometry(gpatches[i],i+1,myVtf,n);
  FEmodel = generateFEModel(gpatches,dims,n);

  bool ok = true;
  int block = 0;
  double time=0;
  for (int i = 0; i <= levels && ok; i += skip) {
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
                                    slist,myVtf);
              pos = end+1;
            }
          }
          else {
            ok &= writeFieldPatch(vec,it->components,
                                  *patches[pit->first][j],
                                  FEmodel[j],j+1,
                                  block,it->name,vlist,slist,myVtf);
          }
        }
      }
    }
    writeFieldBlocks(vlist,slist,myVtf,i+1);

    if (ok)
      myVtf.writeState(i+1,"Time %g",time,0);
    else
      return 3;
    pit = processlist.begin();
    time += pit->second.begin()->timestep*skip;
  }
  hdf.closeFile(levels,true);

  return 0;
}
