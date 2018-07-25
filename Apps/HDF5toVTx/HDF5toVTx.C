// $Id$
//==============================================================================
//!
//! \file HDF5toVTx.C
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
#include "ASM1D.h"
#include "ASM2D.h"
#include "ASM3D.h"
#include "ASMbase.h"
#include "ElementBlock.h"
#include "ProcessAdm.h"
#include "VTU.h"
#include <sstream>
#include <cstdlib>

bool isLR = false;

//! \brief Maps from basis name -> fields
typedef std::map< std::string,std::vector<XMLWriter::Entry> > ProcessList;

//! \brief Struct encapsulation information for a given output field
struct VTFFieldInfo {
  std::string Type;        //!< The type of the field (displacement, field)
  std::string Name;        //!< The name of the field
  std::vector<int> Blocks; //!< The block numbers in the file
};

//! \brief Maps from field description -> (type, VTF blocks)
typedef std::map<std::string,VTFFieldInfo> VTFList;

//! \brief Maps name<> designated ID block in VTF
typedef std::map<std::string,int> VTFFieldBlocks;

//! \brief Struct encapsulating information for a given basis
struct BasisInfo {
  std::vector<ASMbase*>   Patch;     //!< Vector with spline bases
  std::vector<RealArray*> FakeModel; //!< Vector with fake FE model evaluation points
  int                     StartPart; //!< Starting part for fake FE model
};

//! \brief Maps from (basis name -> basis info)
typedef std::map<std::string,BasisInfo> PatchMap;


//! \brief Read a basis from HDF5 into a vector of patch objects
//! \param result The resulting vector of patch objects
//! \param name The name of the basis
//! \param patches The number of patches to read
//! \param hdf The HDF5 reader to read from
//! \param dim The dimensionality of the basis
//! \param level The level the basis is stored at in the file
bool readBasis (std::vector<ASMbase*>& result, const std::string& name,
                int patches, HDF5Writer& hdf, int dim, int level)
{
  result.clear();
  if (dim < 1 || dim > 3)
  {
    std::cerr <<" *** readBasis: Invalid patch dimension "<< dim << std::endl;
    return false;
  }

  result.reserve(patches);
  ASM::Discretization ptype;
  for (int i=0;i<patches;++i) {
    std::stringstream geom, basis;
    geom << '/' << level << "/basis/";
    geom << name;
    geom << "/";
    geom << i+1;
    std::string out;
    hdf.readString(geom.str(),out);
    ptype = out.substr(0,10) == "# LRSPLINE" ? ASM::LRSpline : ASM::Spline;
    isLR = ptype == ASM::LRSpline;
    int gdim = 0;
    if (ptype == ASM::LRSpline)
      gdim = out.substr(11,7) == "SURFACE" ? 2 : 3;
    else {
      int sdim = atoi(out.substr(0,3).c_str());
      if (sdim == 700) sdim = 300;
      gdim = sdim / 100;
    }
    if (gdim != dim) {
      std::cerr << "Basis dimensionality for " << name
                << " does not match dimension, ignoring" << std::endl;
      continue;
    }
    basis << out;
    if (dim == 1)
      result.push_back(ASM1D::create(ptype));
    else if (dim == 2)
      result.push_back(ASM2D::create(ptype));
    else
      result.push_back(ASM3D::create(ptype));
    result.back()->read(basis);
    result.back()->generateFEMTopology();
  }

  return true;
}


//! \brief Write a field to VTF/VTU file
//! \param locvec The patch-level basis coefficients
//! \param components Number of components in field
//! \param patch The patch the field is defined on
//! \param model The evaluation points for this part of the fake FE model
//! \param geomID The ID associated with this patch
//! \param name Name of field
//! \param vlist List of vector fields stored in VTF/VTU
//! \param slist List of scalar fields stored in VTF/VTU
//! \param myVtf The VTF/VTU file to write to
bool writeFieldPatch(const Vector& locvec, int components,
                     ASMbase& patch, RealArray* model, int geomID, int& nBlock,
                     const std::string& name, VTFList& vlist, VTFList& slist,
                     VTF& myVtf, const std::string& description, const std::string& type)
{
  if (dynamic_cast<VTU*>(&myVtf) && type == "displacement")
    return true;
  Matrix field;
  if (!patch.evalSolution(field, locvec, model))
    return false;

  if (components > 1 || type == "eigenmodes") {
    if (!myVtf.writeVres(field,++nBlock,geomID,components))
      return false;
    else {
      vlist[description].Type = type;
      vlist[description].Name = name;
      vlist[description].Blocks.push_back(nBlock);
    }
    if (type == "displacement") {
      vlist["combined displacement"].Type = type;
      vlist["combined displacement"].Name = "combined displacements";
      vlist["combined displacement"].Blocks.push_back(nBlock);
    }
  }

  if (type == "eigenmodes")
    return true;

  for (size_t j = 0; j < field.rows(); j++) {
    std::string nam = name;
    if (field.rows() > 1) {
      nam += "_";
      nam += (char)('x'+j);
    }
    if (!myVtf.writeNres(field.getRow(1+j),++nBlock,geomID))
      return false;
    else {
      slist[nam].Name = nam;
      slist[nam].Blocks.push_back(nBlock);
    }
  }

  return true;
}


//! \brief Write a per-element field to VTF/VTU file
//! \param patch The patch the field is defined on
//! \param grid The (fake) finite elements associated with this patch
//! \param geomID The ID associated with this patch
//! \param nBlock Running VTF block counter
//! \param name Name of field
//! \param description Description of field
//! \param elist List of per-element fields stored in VTF/VTU
//! \param myVtf The VTF/VTU file to write to
bool writeElmPatch(const Vector& locvec,
                   ASMbase& patch, const ElementBlock* grid,
                   int geomID, int& nBlock, const std::string& description,
                   const std::string& name, VTFList& elist, VTF& myVtf)
{
  Matrix field(1,locvec.size());
  field.fill(locvec.data(),locvec.size());
  if (grid->getNoElms() > field.cols())
  {
    // Expand the element result array
    Matrix efield(field);
    field.resize(field.rows(),grid->getNoElms());
    for (size_t j = 1; j <= field.cols(); j++)
      field.fillColumn(j,efield.getColumn(grid->getElmId(j)));
  }

  if (!myVtf.writeEres(field.getRow(1),++nBlock,geomID))
    return false;
  else {
    elist[description].Type = "knotspan";
    elist[description].Name = name;
    elist[description].Blocks.push_back(nBlock);
  }

  return true;
}


//! \brief Write field blocks to VTF/VTU file
//! \param vlist List of vector fields
//! \param slist List of scalar fields
//! \param myvtf The VTF/VTU file to write to
//! \param iStep The level in file to write
void writeFieldBlocks(VTFList& vlist, VTFList& slist, VTF& myvtf,
                      int iStep, VTFFieldBlocks& fieldBlocks)
{
  static int idBlock = 20;
  for (VTFList::iterator it = vlist.begin(); it != vlist.end(); ++it) {
    if (fieldBlocks.find(it->second.Name) == fieldBlocks.end())
      fieldBlocks[it->second.Name] = idBlock++;
    if (it->second.Type == "displacement")
      myvtf.writeDblk(it->second.Blocks,it->second.Name.c_str(),fieldBlocks[it->second.Name],iStep);
    else if (it->second.Type == "eigenmodes")
      myvtf.writeDblk(it->second.Blocks, "Mode Shape", fieldBlocks[it->second.Name], iStep);
    else
      myvtf.writeVblk(it->second.Blocks,it->second.Name.c_str(),fieldBlocks[it->second.Name],iStep);
  }
  for (VTFList::iterator it = slist.begin(); it != slist.end(); ++it) {
    if (fieldBlocks.find(it->second.Name) == fieldBlocks.end())
      fieldBlocks[it->second.Name] = idBlock++;
    myvtf.writeSblk(it->second.Blocks,
                    it->second.Name.c_str(),fieldBlocks[it->second.Name],iStep);
  }
}


//! \brief Write a fake FE model to VTU/VTF file
//! \param patch The patch the fake FE model is associated with
//! \param id The ID associated with this patch
//! \param myVtf The VTF/VTU file to write to
//! \param nViz The number of visualization points per knot-span
//! \param block Running VTF block counter
void writePatchGeometry(ASMbase* patch, int id, VTF& myVtf, int* nViz, int block)
{
  std::stringstream str;
  str << "Patch " << id;

  size_t nd = patch->getNoParamDim();
  ElementBlock* lvb = new ElementBlock(nd == 3 ? 8 : (nd == 2 ? 4 : 2));
  patch->tesselate(*lvb,nViz);
  myVtf.writeGrid(lvb,str.str().c_str(),block);
}


//! \brief Generate evaluation points for fake FE model
//! \param result The resulting evaluation points
//! \param patches The patches to generate the fake FE model for
//! \param dims The dimensionality of the patches
//! \param nViz The number of visualization points per knot-span
void generateFEModel (std::vector<RealArray*>& result,
                      const std::vector<ASMbase*>& patches, int dims, int* nViz)
{
  result.clear();
  result.reserve(patches.size());
  for (size_t i=0;i<patches.size();++i) {
    RealArray* gpar = new RealArray[dims];
    for (int k=0;k<dims;++k) {
      if (dims == 1) {
        ASM1D* patch = dynamic_cast<ASM1D*>(patches[i]);
        if (patch) patch->getGridParameters(gpar[k],nViz[k]-1);
      }
      if (dims == 2) {
        ASM2D* patch = dynamic_cast<ASM2D*>(patches[i]);
        if (patch) patch->getGridParameters(gpar[k],k,nViz[k]-1);
      }
      else if (dims == 3) {
        ASM3D* patch = dynamic_cast<ASM3D*>(patches[i]);
        if (patch) patch->getGridParameters(gpar[k],k,nViz[k]-1);
      }
    }
    result.push_back(gpar);
  }
}


//! \brief Free up the memory allocated in the PatchMap, making sure to avoid double free's
//! \param map The map to free up
void freePatchMap(PatchMap& map)
{
  std::vector<RealArray*> free;
  for (PatchMap::iterator it = map.begin(); it != map.end(); ++it) {
    for (size_t i = 0; i < it->second.Patch.size(); ++i)
      delete it->second.Patch[i];
    free.insert(free.end(), it->second.FakeModel.begin(), it->second.FakeModel.end());
  }

  std::vector<RealArray*>::iterator it = std::unique(free.begin(), free.end());
  for (std::vector<RealArray*>::iterator it2 = free.begin(); it2 != it; ++it2)
    delete *it2;

  map.clear();
}


//! \brief Setup up for our basis.
//! \details Bases with the same size is collapsed to one, while different sized bases are appended
//! \param plist The process list with bases and fields
//! \param level The file level to load bases from
//! \param hdf The HDF5 file reader to use
//! \param dims The dimensionality of the basis
//! \param n The number of points in each direction per knot span in the tesselation
//! \param vtf The VTF/VTU file to write to
//! \param block Running VTF block counter
//! \param vtflevel The time level / load case in the VTF file
PatchMap setupPatchMap(const ProcessList& plist, int level, HDF5Writer& hdf, int dims, int* n,
                       VTF& vtf, int& block, int vtflevel)
{
  vtf.clearGeometryBlocks();
  int start=1;
  std::map<int, BasisInfo> created;
  PatchMap result;
  for (ProcessList::const_iterator it  = plist.begin();
                                   it != plist.end(); ++it) {
    if (it->first == "nodalforces")
      continue;
    readBasis(result[it->first].Patch, it->second[0].basis, it->second[0].patches, hdf, dims, level);
    std::map<int, BasisInfo>::iterator loc = created.find(it->second[0].patches);
    // tesselate some more patches
    if (loc == created.end()) {
      generateFEModel(created[it->second[0].patches].FakeModel, result[it->first].Patch, dims, n);
      loc = created.find(it->second[0].patches);
      loc->second.StartPart = start;
      for (size_t l=0;l<result[it->first].Patch.size();++l)
        writePatchGeometry(result[it->first].Patch[l], start+l, vtf, n, ++block);
      start += it->second[0].patches;
    }

    result[it->first].StartPart = loc->second.StartPart;
    result[it->first].FakeModel = loc->second.FakeModel;
  }

  return result;
}


int main (int argc, char** argv)
{
  int format = 1;
  int n[3] = { 2, 2, 2 };
  int dims = 3;
  int skip=1;
  int start=0;
  int end=-1;
  bool last=false;
  char* infile = 0;
  char* vtffile = 0;
  float starttime = -1, endtime = -1;

  for (int i = 1; i < argc; i++)
    if (!strcmp(argv[i],"-format") && i < argc-1) {
      if (!strcasecmp(argv[++i],"ascii"))
        format = 0;
      else if (!strcasecmp(argv[i],"binary"))
        format = 1;
      else
        format = atoi(argv[i]);
    }
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
    else if (!infile)
      infile = argv[i];
    else if (!vtffile)
      vtffile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile) {
    std::cout <<"usage: "<< argv[0]
              <<" <inputfile> [<vtffile>|<vtufile>] [-nviz <nviz>] \n"
              << "[-ndump <ndump>] [-last] [-start <level>] [-end <level>]\n"
              << "[-starttime <time>] [-endtime <time>] [-1D|-2D]\n"
              << "[-format <0|1|ASCII|BINARY>]\n";
    return 0;
  }
  else if (!vtffile)
    vtffile = infile;

  std::cout <<"\n >>> IFEM HDF5 to VT[F|U] converter <<<"
            <<"\n ==================================\n"
            <<"\nInput file: " << infile;

  std::cout <<"\nOutput file: "<< vtffile
            <<"\nNumber of visualization points: "
            << n[0] <<" "<< n[1] << " " << n[2] << std::endl;

  VTF* myVtf;
  if (strstr(vtffile,".vtf"))
    myVtf = new VTF(vtffile,format);
  else
    myVtf = new VTU(vtffile,last?1:0);

  // Process XML - establish fields and collapse bases
  PatchMap patches;

  HDF5Writer hdf(strtok(infile,"."),ProcessAdm(),true,true);
  XMLWriter xml(infile,ProcessAdm());
  xml.readInfo();

  int levels = xml.getLastTimeLevel();
  std::cout <<"Reading "<< infile <<": Time levels = "<< levels << std::endl;

  const std::vector<XMLWriter::Entry>& entry = xml.getEntries();
  std::vector<XMLWriter::Entry>::const_iterator it;

  ProcessList processlist;
  for (it = entry.begin(); it != entry.end(); ++it) {
    if (!it->basis.empty() && it->type != "restart") {
      processlist[it->basis].push_back(*it);
      std::cout << it->name <<"\t"<< it->description <<"\tnc="<< it->components
                <<"\t"<< it->basis << std::endl;
    }
    if (it->type == "eigenmodes") {
      levels = it->components-1;
      processlist[it->basis].back().components = 1;
    }
    if (it->type == "nodalforces")
      processlist["nodalforces"].push_back(*it);
  }

  if (processlist.empty()) {
    std::cout << "No fields to process, bailing" << std::endl;
    exit(1);
  }

  ProcessList::const_iterator pit = processlist.begin();

  double time = 0.0;

  // setup step boundaries and initial time
  if (starttime > 0)
    start = (int)(floor(starttime/pit->second.begin()->timestep));
  if (endtime > 0)
    end = int(endtime/pit->second.begin()->timestep+0.5f);
  if (end == -1)
    end = levels;
  time=last?end  *pit->second.begin()->timestep:
            start*pit->second.begin()->timestep;

  bool ok = true;
  int block = 0;
  VTFFieldBlocks fieldBlocks;
  int k = 1;
  for (int i = last?end:start; i <= end && ok; i += skip) {
    if (levels > 0) {
      if (processlist.begin()->second.begin()->timestep > 0) {
        hdf.readDouble(i,"timeinfo","SIMbase-1",time);
        std::cout <<"Time level "<< i;
        std::cout << " (t=" << time << ")";
      } else
        std::cout << "Step " << i+1;

      std::cout << std::endl;
    }
    VTFList vlist, slist;
    bool geomWritten=false;

    if ((isLR && hdf.hasGeometries(i)) || patches.empty()) {
      patches = setupPatchMap(processlist, hdf.hasGeometries(i)?i:0, hdf, dims, n, *myVtf, block, k);
      geomWritten = true;
    }

    for (pit = processlist.begin(); pit != processlist.end(); ++pit) {
      for (it = pit->second.begin(); it != pit->second.end() && ok; ++it) {
        if (it->once && k > 1)
          continue;
        if (pit->first != "nodalforces" && patches[pit->first].Patch.empty()) {
          if (k == 1)
            std::cerr << "Ignoring \"" << it->name << "\", basis not loaded" << std::endl;
          continue;
        }
        std::cout <<"Reading \""<< it->name <<"\""<< std::endl;
        if (pit->first == "nodalforces") {
          Vector vec;
          hdf.readVector(i, it->name, -1, vec);
          std::vector<Vec3Pair> pts(vec.size()/6);
          for (size_t j=0;j<vec.size()/6;++j) {
            for (int l=0;l<3;++l) {
              pts[j].first[l] = vec[6*j+l];
              pts[j].second[l] = vec[6*j+l+3];
            }
          }
          int geoBlck=-1;
          ok = myVtf->writeVectors(pts,geoBlck,++block,it->name.c_str(),k);
          continue;
        }
        for( int j=0;j<pit->second[0].patches;++j) {
          Vector vec;
          ok = hdf.readVector(it->once?0:i,it->name,j+1,vec);

          if (it->name.find('+') != std::string::npos) {
            /*
            Temporary hack to split a vector into scalar fields.
            The big assumption here is that the individual scalar names
            are separated by '+'-characters in the vector field name
            */
            Matrix tmp(it->components,vec.size()/it->components);
            tmp.fill(vec.ptr());
            size_t pos = 0;
            size_t fp = it->name.find('+');
            std::string prefix;
            if (fp != std::string::npos) {
              size_t fs = it->name.find(' ');
              if (fs < fp) {
                prefix = it->name.substr(0,fs+1);
                pos = fs+1;
              }
            }
            for (size_t r = 1; r <= tmp.rows() && pos < it->name.size(); r++) {
              size_t end = it->name.find('+',pos);

              ok &= writeFieldPatch(tmp.getRow(r),1, *patches[pit->first].Patch[j],
                                    patches[pit->first].FakeModel[j],
                                    patches[pit->first].StartPart+j,
                                    block, prefix+it->name.substr(pos,end-pos),
                                    vlist, slist, *myVtf, it->description, it->type);
              pos = end+1;
            }
          }
          else {
            if (it->type == "knotspan") {
              ok &= writeElmPatch(vec,*patches[pit->first].Patch[j],myVtf->getBlock(j+1),
                                  patches[pit->first].StartPart+j,block,
                                  it->description, it->name, slist, *myVtf);
            } else if (it->type == "eigenmodes") {
              ok &= writeFieldPatch(vec,it->components,
                                    *patches[pit->first].Patch[j],
                                    patches[pit->first].FakeModel[j],
                                    patches[pit->first].StartPart+j,
                                    block,it->name,vlist,slist,*myVtf, it->description, it->type);
            } else {
              ok &= writeFieldPatch(vec,it->components,
                                    *patches[pit->first].Patch[j],
                                    patches[pit->first].FakeModel[j],
                                    patches[pit->first].StartPart+j,
                                    block,it->name,vlist,slist,*myVtf, it->description, it->type);
            }
          }
        }
      }
    }
    if (geomWritten)
      myVtf->writeGeometryBlocks(k);
    writeFieldBlocks(vlist,slist,*myVtf,k,fieldBlocks);

    if (!ok)
      return 3;

    bool res;
    if (processlist.begin()->second.begin()->type == "eigenmodes") {
      double val;
      bool freq=false;
      if (!hdf.readDouble(i, "1", "eigenval", val)) {
        freq = true;
        hdf.readDouble(i, "1", "eigenfrequency", val);
      }
      res=myVtf->writeState(k++, freq?"Frequency %g" : "Eigenvalue %g", val, 1);
    } else if (processlist.begin()->second.begin()->timestep > 0)
      res=myVtf->writeState(k++,"Time %g",time,0);
    else {
      double foo = k;
      res=myVtf->writeState(k++,"Step %g", foo, 0);
    }
    if (!res) {
      std::cerr << "Error writing state" << std::endl;
      return 4;
    }

    pit = processlist.begin();
    time += pit->second.begin()->timestep*skip;
  }
  hdf.closeFile(levels,true);
  delete myVtf;

  return 0;
}
