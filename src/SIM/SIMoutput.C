// $Id$
//==============================================================================
//!
//! \file SIMoutput.C
//!
//! \date Sep 05 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Sub-class with functionality for result output to VTF and terminal.
//!
//==============================================================================

#include "SIMoutput.h"
#include "SIMoptions.h"
#include "ASMbase.h"
#include "SAM.h"
#include "IntegrandBase.h"
#include "TensorFunction.h"
#include "AnaSol.h"
#include "Vec3Oper.h"
#include "VTF.h"
#include "ElementBlock.h"
#include "Functions.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml2.h"
#include <array>
#include <fstream>
#include <iomanip>
#include <numeric>


SIMoutput::SIMoutput (IntegrandBase* itg) : SIMinput(itg)
{
  myPrec = 3;
  myPtSize = 0.0;
  myGeomID = 0;
  myVtf = nullptr;
  logRpMap = false;
  idxGrid = -1;
}


SIMoutput::~SIMoutput ()
{
  if (myVtf) delete myVtf;

  for (std::pair<const std::string,RealFunc*>& func : myAddScalars)
    delete func.second;
}


void SIMoutput::clearProperties ()
{
  myPoints.clear();
  for (std::pair<const std::string,RealFunc*>& func : myAddScalars)
    delete func.second;
  myAddScalars.clear();
  this->SIMinput::clearProperties();
}


void SIMoutput::setPointResultFile (const std::string& filename, bool dumpCoord)
{
  if (filename.empty()) return;

  if (!myPoints.empty() && myPoints.back().first.empty())
    myPoints.back().first = filename; // points have already been defined
  else
    myPoints.push_back(std::make_pair(filename,ResPointVec()));
  std::string& myPtFile = myPoints.back().first;

  // Append extension .dat if nothing yet
  if (myPtFile.find_last_of('.') == std::string::npos)
    myPtFile.append(".dat");
  IFEM::cout <<"\tOutput file: "<< myPtFile << std::endl;

  // Append _coord to filename if coordinate dump is requested
  if (dumpCoord && myPtFile.find("_coord.") == std::string::npos)
    myPtFile.insert(myPtFile.find_last_of('.'),"_coord");

  // Append _p<PID> to filename unless on processor 0.
  // We use a non-empty myPatches to signal partitioning
  // as adm.dd.isPartitioned() is not yet initialized.
  if (nProc > 1 && !myPatches.empty())
  {
    char cPid[8];
    sprintf(cPid,"_p%04d",myPid);
    myPtFile.insert(myPtFile.find_last_of('.'),cPid);
  }
}


bool SIMoutput::parseOutputTag (const tinyxml2::XMLElement* elem)
{
  IFEM::cout <<"  Parsing <"<< elem->Value() <<">"<< std::endl;

  const char* funcval = utl::getValue(elem,"function");
  if (funcval)
  {
    std::string name;
    utl::getAttribute(elem,"name",name);
    if (name.empty()) return false;

    std::string type("expression");
    utl::getAttribute(elem,"type",type);
    IFEM::cout <<"\tAdditional function '"<< name <<"' ("<< type <<")";
    RealFunc* f = utl::parseRealFunc(funcval,type);
    if (!f) return false;
    IFEM::cout << std::endl;

    myAddScalars[name] = f;
    return true;
  }
  else if (strcasecmp(elem->Value(),"resultpoints"))
    return this->SIMinput::parseOutputTag(elem);

  // Parse the result point specifications.
  // Can either be explicit points, lines or a grid.
  bool newGroup = true;
  const tinyxml2::XMLElement* point = elem->FirstChildElement("point");
  for (int i = 1; point; i++, point = point->NextSiblingElement())
  {
    int patch = 0;
    ResultPoint thePoint;
    if (utl::getAttribute(point,"patch",patch) && patch > 0)
      thePoint.patch = patch;
    IFEM::cout <<"\tPoint "<< i <<": P"<< thePoint.patch;
    if (utl::getAttribute(point,"node",thePoint.inod) && thePoint.inod > 0)
      IFEM::cout <<" node = "<< thePoint.inod;
    else
      IFEM::cout <<" xi =";
    if (utl::getAttribute(point,"u",thePoint.u[0]))
      IFEM::cout <<' '<< thePoint.u[0];
    if (utl::getAttribute(point,"v",thePoint.u[1]))
      IFEM::cout <<' '<< thePoint.u[1];
    if (utl::getAttribute(point,"w",thePoint.u[2]))
      IFEM::cout <<' '<< thePoint.u[2];
    IFEM::cout << std::endl;
    if (newGroup)
      myPoints.push_back({"",{thePoint}});
    else
      myPoints.back().second.push_back(thePoint);
    newGroup = false;
  }

  const tinyxml2::XMLElement* line = elem->FirstChildElement("line");
  for (int j = 1; line; j++, line = line->NextSiblingElement())
  {
    int patch = 0;
    ResultPoint thePoint;
    if (utl::getAttribute(line,"patch",patch) && patch > 0)
      thePoint.patch = patch;

    double u0[3], u1[3];
    if (!utl::getAttribute(line,"u0",u0[0])) u0[0] = 0.0;
    if (!utl::getAttribute(line,"v0",u0[1])) u0[1] = 0.0;
    if (!utl::getAttribute(line,"w0",u0[2])) u0[2] = 0.0;
    if (!utl::getAttribute(line,"u1",u1[0])) u1[0] = u0[0];
    if (!utl::getAttribute(line,"v1",u1[1])) u1[1] = u0[1];
    if (!utl::getAttribute(line,"w1",u1[2])) u1[2] = u0[2];
    int npt = line->FirstChild() ? atoi(line->FirstChild()->Value()) : 2;
    if (u0[0] == u1[0] && u0[1] == u1[1] && u0[2] == u1[2]) npt = 1;

    memcpy(thePoint.u,u0,3*sizeof(double));
    if (newGroup)
      myPoints.push_back({"",{thePoint}});
    else
      myPoints.back().second.push_back(thePoint);
    newGroup = false;

    for (int i = 1; i < npt-1; i++)
    {
      double xi = double(i)/double(npt-1);
      for (int d = 0; d < 3; d++)
        thePoint.u[d] = u0[d]*(1.0-xi) + u1[d]*xi;
      myPoints.back().second.push_back(thePoint);
    }

    memcpy(thePoint.u,u1,3*sizeof(double));
    myPoints.back().second.push_back(thePoint);

    IFEM::cout <<"\tLine "<< j <<": P"<< thePoint.patch
               <<" npt = "<< npt <<" xi =";
    for (int d = 0; d < 3; d++)
    {
      IFEM::cout <<' '<< u0[d];
      if (u1[d] != u0[d])
        IFEM::cout <<'-'<< u1[d];
    }
    IFEM::cout << std::endl;
  }

  const tinyxml2::XMLElement* grid = elem->FirstChildElement("grid");
  if (!newGroup)
    grid = nullptr; // Don't mix grid output with other lines or points
  else if (grid)
    idxGrid = myPoints.size();

  for (int g = 1; grid; g++, grid = grid->NextSiblingElement())
  {
    int patch = 0;
    ResultPoint thePoint;
    if (utl::getAttribute(grid,"patch",patch) && patch > 0)
      thePoint.patch = patch;

    // Get grid bounds
    Vec3 X0, X1;
    bool center = false;
    bool cartesianGrid = (utl::getAttribute(grid,"X0",X0) &&
                          utl::getAttribute(grid,"X1",X1));
    double u0[3], u1[3];
    if (cartesianGrid)
      utl::getAttribute(grid,"center",center);
    else
    {
      if (!utl::getAttribute(grid,"u0",u0[0])) u0[0] = 0.0;
      if (!utl::getAttribute(grid,"v0",u0[1])) u0[1] = 0.0;
      if (!utl::getAttribute(grid,"w0",u0[2])) u0[2] = 0.0;
      if (!utl::getAttribute(grid,"u1",u1[0])) u1[0] = u0[0];
      if (!utl::getAttribute(grid,"v1",u1[1])) u1[1] = u0[1];
      if (!utl::getAttribute(grid,"w1",u1[2])) u1[2] = u0[2];
      X0 = Vec3(u0);
      X1 = Vec3(u1);
    }

    // Get grid resolution
    int npt[3] = { 10, 10, 10 };
    for (int d = 0; d < 3; d++)
      if (X0[d] == X1[d]) npt[d] = 1;
    if (grid->FirstChild())
    {
      char* sval = strdup(grid->FirstChild()->Value());
      char* cval = strtok(sval," ");
      for (int d = 0; d < 3 && cval; cval = strtok(nullptr," "))
      {
        while (d < 3 && npt[d] == 1) ++d;
        if (d < 3) npt[d++] = atoi(cval);
      }
      free(sval);
    }

    if (newGroup)
      myPoints.push_back({"",{}});
    newGroup = false;

    if (cartesianGrid)
    {
      thePoint.inod = -npt[0];
      Vec3 dX = X1 - X0;
      for (int d = 0; d < 3; d++)
        if (center && npt[d] > 1)
        {
          dX[d] /= npt[d];
          X0[d] += 0.5*dX[d];
        }
        else if (npt[d] > 2)
          dX[d] /= npt[d]-1;

      for (int k = 0; k < npt[2]; k++)
      {
        thePoint.X.z = X0.z + dX.z*k;
        for (int j = 0; j < npt[1]; j++)
        {
          thePoint.X.y = X0.y + dX.y*j;
          for (int i = 0; i < npt[0]; i++)
          {
            thePoint.X.x = X0.x + dX.x*i;
            myPoints.back().second.push_back(thePoint);
          }
        }
      }

      IFEM::cout <<"\tGrid "<< g <<": P"<< thePoint.patch
                 <<" npt = "<< npt[0]*npt[1]*npt[2] <<" X = "
                 << myPoints.back().second.front().X <<" - "
                 << myPoints.back().second.back().X;
    }
    else
    {
      for (int k = 0; k < npt[2]; k++)
      {
        double zeta = npt[2] > 1 ? double(k)/double(npt[2]-1) : 0.0;
        thePoint.u[2] = u0[2]*(1.0-zeta) + u1[2]*zeta;
        for (int j = 0; j < npt[1]; j++)
        {
          double eta = npt[1] > 1 ? double(j)/double(npt[1]-1) : 0.0;
          thePoint.u[1] = u0[1]*(1.0-eta) + u1[1]*eta;
          for (int i = 0; i < npt[0]; i++)
          {
            double xi = npt[0] > 1 ? double(i)/double(npt[0]-1) : 0.0;
            thePoint.u[0] = u0[0]*(1.0-xi) + u1[0]*xi;
            myPoints.back().second.push_back(thePoint);
          }
        }
      }

      IFEM::cout <<"\tGrid "<< g <<": P"<< thePoint.patch
                 <<" npt = "<< npt[0]*npt[1]*npt[2] <<" xi =";
      for (int d = 0; d < 3; d++)
      {
        IFEM::cout <<' '<< u0[d];
        if (u1[d] != u0[d])
          IFEM::cout <<'-'<< u1[d];
      }
    }
    IFEM::cout << std::endl;
  }

  utl::getAttribute(elem,"printmapping",logRpMap);
  utl::getAttribute(elem,"precision",myPrec);
  utl::getAttribute(elem,"vtfsize",myPtSize);

  std::string fname;
  if (utl::getAttribute(elem,"file",fname))
    this->setPointResultFile(fname,elem->FirstChildElement("dump_coordinates"));

  return true;
}


bool SIMoutput::parse (char* keyWord, std::istream& is)
{
  if (strncasecmp(keyWord,"RESULTPOINTS",12))
    return this->SIMinput::parse(keyWord,is);

  char* cline = keyWord+12;
  int nres = atoi(strtok(cline," "));
  if (nres > 0)
  {
    IFEM::cout <<"\nNumber of result points: "<< nres <<"\n";
    if ((cline = strtok(nullptr," ")))
      this->setPointResultFile(cline);
  }

  myPoints = {{"",ResPointVec(nres)}};
  for (int i = 0; i < nres && (cline = utl::readLine(is)); i++)
  {
    ResultPoint& thePoint = myPoints.back().second[i];
    thePoint.patch = atoi(strtok(cline," "));
    IFEM::cout <<"\tPoint "<< i+1 <<": P"<< thePoint.patch <<" xi =";
    for (int j = 0; j < 3 && (cline = strtok(nullptr," ")); j++)
    {
      thePoint.u[j] = atof(cline);
      IFEM::cout <<' '<< thePoint.u[j];
    }
    IFEM::cout << std::endl;
  }

  return true;
}


void SIMoutput::preprocessResultPoints ()
{
  for (ResPtPair& rptp : myPoints)
    this->preprocessResPtGroup(rptp.first,rptp.second);
}


void SIMoutput::preprocessResPtGroup (std::string& ptFile, ResPointVec& points)
{
  if (points.empty()) return;

  // Check if we have Cartesian grid input
  size_t nX = 0;
  IntMat pointMap;
  if (points.front().inod < 0 && points.back().inod < 0)
  {
    nX = -points.front().inod;
    size_t nY = points.size() / nX;
    if (nX*nY == points.size())
      pointMap.resize(nY,IntVec(nX,0));
    else
      nX = 0;
  }

  int iPoint = 0;
  size_t iX = 0, iY = 0;
  for (ResPointVec::iterator pit = points.begin(); pit != points.end();)
  {
    bool pointIsOK = false;
    ASMbase* pch = this->getPatch(pit->patch,true);
    if (pch && !pch->empty())
    {
      if ((pointIsOK = pit->inod > 0)) // A nodal number is specified
        pit->X = pch->getCoord(pit->inod);

      else if (pit->inod < 0) // A spatial point is specified
        pointIsOK = fabs(pch->findPoint(pit->X,pit->u)) < 1.0e-3;

      else // A parametric point is specified
        pointIsOK = (pit->inod = pch->evalPoint(pit->u,pit->u,pit->X)) >= 0;
    }

    if (pointIsOK)
    {
      (pit++)->npar = pch->getNoParamDim();
      ++iPoint;
      if (nX > 0)
        pointMap[iY][iX] = iPoint;
    }
    else
      pit = points.erase(pit);

    if (++iX == nX)
    {
      ++iY;
      iX = 0;
    }
  }
  if (iPoint == 0) return; // No valid result points

  if (logRpMap || msgLevel > 1)
  {
    iPoint = 0;
    for (const ResultPoint& pt : points)
    {
      if (++iPoint == 1) IFEM::cout <<'\n';
      IFEM::cout <<"Result point #"<< iPoint <<": patch #"<< pt.patch;
      switch (pt.npar) {
      case 1: IFEM::cout <<" u="; break;
      case 2: IFEM::cout <<" (u,v)=("; break;
      case 3: IFEM::cout <<" (u,v,w)=("; break;
      }
      IFEM::cout << pt.u[0];
      for (unsigned char c = 1; c < pt.npar; c++)
        IFEM::cout <<','<< pt.u[c];
      if (pt.npar > 1) IFEM::cout <<')';
      if (pt.inod > 0) IFEM::cout <<", node #"<< pt.inod;
      if (pt.inod > 0 && myModel.size() > 1)
      {
        ASMbase* pch = this->getPatch(pt.patch,true);
        IFEM::cout <<", global #"<< pch->getNodeID(pt.inod);
      }
      IFEM::cout <<", X = "<< pt.X << std::endl;
    }
  }

  if (!pointMap.empty())
  {
    // Write point mapping to file
    size_t idot = ptFile.find_last_of('.');
    std::string mapFile = ptFile.substr(0,idot) + ".map";
    std::ofstream os(mapFile,std::ios::out);
    for (const IntVec& row : pointMap)
    {
      for (int idx : row) os <<" "<< idx;
      os <<"\n";
    }
  }

  size_t icoord = ptFile.find("_coord.");
  if (icoord == std::string::npos) return;

  // Dump of result point coordinates to file is requested.
  // Formatted output, use scientific notation with fixed field width.
  std::streamsize precision = 3;
  std::streamsize flWidth = 8 + precision;
  std::ofstream os(ptFile,std::ios::out);
  os.flags(std::ios::scientific | std::ios::right);
  os.precision(precision);

  for (const ResultPoint& pt : points)
  {
    os << 0.0 <<" "; // dummy time
    for (unsigned char k = 0; k < pt.npar; k++)
      os << std::setw(flWidth) << pt.X[k];
    os << std::endl;
  }

  ptFile.erase(icoord,6); // erasing "_coord" from the file name
}


bool SIMoutput::merge (SIMbase* that, const std::map<int,int>* old2new,
                       int poff)
{
  if (!this->SIMbase::merge(that,old2new,poff))
    return false;
  else if (poff < 1)
    return true;

  int iPoint = 0;
  for (ResPtPair& rptp : static_cast<SIMoutput*>(that)->myPoints)
    for (ResultPoint& rp : rptp.second)
    {
      if (++iPoint == 1) IFEM::cout <<'\n';
      ASMbase* pch = that->getPatch(rp.patch,true);
      if (pch && rp.inod > 0)
        IFEM::cout <<"Result point #"<< iPoint <<": patch #"<< rp.patch
                   <<",  local node "<< rp.inod
                   <<"  global node "<< pch->getNodeID(rp.inod) << std::endl;
      rp.patch += poff;
    }

  return true;
}


ElementBlock* SIMoutput::tesselatePatch (size_t pidx) const
{
  if (pidx >= myModel.size())
    return nullptr;

  size_t nd = myModel[pidx]->getNoParamDim();
  ElementBlock* lvb = new ElementBlock(nd == 3 ? 8 : (nd == 2 ? 4 : 2));
  if (myModel[pidx]->tesselate(*lvb,opt.nViz))
    return lvb;

  delete lvb;
  return nullptr;
}


bool SIMoutput::writeGlvG (int& nBlock, const char* inpFile, bool doClear)
{
  if (adm.dd.isPartitioned() && adm.getProcId() != 0)
    return true; // write VTF-file only for procId=0 when domain decomposition

  // Lambda function for creating a new VTF-file name.
  auto&& getVTFname = [this](const char* inpFile)
  {
    char* vtfName = new char[strlen(inpFile)+10];
    strtok(strcpy(vtfName,inpFile),".");
    if (!adm.dd.isPartitioned() && nProc > 1)
      sprintf(vtfName+strlen(vtfName),"_p%04d",myPid);
#if HAS_VTFAPI == 2
    strcat(vtfName,".vtfx");
#else
    strcat(vtfName,".vtf");
#endif
    return vtfName;
  };

  if (inpFile)
  {
    if (myVtf)
    {
      std::cerr <<" *** SIMoutput::writeGlvG: Logic error,"
                <<" VTF file is already opened."<< std::endl;
      return false;
    }

    // Open a new VTF-file
    const char* fName = getVTFname(opt.vtf.empty() ? inpFile : opt.vtf.c_str());
    IFEM::cout <<"\nWriting VTF-file "<< fName << std::endl;
    myVtf = new VTF(fName,opt.format);
    delete[] fName;
  }
  else if (!myVtf)
  {
    std::cerr <<" *** SIMoutput::writeGlvG: Logic error,"
              <<" VTF file is not opened yet."<< std::endl;
    return false;
  }
  else if (doClear)
    myVtf->clearGeometryBlocks();

  ElementBlock* lvb;
  char pname[32];

  // Convert and write model geometry
  size_t pidx = 0;
  for (const ASMbase* pch : myModel)
  {
    ++pidx;
    if (pch->empty())
      continue; // skip empty patches

    if (!(lvb = this->tesselatePatch(pidx-1)))
      return false;

    if (msgLevel > 1)
      IFEM::cout <<"Writing geometry for patch "
                 << pch->idx+1 <<" ("<< lvb->getNoNodes() <<")"<< std::endl;

    sprintf(pname,"Patch %zu",pch->idx+1);
    if (!myVtf->writeGrid(lvb,pname,++nBlock))
      return false;
  }

  // Additional geometry for immersed boundaries
  for (const ASMbase* pch : myModel)
    if ((lvb = pch->immersedGeometry()))
    {
      sprintf(pname,"Immersed boundary %zu",pch->idx+1);
      if (!myVtf->writeGrid(lvb,pname,++nBlock))
        return false;
    }

  // Additional geometry for result point visualization
  if (myPtSize > 0.0 && !myPoints.empty())
  {
    lvb = new ElementBlock(8);
    for (const ResPtPair& points : myPoints)
      for (const ResultPoint& pt : points.second)
        lvb->merge(CubeBlock(pt.X,myPtSize));

    if (!myVtf->writeGrid(lvb,"Result points",++nBlock))
      return false;
  }

  // Additional geometry for the extra points
  if (myPtSize > 0.0 && !myTopPts.empty())
  {
    lvb = new ElementBlock(8);
    for (const IdxVec3& pt : myTopPts)
      lvb->merge(CubeBlock(pt.second,myPtSize));

    if (!myVtf->writeGrid(lvb,"Extra points",++nBlock))
      return false;
  }

  // Do not write the geometry blocks to file yet, method VTF::writeVectors()
  // might create an additional block for the point vectors
  return true;
}


bool SIMoutput::writeGlvBC (int& nBlock, int iStep) const
{
  if (!myVtf)
    return true;

  Matrix field;
  std::array<IntVec,6> dID;

  size_t j, n;
  int geomID = myGeomID;
  for (const ASMbase* pch : myModel)
  {
    if (pch->empty()) continue; // skip empty patches

    geomID++;
    size_t nbc = pch->getNoFields(1);
    int nNodes = pch->getNoNodes(-1);
    for (n = 2; n <= pch->getNoBasis(); n++)
      nNodes -= pch->getNoNodes(n);
    Matrix bc(nbc,nNodes);
    std::array<int,6> flag{0,0,0,0,0,0};
    ASMbase::BCVec::const_iterator bit;
    for (bit = pch->begin_BC(); bit != pch->end_BC(); ++bit)
      if ((n = pch->getNodeIndex(bit->node,true)) && n <= bc.cols())
      {
        if (!bit->CX && nbc > 0) bc(1,n) = flag[0] = 1;
        if (!bit->CY && nbc > 1) bc(2,n) = flag[1] = 1;
        if (!bit->CZ && nbc > 2) bc(3,n) = flag[2] = 1;
        if (!bit->RX && nbc > 3) bc(4,n) = flag[3] = 1;
        if (!bit->RY && nbc > 4) bc(5,n) = flag[4] = 1;
        if (!bit->RZ && nbc > 5) bc(6,n) = flag[5] = 1;
      }

    if (std::accumulate(flag.begin(),flag.end(),0) == 0)
      continue; // nothing on this patch

    if (msgLevel > 1)
      IFEM::cout <<"Writing boundary conditions for patch "
                 << pch->idx+1 << std::endl;

    if (!pch->evalSolution(field,bc,opt.nViz,nbc))
      return false;

    // The BC fields should either be 0.0 or 1.0
    if (opt.nViz[0] > 2 || opt.nViz[1] > 2 || opt.nViz[2] > 2)
      for (j = 1; j <= 6; j++)
        if (flag[j-1])
          for (n = 1; n <= field.cols(); n++)
            if (field(j,n) < 0.9999) field(j,n) = 0.0;

    for (j = 0; j < 6; j++)
      if (flag[j])
      {
        if (!myVtf->writeNres(field.getRow(1+j),++nBlock,geomID))
          return false;
        dID[j].push_back(nBlock);
      }
  }

  const char* label[6] = {
    "fix_X", "fix_Y", "fix_Z", "fix_RX", "fix_RY", "fix_RZ"
  };

  for (j = 0; j < 6; j++)
    if (!dID[j].empty())
      if (!myVtf->writeSblk(dID[j],label[j],1+j,iStep))
        return false;

  return true;
}


bool SIMoutput::writeGlvNo (int& nBlock, int iStep, int idBlock) const
{
  Vector elm(this->getNoElms(false));
  std::iota(elm.begin(),elm.end(),1.0);
  if (!this->writeGlvE(elm,iStep,nBlock,"Global element numbers",idBlock++))
    return false;

  Vector scl(this->getNoNodes());
  std::iota(scl.begin(),scl.end(),1.0);
  if (!this->writeGlvS(scl,"Global node numbers",iStep,nBlock,idBlock++))
    return false;

  scl.clear();
  for (int n : myLoc2Glb) scl.push_back(n);
  if (!this->writeGlvS(scl,"Original node numbers",iStep,nBlock,idBlock++))
    return false;
  else if (myDegenElm.empty())
    return true;

  scl.fill(0.0);
  scl.resize(this->getNoNodes(),0.0);
  for (const std::pair<const int,int>& degen : myDegenElm)
    scl[degen.first-1] = degen.second;
  return this->writeGlvS(scl,"Collapsed triangles",iStep,nBlock,idBlock);
}


bool SIMoutput::writeGlvT (int iStep, int& geoBlk, int& nBlock) const
{
  if (myVtf && myProblem->hasTractionValues())
  {
    if (msgLevel > 1)
      IFEM::cout <<"Writing boundary tractions"<< std::endl;
    return myProblem->writeGlvT(myVtf,iStep,geoBlk,nBlock);
  }

  return true;
}


bool SIMoutput::writeGlvV (const Vector& vec, const char* fieldName,
                           int iStep, int& nBlock, int idBlock, int ncmp) const
{
  if (vec.empty() || !myVtf)
    return true;

  Matrix field;
  Vector lovec;
  IntVec vID;

  bool empty = false;
  int geomID = myGeomID;
  for (const ASMbase* pch : myModel)
  {
    if (!this->extractNodeVec(vec,lovec,pch,ncmp,empty))
      return false;
    else if (pch->empty())
      continue; // skip empty patches

    if (msgLevel > 1)
      IFEM::cout <<"Writing vector field for patch "<< pch->idx+1 << std::endl;

    if (!pch->evalSolution(field,lovec,opt.nViz))
      return false;

    if (!myVtf->writeVres(field,++nBlock,++geomID,this->getNoSpaceDim()))
      return false;

    vID.push_back(nBlock);
  }

  return myVtf->writeVblk(vID,fieldName,idBlock,iStep);
}


/*!
  This method assumes the scalar field is attached to the first basis
  if we are using a mixed basis.
*/

bool SIMoutput::writeGlvS (const Vector& scl, const char* fieldName,
                           int iStep, int& nBlock, int idBlock) const
{
  if (scl.empty() || !myVtf)
    return true;

  const bool piolaMapping = myProblem ?
    myProblem->getIntegrandType() & Integrand::PIOLA_MAPPING : false;

  Matrix field;
  Vector lovec;
  IntVec sID;

  int basis  = 1;
  int geomID = myGeomID;
  for (const ASMbase* pch : myModel)
  {
    if (pch->empty()) continue; // skip empty patches

    if (msgLevel > 1)
      IFEM::cout <<"Writing scalar field \""<< fieldName
                 <<"\" for patch "<< pch->idx+1 << std::endl;

    int ncmp = scl.size() / this->getNoNodes(basis);
    this->extractPatchSolution(scl,lovec,pch,ncmp,basis);
    if (!pch->evalSolution(field,lovec,opt.nViz,ncmp,piolaMapping))
      return false;

    if (!myVtf->writeNres(field,++nBlock,++geomID))
      return false;

    sID.push_back(nBlock);
  }

  return myVtf->writeSblk(sID,fieldName,idBlock,iStep);
}


bool SIMoutput::writeGlvS (const Vector& psol, int iStep, int& nBlock,
                           double time, const char* pvecName,
                           int idBlock, int psolComps)
{
  idBlock = this->writeGlvS1(psol,iStep,nBlock,time,
                             pvecName,idBlock,psolComps);

  if (idBlock < 0)
    return false;
  else if (!this->writeAddFuncs(iStep,nBlock,50,time))
    return false;
  else if (idBlock == 0 || opt.pSolOnly)
    return true;

  return this->writeGlvS2(psol,iStep,nBlock,time,idBlock,psolComps);
}


/*!
  This method writes only the primary solution field to the VTF-file.
  If analytical solution fields are available, those fields are written as well.

  The way the solution is written depends on whether it is a scalar or vector
  field, and on the input parameters \a pvecName and \a psolComps, as follows:

  - If the primary solution is a vector field, it is written as a deformation
    (labelled "Solution") if \a pvecName is null or \a psolComps is negative,
    otherwise it is written as a named vector field (no deformation plot).
  - If the primary solution is a scalar field and \a pvecName is null,
    the field value is interpreted as a deformation along the global Z-axis.
  - If the primary solution is a vector field, each vector component is in
    addition written as a scalar field. If \a scalarOnly is \e true, it is only
    written as scalar field components (no deformation or vector field output).

  The scalar field components are labelled \a pvecName_(i)
  where (i) is in (x,y,z,rx,ry,rz), if \a pvecName is not null and
  \a psolComps is positive or the \ref myProblem member is null.
  Otherwise, their names are obtained by invoking the method
  IntegrandBase::getField1Name() of the \ref myProblem member,
  prefixed by \a pVecName if the latter is not null and \a psolComps < 0.
*/

int SIMoutput::writeGlvS1 (const Vector& psol, int iStep, int& nBlock,
                           double time, const char* pvecName,
                           int idBlock, int psolComps, bool scalarOnly)
{
  if (psol.empty() || !myVtf)
    return 0; // no primary solution

  const bool piolaMapping = myProblem ?
    myProblem->getIntegrandType() & Integrand::PIOLA_MAPPING : false;

  size_t nf     = scalarOnly ? 1 : this->getNoFields();
  size_t nVcomp = nf > this->getNoSpaceDim() ? this->getNoSpaceDim() : nf;
  bool haveXsol = false;
  if (mySol && psolComps >= 0)
  {
    if (nf == 1)
      haveXsol = mySol->getScalarSol() != nullptr;
    else
      haveXsol = mySol->getVectorSol() != nullptr;
  }

  std::array<IntVec,2> vID;
  std::vector<IntVec> sID, xID;
  sID.reserve(nf);
  if (haveXsol) xID.reserve(nf);

  Matrix field;
  Vector lovec;

  bool empty = false;
  int geomID = myGeomID;
  for (const ASMbase* pch : myModel)
  {
    if (!this->extractNodeVec(psol,lovec,pch,psolComps,empty))
      return false;
    else if (pch->empty())
      continue; // skip empty patches

    if (msgLevel > 1)
      IFEM::cout <<"Writing primary solution for patch "
                 << pch->idx+1 << std::endl;

    // Evaluate primary solution variables

    if (!pch->evalSolution(field,lovec,opt.nViz,0,piolaMapping))
      return -1;

    pch->filterResults(field,myVtf->getBlock(++geomID));

    if (!scalarOnly && (nVcomp > 1 || !pvecName))
    {
      // Output as vector field
      if (!myVtf->writeVres(field,++nBlock,geomID,nVcomp))
        return -2;

      vID[0].push_back(nBlock);
    }

    if (myProblem) // Compute application-specific primary solution quantities
      myProblem->primaryScalarFields(field);

    if (!this->writeScalarFields(field,geomID,nBlock,sID))
      return -3;

    if (haveXsol)
    {
      if (msgLevel > 1)
        IFEM::cout <<"Writing exact solution for patch "
                   << pch->idx+1 << std::endl;

      mySol->initPatch(pch->idx);

      // Evaluate exact primary solution

      const ElementBlock* grid = myVtf->getBlock(geomID);
      Vec3Vec::const_iterator cit = grid->begin_XYZ();
      field.fill(0.0);
      if (nf == 1) // Scalar solution
      {
        const RealFunc& pSol = *mySol->getScalarSol();
        for (size_t j = 1; cit != grid->end_XYZ(); j++, ++cit)
          field(1,j) = pSol(Vec4(*cit,time,grid->getParam(j-1)));
      }
      else
      {
        const VecFunc& pSol = *mySol->getVectorSol();
        for (size_t j = 1; cit != grid->end_XYZ(); j++, ++cit)
          field.fillColumn(j,pSol(Vec4(*cit,time,grid->getParam(j-1))).ptr());
        if (mySol->getScalarSol())
        {
          const RealFunc* psSol;
          size_t i = 0, r = pSol.dim();
          while ((psSol = mySol->getScalarSol(i++)) && r++ < field.rows())
          {
            cit = grid->begin_XYZ();
            const RealFunc& sSol = *psSol;
            for (size_t j = 1; cit != grid->end_XYZ(); j++, ++cit)
              field(r,j) = sSol(Vec4(*cit,time,grid->getParam(j-1)));
          }
        }
      }

      if (!scalarOnly && (nVcomp > 1 || !pvecName))
      {
        // Output as vector field
        if (!myVtf->writeVres(field,++nBlock,geomID,nVcomp))
          return -2;

        vID[1].push_back(nBlock);
      }

      if (!this->writeScalarFields(field,geomID,nBlock,xID))
        return -3;
    }
  }

  // Write result block identifications

  size_t i;
  bool ok = true;
  std::string pname(pvecName ? pvecName : "Solution");
  for (i = 0; i < 2 && ok; i++)
    if (!vID[i].empty())
    {
      std::string vname(i == 1 ? "Exact " + pname : pname);
      if (pvecName && psolComps >= 0)
        ok = myVtf->writeVblk(vID[i],vname.c_str(),idBlock+i,iStep);
      else
        ok = myVtf->writeDblk(vID[i],vname.c_str(),idBlock+i,iStep);
    }

  int nbc = this->getNoFields(1);
  if (idBlock <= nbc) idBlock = nbc+1; // since we might have written BCs above

  std::vector<std::string> xname;
  if (haveXsol) xname.reserve(nf);
  if (nf > 1) pname += "_w";
  for (i = 0; i < sID.size() && !sID[i].empty() && ok; i++)
  {
    if (myProblem && (!pvecName || psolComps <= 0))
    {
      pname = myProblem->getField1Name(i);
      if (pvecName && psolComps < 0)
        pname = std::string(pvecName) + " " + pname;
    }
    else if (i > 0 && nVcomp > 1 && i%nVcomp == 0)
    {
      pname.back() = 'r';
      pname += 'x';
    }
    else if (nf > 1)
      ++pname.back();
    ok = myVtf->writeSblk(sID[i],pname.c_str(),idBlock++,iStep);
    if (haveXsol) xname.push_back("Exact " + pname);
  }

  for (i = 0; i < xname.size() && i < xID.size() && !xID[i].empty() && ok; i++)
    ok = myVtf->writeSblk(xID[i],xname[i].c_str(),idBlock++,iStep);

  return ok ? idBlock : -4;
}


/*!
  This method writes only the secondary solution fields to the VTF-file.
  If analytical solution fields are available, those fields are written as well.
*/

bool SIMoutput::writeGlvS2 (const Vector& psol, int iStep, int& nBlock,
                            double time, int idBlock, int psolComps)
{
  if (psol.empty() || !myVtf)
    return true; // no primary solution
  else if (!myProblem)
    return false; // no integrand (should not happen at this point)
  else if (myProblem->getNoSolutions() < 1)
    return true; // no patch-level primary solution

  size_t nf = myProblem->getNoFields(2);
  if (nf < 1) return true; // no secondary solution

  bool haveAsol = false;
  if (mySol)
  {
    if (this->getNoFields() == 1)
      haveAsol = mySol->hasScalarSol() > 1;
    else
      haveAsol = mySol->hasVectorSol() > 1;
  }

  size_t nProj = (opt.discretization == ASM::Spline ||
                  opt.discretization == ASM::SplineC1) &&
    opt.project.find(SIMoptions::GLOBAL) != opt.project.end() ? nf : 0;

  std::array<IntVec,2> vID;
  std::vector<IntVec> sID;
  sID.reserve(nf+nProj);

  Matrix field, pdir;
  Vector lovec;

  size_t i, j;
  bool empty = false;
  int geomID = myGeomID;
  for (const ASMbase* pch : myModel)
  {
    if (!this->extractNodeVec(psol,myProblem->getSolution(),
                              pch,psolComps,empty))
      return false;
    else if (pch->empty())
      continue; // skip empty patches

    if (msgLevel > 1)
      IFEM::cout <<"Writing secondary solution for patch "
                 << pch->idx+1 << std::endl;

    // Direct evaluation of secondary solution variables

    myProblem->initResultPoints(time,true); // include principal stresses
    if (!this->initPatchForEvaluation(pch->idx+1))
      return false;

    if (!pch->evalSolution(field,*myProblem,opt.nViz))
      return false;

    size_t k = 0;
    pch->filterResults(field,myVtf->getBlock(++geomID));
    if (!this->writeScalarFields(field,geomID,nBlock,sID,&k))
      return false;

    // Write principal directions, if any, as vector fields

    size_t nPoints = field.cols();
    for (j = 0; j < 2 && myProblem->getPrincipalDir(pdir,nPoints,j+1); j++)
    {
      pch->filterResults(pdir,myVtf->getBlock(geomID));
      if (!myVtf->writeVres(pdir,++nBlock,geomID,this->getNoSpaceDim()))
        return -2;

      vID[j].push_back(nBlock);
    }

    if (nProj > 0)
    {
      // Projection of secondary solution variables (tensorial splines only)

      myProblem->initResultPoints(time);
      if (!pch->evalSolution(field,*myProblem,opt.nViz,'D'))
        return false;

      pch->filterResults(field,myVtf->getBlock(geomID));
      if (!this->writeScalarFields(field,geomID,nBlock,sID,&k))
        return false;
    }

    if (haveAsol)
    {
      if (msgLevel > 1)
        IFEM::cout <<"Writing exact solution for patch "
                   << pch->idx+1 << std::endl;

      mySol->initPatch(pch->idx);

      // Evaluate analytical solution variables

      const ElementBlock* grid = myVtf->getBlock(geomID);
      Vec3Vec::const_iterator cit = grid->begin_XYZ();
      for (j = 1; cit != grid->end_XYZ() && haveAsol; j++, ++cit)
      {
        Vec4 Xt(*cit,time,grid->getParam(j-1));
        if (mySol->hasScalarSol() == 3 || mySol->hasVectorSol() == 3)
          haveAsol = myProblem->evalSol(lovec,*mySol->getStressSol(),Xt);
        else if (this->getNoFields() == 1)
          haveAsol = myProblem->evalSol(lovec,*mySol->getScalarSecSol(),Xt);
        else
          haveAsol = myProblem->evalSol(lovec,*mySol->getVectorSecSol(),Xt);
        if (haveAsol && j == 1)
          field.resize(lovec.size(), field.cols());
        if (haveAsol)
          field.fillColumn(j,lovec);
      }

      if (haveAsol)
        if (!this->writeScalarFields(field,geomID,nBlock,sID,&k))
          return false;
    }
  }

  // Write result block identifications

  bool ok = true;
  std::string vname("Principal direction P1");
  for (i = 0; i < 2 && ok; i++, vname[vname.size()-1]++)
    if (!vID[i].empty())
      ok = myVtf->writeVblk(vID[i],vname.c_str(),idBlock+i,iStep);

  const char* prefix = haveAsol ? "FE" : nullptr;
  for (i = j = 0; i < nf && j < sID.size() && !sID[j].empty() && ok; i++)
    ok = myVtf->writeSblk(sID[j++],
                          myProblem->getField2Name(i,prefix).c_str(),
                          idBlock++,iStep);

  for (i = 0; i < nProj && j < sID.size() && !sID[j].empty() && ok; i++)
    ok = myVtf->writeSblk(sID[j++],
                          myProblem->getField2Name(i,"Projected").c_str(),
                          idBlock++,iStep);

  if (haveAsol)
    for (i = 0; j < sID.size() && !sID[j].empty() && ok; i++)
      ok = myVtf->writeSblk(sID[j++],
                            myProblem->getField2Name(i,"Exact").c_str(),
                            idBlock++,iStep);

  return ok;
}


bool SIMoutput::eval2ndSolution (const Vector& psol, double time, int psolComps)
{
  if (psol.empty() || !myVtf)
    return true; // no primary solution
  else if (!myProblem)
    return false; // no integrand (should not happen at this point)
  else if (myProblem->getNoSolutions() < 1)
    return true; // no patch-level primary solution
  else if (myProblem->getNoFields(2) < 1 || opt.pSolOnly)
    return true; // no secondary solution

  myProblem->initResultPoints(time);

  Matrix field;
  bool empty = false;
  for (const ASMbase* pch : myModel)
    if (!this->extractNodeVec(psol,myProblem->getSolution(),
                              pch,psolComps,empty))
      return false;
    else if (pch->empty())
      continue; // skip empty patches
    else if (!this->initPatchForEvaluation(pch->idx+1))
      return false;
    else if (!pch->evalSolution(field,*myProblem,opt.nViz))
      return false;

  return true;
}


/*!
  If \a iStep is zero (or negative), this method only evaluates the projected
  solution at the visualization points and updates the \a maxVal array,
  without writing data to the VTF file.
*/

bool SIMoutput::writeGlvP (const Vector& ssol, int iStep, int& nBlock,
                           int idBlock, const char* prefix,
                           std::vector<PointValues>* maxVal)
{
  if (ssol.empty() || !myVtf)
    return true; // no projected solution

  // Lambda function for updating (patch-wise) maximum result values.
  auto&& updateMaxVal = [](PointValues& maxVal, double res,
                           size_t pidx, const Vec3& Xpt)
  {
    if (maxVal.size() == 1)
      pidx = 0; // only compute global maximum values
    else if (pidx >= maxVal.size())
      return; // patch index out of range

    if (fabs(res) > fabs(maxVal[pidx].second))
      maxVal[pidx] = std::make_pair(Xpt,res);
  };

  std::vector<IntVec> sID;
  size_t nComp = myProblem->getNoFields(2);
  if (iStep > 0) sID.reserve(nComp);

  Matrix field;
  Vector lovec;

  Vector::const_iterator ssolIt = ssol.begin();
  int geomID = myGeomID;
  for (const ASMbase* pch : myModel)
  {
    if (pch->empty()) continue; // skip empty patches

    if (iStep > 0 && msgLevel > 1)
      IFEM::cout <<"Writing projected solution for patch "
                 << pch->idx+1 << std::endl;

    if (this->fieldProjections())
    {
      size_t nval = nComp*pch->getNoProjectionNodes();
      lovec = RealArray(ssolIt,ssolIt+nval);
      ssolIt += nval;
    }
    else
      this->extractPatchSolution(ssol,lovec,pch,nComp,1);

    // Evaluate the solution variables at the visualization points
    if (!pch->evalProjSolution(field,lovec,opt.nViz,nComp))
      return false;

    size_t j = 0;
    ++geomID;

    // Write out to VTF-file as scalar fields
    if (iStep > 0 && !this->writeScalarFields(field,geomID,nBlock,sID,&j))
      return false;
    else if (!maxVal)
      continue;

    // Update extremal values
    const ElementBlock* grid = myVtf->getBlock(geomID);
    for (j = 0; j < maxVal->size() && j < field.rows(); j++)
    {
      size_t indx = 0;
      double cmax = field.getRow(1+j).normInf(indx,1,true);
      if (indx > 0 && indx <= grid->getNoNodes())
        updateMaxVal((*maxVal)[j],cmax,pch->idx,grid->getCoord(indx-1));
    }
  }

  // Write result block identifications
  bool ok = true;
  for (size_t j = 0; j < sID.size() && !sID[j].empty() && ok; j++)
    ok = myVtf->writeSblk(sID[j],
                          myProblem->getField2Name(j,prefix).c_str(),
                          ++idBlock,iStep);

  return ok;
}


bool SIMoutput::writeScalarFields (const Matrix& field, int geomID,
                                   int& nBlock, std::vector<IntVec>& sID,
                                   size_t* i)
{
  size_t nS = 0;
  size_t& k = i ? *i : nS;
  for (size_t j = 1; j <= field.rows(); j++, k++)
    if (!myVtf->writeNres(field.getRow(j),++nBlock,geomID))
      return false;
    else if (k < sID.size())
      sID[k].push_back(nBlock);
    else
      sID.push_back({nBlock});

  return true;
}


bool SIMoutput::writeGlvF (const RealFunc& f, const char* fname,
                           int iStep, int& nBlock, int idBlock, double time)
{
  if (!myVtf)
    return true;

  IntVec sID;

  int geomID = myGeomID;
  for (const ASMbase* pch : myModel)
  {
    if (pch->empty()) continue; // skip empty patches

    if (msgLevel > 1)
      IFEM::cout <<"Writing function '"<< fname <<"' for patch "
                 << pch->idx+1 << std::endl;

    if (!myVtf->writeNfunc(f,time,++nBlock,++geomID))
      return false;

    sID.push_back(nBlock);
  }

  return myVtf->writeSblk(sID,fname,idBlock,iStep);
}


bool SIMoutput::writeGlvStep (int iStep, double value, int itype)
{
  if (myVtf)
    myVtf->writeGeometryBlocks(iStep);
  else
    return true;

  if (itype < 0)
    return true;
  else if (itype == 0)
    return myVtf->writeState(iStep,"Time %g",value,itype);
  else
    return myVtf->writeState(iStep,"Step %g",value,itype);
}


bool SIMoutput::writeGlvM (const Mode& mode, bool freq, int& nBlock)
{
  if (mode.eigVec.empty() || !myVtf)
    return true; // no eigen modes

  if (msgLevel > 1)
    IFEM::cout <<"Writing eigenvector for Mode "<< mode.eigNo << std::endl;

  Vector displ;
  Matrix field;
  IntVec vID;

  bool empty = false;
  int geomID = myGeomID;
  for (const ASMbase* pch : myModel)
  {
    if (!this->extractNodeVec(mode.eigVec,displ,pch,0,empty))
      return false;
    else if (pch->empty())
      continue; // skip empty patches

    if (myModel.size() > 1 && msgLevel > 1)
      IFEM::cout <<"."<< std::flush;

    if (!pch->evalSolution(field,displ,opt.nViz))
      return false;

    if (!myVtf->writeVres(field,++nBlock,++geomID))
      return false;

    vID.push_back(nBlock);
  }
  if (myModel.size() > 1 && msgLevel > 1)
    IFEM::cout << std::endl;

  int idBlock = 10;
  if (!myVtf->writeDblk(vID,"Mode Shape",idBlock,mode.eigNo))
    return false;

  return myVtf->writeState(mode.eigNo, freq ? "Frequency %g" : "Eigenvalue %g",
                           mode.eigVal, 1);
}


/*!
  If \a dualPrefix is not null, it is assumed that the \a norms were
  computed from the dual solution and therefore labelled as such.
  In addition to the \a norms, also the \a dualField function value is written
  as a scalar field to the VTF-file in that case, evaluated at the centre of
  each visualization element since this function is typically discontinuous.
*/

bool SIMoutput::writeGlvN (const Matrix& norms, int iStep, int& nBlock,
                           const std::vector<std::string>& prefix,
                           int idBlock, const char* dualPrefix)
{
  if (norms.empty() || !myVtf)
    return true; // no element norms

  NormBase* norm = myProblem->getNormIntegrand(mySol);

  // Lambda function telling whether a norm quantity should be saved or not.
  auto&& writeNorm = [norm,dualPrefix](size_t iGroup, size_t iNorm)
  {
    if (iGroup == 1 && dualPrefix)
      return iNorm <= 2; // hack: dual refinement indicators as norm #2
    else
      return norm->hasElementContributions(iGroup,iNorm);
  };

  size_t idxW = 0;
  size_t nrow = norms.rows();
  if (dualField && dualPrefix) idxW = ++nrow;
  std::vector<IntVec> sID(nrow);
  Matrix field;

  size_t i, j, k, l;
  int geomID = myGeomID;
  for (const ASMbase* pch : myModel)
  {
    if (pch->empty()) continue; // skip empty patches

    if (msgLevel > 1)
      IFEM::cout <<"Writing element norms for patch "
                 << pch->idx+1 << std::endl;

    const ElementBlock* grid = myVtf->getBlock(++geomID);
    pch->extractElmRes(norms,field);

    size_t ncol = grid->getNoElms();
    if (ncol > field.cols() || nrow > field.rows())
    {
      // Expand the element result array
      Matrix efield(field);
      field.resize(nrow,ncol);
      for (j = 1; j <= ncol; j++)
        field.fillColumn(j,efield.getColumn(grid->getElmId(j)));
      if (idxW && dualField->initPatch(pch->idx))
        for (j = 1; j <= ncol; j++)
          field(idxW,j) = dualField->getScalarValue(grid->getCenter(j));
    }

    for (k = 0, i = j = l = 1; i <= nrow; i++, l++)
    {
      if (l > norm->getNoFields(j))
        l = 1, ++j;

      if (writeNorm(j,l) || i == idxW)
      {
        if (!myVtf->writeEres(field.getRow(i),++nBlock,geomID))
          return false;

        sID[k++].push_back(nBlock);
      }
    }
  }

  std::string normName;
  for (k = 0, i = j = l = 1; k < sID.size() && !sID[k].empty(); i++, l++)
  {
    if (l > norm->getNoFields(j))
      l = 1, ++j;

    if (writeNorm(j,l) || i == idxW)
    {
      if (i == idxW && dualPrefix)
        normName = std::string(dualPrefix) + " extraction function";
      else if (j == 1 && dualPrefix)
        normName = norm->getName(j,l,dualPrefix);
      else if (j > 1 && j-2 < prefix.size())
        normName = norm->getName(j,l,prefix[j-2].c_str());
      else
        normName = norm->getName(j,l);

      if (!myVtf->writeSblk(sID[k++],normName.c_str(),++idBlock,iStep,true))
        return false;
    }
  }

  delete norm;
  return true;
}


bool SIMoutput::writeGlvE (const Vector& vec, int iStep, int& nBlock,
                           const char* name, int idBlock,
                           bool internalOrder) const
{
  if (!myVtf)
    return true;

  Vector field;
  IntVec sID;

  int geomID = myGeomID;
  for (const ASMbase* pch : myModel)
  {
    if (pch->empty()) continue; // skip empty patches

    if (msgLevel > 1)
      IFEM::cout <<"Writing element field '"<< name <<"' for patch "
                 << pch->idx+1 << std::endl;

    pch->extractElmRes(vec,field,internalOrder);
    if (!myVtf->writeEres(field,++nBlock,++geomID))
      return false;

    sID.push_back(nBlock);
  }

  return myVtf->writeSblk(sID,name,idBlock,iStep,true);
}


void SIMoutput::closeGlv ()
{
  delete myVtf;
  myVtf = nullptr;
}


bool SIMoutput::dumpMatlabGrid (std::ostream& os, const std::string& name,
                                const std::vector<std::string>& sets,
                                double scale) const
{
  // Write function definition
  os <<"function [Node, Element";
  for (const std::string& setname : sets) os <<", "<< setname;
  os <<"]="<< name << std::endl;

  // Write out nodes for the specified topology sets
  ASMbase* pch;
  TopologySet::const_iterator tit;
  for (const std::string& setname : sets)
  {
    std::set<int> nodeSet;
    if ((tit = myEntitys.find(setname)) != myEntitys.end())
      for (const TopItem& top : tit->second)
        if ((pch = this->getPatch(top.patch)))
          if (top.idim+1 == pch->getNoParamDim())
          {
            IntVec nodes;
            pch->getBoundaryNodes(top.item,nodes);
            for (int n : nodes) nodeSet.insert(n);
          }

    os <<"\n  "<< setname <<"=[";
    int count = 0;
    for (int node : nodeSet)
      if (++count == 1)
        os << node;
      else if (count%16 == 1)
        os <<",...\n    "<< node;
      else
        os <<", "<< node;
    os <<"];"<< std::endl;
  }

  // Write out all nodes
  os <<"\n  Node=[1";
  size_t nnod = this->getNoNodes();
  for (size_t inod = 1; inod <= nnod; inod++)
  {
    Vec3 X = this->getNodeCoord(inod);
    if (inod > 1) os <<"\n        "<< inod;
    for (size_t d = 0; d < nsd; d++)
      os <<", "<< scale*X[d];
  }
  os <<"];\n  Node=Node(:,2:"<< 1+nsd <<");"<< std::endl;

  // Write out all elements
  os <<"\n  Element=[1";
  size_t nen = 0, nel = this->getNoElms(true);
  for (size_t iel = 1; iel <= nel; iel++)
  {
    IntVec nodes;
    if (!this->getElmNodes(nodes,iel))
      return false;
    else if (nodes.size() == 4)
      std::swap(nodes[2],nodes[3]);
    if (iel > 1) os <<"\n           "<< iel;
    for (int n : nodes) os <<", "<< n;
    if (nen == 0)
      nen = nodes.size();
    else if (nen != nodes.size() || nen > 4)
    {
      std::cerr <<" *** SIMoutput::dumpMatlabGrid:"
                <<" Only linear standard elements are supported."<< std::endl;
      return false;
    }
  }
  os <<"];\n  Element=Element(:,2:"<< 1+nen <<");"<< std::endl;

  return true;
}


bool SIMoutput::dumpGeometry (std::ostream& os) const
{
  for (const ASMbase* pch : myModel)
    if (!pch->empty() && !pch->write(os))
      return false;

  return true;
}


void SIMoutput::dumpPrimSol (const Vector& psol, utl::LogStream& os,
                             bool withID) const
{
  if (psol.empty() || !mySam)
    return;

  for (const ASMbase* pch : myModel)
  {
    if (pch->empty()) continue; // skip empty patches

    Vector patchSol;
    pch->extractNodalVec(psol,patchSol,mySam->getMADOF());

    if (withID)
    {
      if (myModel.size() > 1)
        os <<"\n# Patch: "<< pch->idx+1;
      os <<"\n# inod/gnod\tNodal Coordinates\tSolution\n";
    }

    size_t ip = 0, nnod = pch->getNoNodes();
    for (size_t j = 1; j <= nnod; j++)
    {
      unsigned char nndof = pch->getNodalDOFs(j);
      if (nndof == 0)
        continue;
      else if (withID)
        os << j <<" "<< pch->getNodeID(j) <<"\t\t"<< pch->getCoord(j) <<"\t\t";

      os << utl::trunc(patchSol[ip++]);
      for (unsigned char k = 1; k < nndof; k++)
        os <<' '<< utl::trunc(patchSol[ip++]);
      os <<'\n';
    }
  }

  os.flush();
}


bool SIMoutput::dumpSolution (const Vector& psol, utl::LogStream& os) const
{
  if (psol.empty())
    return true;
  else if (!mySam)
    return false;

  Matrix field;

  for (const ASMbase* pch : myModel)
  {
    if (pch->empty()) continue; // skip empty patches

    if (myModel.size() > 1)
      os <<"\n# Patch: "<< pch->idx+1;

    // Extract and write primary solution
    size_t nf = pch->getNoFields(1);
    Vector& patchSol = myProblem->getSolution();
    pch->extractNodalVec(psol,patchSol,mySam->getMADOF());
    for (size_t k = 0; k < nf; k++)
    {
      os << myProblem->getField1Name(k,"# FE");
      for (size_t j = 1; j <= pch->getNoNodes(); j++)
      {
        std::pair<int,int> dofs = mySam->getNodeDOFs(j);
        int idof = dofs.first+k;
        if (idof <= dofs.second)
          os <<"\n"<< utl::trunc(patchSol[idof-1]);
      }
      os << std::endl;
    }

    // Evaluate secondary solution variables
    if (!this->initPatchForEvaluation(pch->idx+1))
      return false;
    if (!pch->evalSolution(field,*myProblem))
      return false;

    // Write the secondary solution
    for (size_t j = 1; j <= field.rows(); j++)
    {
      os << myProblem->getField2Name(j-1,"# FE");
      for (size_t k = 1; k <= field.cols(); k++)
        os <<"\n"<< utl::trunc(field(j,k));
      os << std::endl;
    }
  }

  return true;
}


bool SIMoutput::dumpResults (const Vector& psol, double time,
                             utl::LogStream& os, const ResPointVec& gPoints,
                             bool formatted, std::streamsize precision) const
{
  if (gPoints.empty())
    return true;

  if (adm.dd.isPartitioned() && adm.getProcId() != 0)
    return true; // dump only for procId=0 when domain decomposition

  const RealArray* reactionForces = this->getReactionForces();
  RealFunc*     psolScl = mySol ? mySol->getScalarSol() : nullptr;
  VecFunc*      psolVec = mySol ? mySol->getVectorSol() : nullptr;
  VecFunc*      ssolScl = mySol ? mySol->getScalarSecSol() : nullptr;
  TensorFunc*   ssolVec = mySol ? mySol->getVectorSecSol() : nullptr;
  STensorFunc*  ssolStr = mySol ? mySol->getStressSol() : nullptr;

  size_t nxsol = 0;
  if (psolScl)
    nxsol = 1;
  else if (psolVec)
    nxsol = this->getNoSpaceDim();

  for (const ASMbase* pch : myModel)
  {
    IntVec  points;
    Vec3Vec Xp;
    Matrix  sol1, sol2;
    if (!this->evalResults(psol,gPoints,pch,points,Xp,sol1,sol2))
      return false;
    else if (points.empty())
      continue; // no result points for this patch

    // Formatted output, use scientific notation with fixed field width
    std::streamsize flWidth = 8 + precision;
    std::streamsize oldPrec = os.precision(precision);
    std::ios::fmtflags oldF = os.flags(std::ios::scientific | std::ios::right);
    for (size_t j = 0; j < points.size(); j++)
    {
      if (!formatted)
        os << time <<" ";
      else if (points[j] < 0)
        os <<"  Point #"<< -points[j] <<":\tsol1 =";
      else
      {
        points[j] = pch->getNodeID(points[j]);
        os <<"  Node #"<< points[j] <<":\tsol1 =";
      }

      for (size_t i = 1; i <= sol1.rows(); i++)
        os << std::setw(flWidth) << utl::trunc(sol1(i,j+1));

      Vec3 sol4;
      if (psolScl)
        sol4.x = (*psolScl)(Vec4(Xp[j],time));
      else if (psolVec)
        sol4 = (*psolVec)(Vec4(Xp[j],time));
      if (formatted && nxsol > 0)
        os <<"\n\t\texact1";
      for (size_t k = 0; k < nxsol; k++)
        os << std::setw(flWidth) << utl::trunc(sol4[k]);

      if (opt.discretization >= ASM::Spline)
      {
        int isol = 1;
        for (size_t i = 1; i <= sol2.rows(); i++)
        {
          if (formatted && i%myProblem->getNo2ndSolPerLine() == 1)
            os <<"\n\t\tsol"<< ++isol <<" =";
          os << std::setw(flWidth) << utl::trunc(sol2(i,j+1));
        }

        Vector sol3;
        if (ssolScl)
          sol3 = (*ssolScl)(Vec4(Xp[j],time)).vec(this->getNoSpaceDim());
        else if (ssolVec)
          sol3 = (*ssolVec)(Vec4(Xp[j],time));
        else if (ssolStr)
          sol3 = (*ssolStr)(Vec4(Xp[j],time));
        if (formatted && !sol3.empty())
          os <<"\n\t\texact2";
        for (double s : sol3) os << std::setw(flWidth) << utl::trunc(s);
      }

      if (reactionForces && points[j] > 0)
      {
        Vector rf; // print nodal reaction forces for nodes with prescribed DOFs
        if (mySam->getNodalReactions(points[j],*reactionForces,rf))
        {
          if (formatted)
            os <<"\n\t\treac =";
          for (double f : rf) os << std::setw(flWidth) << utl::trunc(f);
        }
      }
      os << std::endl;
    }
    os.precision(oldPrec);
    os.flags(oldF);
  }

  return true;
}


bool SIMoutput::evalResults (const Vector& psol, const ResPointVec& gPoints,
                             const ASMbase* patch, IntVec& points, Vec3Vec& Xp,
                             Matrix& sol1, Matrix& sol2) const
{
  points.clear(); Xp.clear();
  if (gPoints.empty() || patch->empty())
    return true;

  // Find all evaluation points within this patch, if any
  std::array<RealArray,3> params;

  size_t jPoint = 0;
  for (const ResultPoint& pt : gPoints)
  {
    ++jPoint;
    if (pt.patch == patch->idx+1)
    {
      if (opt.discretization >= ASM::Spline)
      {
        points.push_back(pt.inod > 0 ? pt.inod : -jPoint);
        for (unsigned short int k = 0; k < pt.npar; k++)
          params[k].push_back(pt.u[k]);
        if (mySol) Xp.push_back(pt.X);
      }
      else if (pt.inod > 0)
      {
        points.push_back(pt.inod);
        if (mySol) Xp.push_back(pt.X);
      }
    }
  }

  if (points.empty())
    return true; // no points in this patch

  // Extract patch-level control/nodal point values of the primary solution
  patch->extractNodalVec(psol,myProblem->getSolution(),mySam->getMADOF(),-2);

  // Lambda function augmenting two matrices with the same number of rows.
  auto&& augment = [](Matrix& A, const Matrix& B)
  {
    if (A.empty())
      A = B;
    else if (!A.augmentCols(B))
    {
      std::cerr <<" *** SIMoutput::evalResults: Incompatible matrices, rows(A)="
                << A.rows() <<" rows(B)="<< B.rows() << std::endl;
      return false;
    }

    return true;
  };

  bool ok = true;
  Matrix tmp;
  if (opt.discretization < ASM::Spline)
    // Extract primary solution variables, for nodal points only
    ok = patch->getSolution(tmp,myProblem->getSolution(),points);
  else
    // Evaluate the primary solution variables
    ok = patch->evalSolution(tmp,myProblem->getSolution(),params.data(),false);

  if (!ok || !augment(sol1,tmp))
    return false;

  if (myProblem->getNoFields(2) < 1)
    return true; // no secondary solution variables

  // Initialize patch for secondary solution evaluation
  if (!this->initPatchForEvaluation(patch->idx+1))
    return false;

  // Evaluate the secondary solution variables
  if (!patch->evalSolution(tmp,*myProblem,params.data(),false))
    return false;

  return augment(sol2,tmp);
}


bool SIMoutput::initPatchForEvaluation (int patchNo) const
{
  int pindx = this->getLocalPatchIndex(patchNo);
  if (pindx < 1) return false;

  LocalSystem::patch = pindx-1;
  this->setPatchMaterial(pindx);
  return this->extractPatchDependencies(myProblem,myModel,pindx-1);
}


bool SIMoutput::dumpVector (const Vector& vsol, const char* fname,
                            utl::LogStream& os, std::streamsize precision) const
{
  if (vsol.empty() || myPoints.empty())
    return true;

  size_t ngNodes = this->getNoNodes();
  size_t nComp = vsol.size() / ngNodes;
  bool useMADOF = nComp*ngNodes != vsol.size();
  if (!useMADOF && nComp == this->getNoFields())
    nComp = 0; // using the number of primary field components

  Matrix sol1;
  Vector lsol;

  bool ok = true;
  std::ofstream* fs = nullptr;
  std::streamsize flWidth = 8 + precision;
  std::streamsize oldPrec;
  std::ios::fmtflags oldFlgs;

  int gCount = 0;
  for (const ResPtPair& rptp : myPoints)
  {
    if (gCount++ == idxGrid)
      continue; // Skip grid output

    if (fname)
    {
      // Formatted output, use scientific notation with fixed field width
      oldPrec = os.precision(precision);
      oldFlgs = os.flags(std::ios::scientific | std::ios::right);
    }
    else if (!rptp.first.empty())
    {
      // Output to a separate file for plotting
      fs = new std::ofstream(rptp.first,std::ios::out);
      oldPrec = fs->precision(precision);
      oldFlgs = fs->flags(std::ios::scientific | std::ios::right);
    }
    else
      continue; // skip output for this point group

    int iPoint = 1;
    for (const ASMbase* pch : myModel)
    {
      if (pch->empty()) continue; // skip empty patches

      // Find all evaluation points within this patch, if any
      std::array<RealArray,3> params;
      IntVec points;
      int jPoint = iPoint;
      for (const ResultPoint& pt : rptp.second)
      {
        if (pt.patch == pch->idx+1)
        {
          if (opt.discretization >= ASM::Spline)
          {
            points.push_back(pt.inod > 0 ? pt.inod : -jPoint);
            for (unsigned short int k = 0; k < pch->getNoParamDim(); k++)
              params[k].push_back(pt.u[k]);
          }
          else if (pt.inod > 0)
            points.push_back(pt.inod);
        }
        ++jPoint;
      }

      if (points.empty()) continue; // no points in this patch

      // Evaluate/extract nodal solution variables
      if (useMADOF)
        pch->extractNodalVec(vsol,lsol,mySam->getMADOF(),-2);
      else
        pch->extractNodeVec(vsol,lsol,nComp);
      if (opt.discretization >= ASM::Spline)
        ok &= pch->evalSolution(sol1,lsol,params.data(),false);
      else
        ok &= pch->ASMbase::getSolution(sol1,lsol,points);

      if (fs) // Single-line output to separate file
        for (size_t j = 0; j < points.size() && ok; j++, iPoint++)
        {
          *fs << iPoint;
          for (size_t i = 1; i <= sol1.rows(); i++)
            *fs << std::setw(flWidth) << utl::trunc(sol1(i,j+1));
          *fs << std::endl;
        }

      else // Formatted output to log stream
        for (size_t j = 0; j < points.size() && ok; j++, iPoint++)
        {
          if (points[j] < 0)
            os <<"  Point #"<< -points[j];
          else
          {
            points[j] = pch->getNodeID(points[j]);
            os <<"  Node #"<< points[j];
          }

          os <<":\t"<< fname <<" =";
          for (size_t i = 1; i <= sol1.rows(); i++)
            os << std::setw(flWidth) << utl::trunc(sol1(i,j+1));

          os << std::endl;
        }
    }

    if (fs)
      delete fs;
    else
    {
      os.precision(oldPrec);
      os.flags(oldFlgs);
    }
  }

  return ok;
}


bool SIMoutput::dumpResults (const Vector& psol, double time,
                             utl::LogStream& os,
                             bool formatted, std::streamsize precision) const
{
  if (psol.empty())
    return true;

  int gCount = 0;
  bool first = true;
  for (const ResPtPair& rptp : myPoints)
    if (gCount++ != idxGrid && (!formatted || rptp.first.empty()))
    {
      if (first)
        myProblem->initResultPoints(time);
      first = false;

      if (!this->dumpResults(psol,time,os,rptp.second,formatted,precision))
        return false;
    }

  return true;
}


bool SIMoutput::savePoints (const Vector& psol, double time, int step) const
{
  if (step < 1 || psol.empty())
    return true;

  int gCount = 0;
  bool first = true;
  for (const ResPtPair& rptp : myPoints)
    if (gCount++ != idxGrid && !rptp.first.empty())
      for (const ResultPoint& resPt : rptp.second)
        if (this->getLocalPatchIndex(resPt.patch) > 0)
        {
          if (first)
            myProblem->initResultPoints(time);
          first = false;

          std::ofstream fs(rptp.first,step > 1 ? std::ios::app : std::ios::out);
          utl::LogStream logs(fs);
          if (!this->dumpResults(psol,time,logs,rptp.second,false,myPrec))
            return false;
          break;
        }

  return true;
}


bool SIMoutput::saveResults (const Vector& psol, double time, int step) const
{
  if (idxGrid < 0 || step < 1 || psol.empty())
    return true;
  else if (idxGrid > static_cast<int>(myPoints.size()))
    return false;

  const std::string& filName = myPoints[idxGrid].first;
  const ResPointVec& gridPts = myPoints[idxGrid].second;
  if (filName.empty() || gridPts.empty())
    return false;

  size_t idot = filName.find_last_of('.');
  std::string fileName = filName.substr(0,idot);
  std::string fileExt  = filName.substr(idot);

  myProblem->initResultPoints(time);

  // Evaluate the solution variables in all grid points
  IntVec  points;
  Vec3Vec Xp;
  Matrix  sol, sol2;
  for (const ASMbase* pch : myModel)
    if (!this->evalResults(psol,gridPts,pch,points,Xp,sol,sol2))
      return false;

  if (!sol2.empty() && !sol.augmentRows(sol2))
    return false;

  // Write one file for each result component
  for (size_t r = 1; r <= sol.rows(); r++)
  {
    std::string fName = fileName + "_" + std::to_string(r) + fileExt;
    std::ofstream fs(fName,step > 1 ? std::ios::app : std::ios::out);
    fs << time;
    for (size_t c = 1; c <= sol.cols(); c++)
      fs <<" "<< sol(r,c);
    fs << '\n';
  }

  return true;
}


bool SIMoutput::extractNodeVec (const Vector& glbVec, Vector& locVec,
                                const ASMbase* patch, int nodalCmps,
                                bool& emptyPatches) const
{
  if (patch->empty())
  {
    emptyPatches = true;
    if (nodalCmps > 0 && nodalCmps != patch->getNoFields())
    {
      std::cerr <<" *** SIMoutput::extractNodeVec: Can't extract "<< nodalCmps
                <<"-nodal component vectors with empty patches."<< std::endl;
      return false; // can't use MADOF since different number of nodal DOFs
    }
  }
  else if (mySam->getNoNodes('X') > 0) // exclude any extraordinary DOFs
    patch->extractNodalVec(glbVec,locVec,mySam->getMADOF(),-2);
  else if (emptyPatches || this->getMDflag() > 0)
    patch->extractNodalVec(glbVec,locVec,mySam->getMADOF());
  else
  {
    if (nodalCmps < 0) nodalCmps = 0;
    patch->extractNodeVec(glbVec,locVec,nodalCmps);
  }

#if SP_DEBUG > 2
  std::cout <<"\nSolution vector for patch "<< patch->idx+1 << locVec;
#endif
  return true;
}


bool SIMoutput::hasPointResultFile () const
{
  for (const ResPtPair& rptp : myPoints)
    if (!rptp.first.empty())
      return true;

  return false;
}


bool SIMoutput::serialize (std::map<std::string,std::string>&) const
{
  std::cerr <<" *** SIMoutput::serialize: Must be implemented in sub-class.\n"
            <<"     Restart not supported for "<< this->getName() << std::endl;
  return false;
}


bool SIMoutput::writeAddFuncs (int iStep, int& nBlock, int idBlock, double time)
{
  for (const std::pair<const std::string,RealFunc*>& func : myAddScalars)
    if (!this->writeGlvF(*func.second, func.first.c_str(), iStep,
                         nBlock, idBlock++, time))
      return false;

  return true;
}


void SIMoutput::addAddFunc (const std::string& name, RealFunc* f)
{
  myAddScalars[name] = f;
}
