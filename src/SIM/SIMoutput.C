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
#include "tinyxml.h"
#include <fstream>
#include <iomanip>
#include <array>


SIMoutput::SIMoutput (IntegrandBase* itg) : SIMinput(itg)
{
  myPrec = 3;
  myGeomID = 0;
  myVtf = nullptr;
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
  this->SIMinput::clearProperties();
}


void SIMoutput::setPointResultFile (const std::string& filename, bool dumpCoord)
{
  if (filename.empty()) return;

  if (!myPoints.empty() && myPoints.back().first.empty())
    myPoints.back().first = filename; // Points have already been defined
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

  // Append _p<PID> to filename unless on processor 0
  if (nProc > 1)
  {
    char cPid[8];
    sprintf(cPid,"_p%04d",myPid);
    myPtFile.insert(myPtFile.find_last_of('.'),cPid);
  }
}


bool SIMoutput::parseOutputTag (const TiXmlElement* elem)
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

    myAddScalars[name] = f;
    return true;
  }
  else if (strcasecmp(elem->Value(),"resultpoints"))
    return this->SIMinput::parseOutputTag(elem);

  bool newGroup = true;
  const TiXmlElement* point = elem->FirstChildElement("point");
  for (int i = 1; point; i++, point = point->NextSiblingElement())
  {
    int patch = 0;
    ResultPoint thePoint;
    if (utl::getAttribute(point,"patch",patch) && patch > 0)
      thePoint.patch = patch;
    IFEM::cout <<"\tPoint "<< i <<": P"<< thePoint.patch <<" xi =";
    if (utl::getAttribute(point,"u",thePoint.u[0]))
      IFEM::cout <<' '<< thePoint.u[0];
    if (utl::getAttribute(point,"v",thePoint.u[1]))
      IFEM::cout <<' '<< thePoint.u[1];
    if (utl::getAttribute(point,"w",thePoint.u[2]))
      IFEM::cout <<' '<< thePoint.u[2];
    IFEM::cout << std::endl;
    if (newGroup)
      myPoints.push_back(std::make_pair("",ResPointVec(1,thePoint)));
    else
      myPoints.back().second.push_back(thePoint);
    newGroup = false;
  }

  const TiXmlElement* line = elem->FirstChildElement("line");
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
      myPoints.push_back(std::make_pair("",ResPointVec(1,thePoint)));
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

  utl::getAttribute(elem,"precision",myPrec);

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

  myPoints.resize(1);
  myPoints.back().second.resize(nres);
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
  for (ResPointVec::iterator p = points.begin(); p != points.end();)
  {
    ASMbase* pch = this->getPatch(p->patch,true);
    if (!pch || pch->empty())
      p = points.erase(p);
    else if ((p->inod = pch->evalPoint(p->u,p->u,p->X)) < 0)
      p = points.erase(p);
    else
      (p++)->npar = pch->getNoParamDim();
  }

  if (points.size() < 5 || msgLevel > 1)
  {
    int iPoint = 0;
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

  size_t icoord = ptFile.find("_coord.");
  if (icoord == std::string::npos) return;

  // Dump of result point coordinates to file is requested.
  // Formatted output, use scientific notation with fixed field width.
  std::streamsize precision = 3;
  std::streamsize flWidth = 8 + precision;
  std::ofstream os(ptFile.c_str(),std::ios::out);
  os.flags(std::ios::scientific | std::ios::right);
  os.precision(precision);

  for (const ResultPoint& pt : points)
  {
    os << 0.0 <<" "; // dummy time
    for (unsigned char k = 0; k < pt.npar; k++)
      os << std::setw(flWidth) << pt.X[k];
    os << std::endl;
  }

  ptFile.erase(icoord,6); // Erase "_coord" from the file name
}


bool SIMoutput::writeGlvG (int& nBlock, const char* inpFile, bool doClear)
{
  if (adm.dd.isPartitioned() && adm.getProcId() != 0)
    return true;

  // Lambda function for creating a new VTF-file name
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
    if (myVtf) return false;

    // Open a new VTF-file
    const char* fName = getVTFname(opt.vtf.empty() ? inpFile : opt.vtf.c_str());
    IFEM::cout <<"\nWriting VTF-file "<< fName << std::endl;
    myVtf = new VTF(fName,opt.format);
    delete[] fName;
  }
  else if (doClear && myVtf)
    myVtf->clearGeometryBlocks();
  else if (!myVtf)
    return false;

  ElementBlock* lvb;
  char pname[32];

  // Convert and write model geometry
  for (const ASMbase* pch : myModel)
  {
    if (pch->empty())
      continue; // skip empty patches

    if (msgLevel > 1)
      IFEM::cout <<"Writing geometry for patch "
                 << pch->idx+1 << std::endl;
    size_t nd = pch->getNoParamDim();
    lvb = new ElementBlock(nd == 3 ? 8 : (nd == 2 ? 4 : 2));
    if (!pch->tesselate(*lvb,opt.nViz))
      return false;

    sprintf(pname,"Patch %zu",pch->idx+1);
    if (!myVtf->writeGrid(lvb,pname,pch->idx+1,++nBlock))
      return false;
  }

  // Additional geometry for immersed boundaries
  for (const ASMbase* pch : myModel)
    if ((lvb = pch->immersedGeometry()))
    {
      sprintf(pname,"Immersed boundary %zu",pch->idx+1);
      if (!myVtf->writeGrid(lvb,pname,nGlPatches+pch->idx+1,++nBlock))
        return false;
    }

  // Do not write the geometry blocks to file yet, writeVectors might create
  // an additional block for the point vectors, see method VTF::writeVectors
  return true;
}


bool SIMoutput::writeGlvBC (int& nBlock, int iStep) const
{
  if (adm.dd.isPartitioned() && adm.getProcId() != 0)
    return true;
  else if (!myVtf)
    return false;

  Matrix field;
  std::array<IntVec,3> dID;

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
    RealArray flag(3,0.0);
    ASMbase::BCVec::const_iterator bit;
    for (bit = pch->begin_BC(); bit != pch->end_BC(); ++bit)
      if ((n = pch->getNodeIndex(bit->node,true)) && n <= bc.cols())
      {
        if (!bit->CX && nbc > 0) bc(1,n) = flag[0] = 1.0;
        if (!bit->CY && nbc > 1) bc(2,n) = flag[1] = 1.0;
        if (!bit->CZ && nbc > 2) bc(3,n) = flag[2] = 1.0;
      }

    if (flag[0]+flag[1]+flag[2] == 0.0)
      continue; // nothing on this patch

    if (msgLevel > 1)
      IFEM::cout <<"Writing boundary conditions for patch "
                 << pch->idx+1 << std::endl;

    if (!pch->evalSolution(field,bc,opt.nViz,nbc))
      return false;

    // The BC fields should either be 0.0 or 1.0
    if (opt.nViz[0] > 2 || opt.nViz[1] > 2 || opt.nViz[2] > 2)
      for (j = 1; j <= 3; j++)
        if (flag[j-1] == 1.0)
          for (n = 1; n <= field.cols(); n++)
            if (field(j,n) < 0.9999) field(j,n) = 0.0;

    for (j = 0; j < 3; j++)
      if (flag[j] == 1.0)
        if (myVtf->writeNres(field.getRow(1+j),++nBlock,geomID))
          dID[j].push_back(nBlock);
        else
          return false;
  }

  const char* label[3] = {
    "fix_X", "fix_Y", "fix_Z"
  };

  for (j = 0; j < 3; j++)
    if (!dID[j].empty())
      if (!myVtf->writeSblk(dID[j],label[j],1+j,iStep))
        return false;

  return true;
}


bool SIMoutput::writeGlvT (int iStep, int& geoBlk, int& nBlock) const
{
  if (adm.dd.isPartitioned() && adm.getProcId() != 0)
    return true;

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
  if (adm.dd.isPartitioned() && adm.getProcId() != 0)
    return true;
  else if (vec.empty())
    return true;
  else if (!myVtf)
    return false;

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
  if (adm.dd.isPartitioned() && adm.getProcId() != 0)
    return true;
  else if (scl.empty())
    return true;
  else if (!myVtf)
    return false;

  Matrix field;
  Vector lovec;
  IntVec sID;

  int basis  = 1;
  int geomID = myGeomID;
  for (const ASMbase* pch : myModel)
  {
    if (pch->empty()) continue; // skip empty patches

    if (msgLevel > 1)
      IFEM::cout <<"Writing scalar field for patch "<< pch->idx+1 << std::endl;

    int ncmp = scl.size() / this->getNoNodes(basis);
    this->extractPatchSolution(scl,lovec,pch,ncmp,basis);
    if (!pch->evalSolution(field,lovec,opt.nViz,ncmp))
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
  The primary solution is written as a deformation plot (labelled "Solution")
  if \a pvecName is null. If the primary solution is a scalar field, the field
  value is in that case interpreted as a deformation along the global Z-axis.
  If the primary solution is a vector field and \a pvecName is not null,
  it is written as a named vector field instead (no deformation plot),
  unless \a psolComps is negative.

  If the primary solution is a vector field, each vector component is written
  as a scalar field in addition. If \a scalarOnly is \e true, it is only
  written as scalar field components (no deformation or vector field output).

  If analytical solution fields are available, those fields are written as well.
*/

int SIMoutput::writeGlvS1 (const Vector& psol, int iStep, int& nBlock,
                           double time, const char* pvecName,
                           int idBlock, int psolComps, bool scalarOnly)
{
  if (adm.dd.isPartitioned() && adm.getProcId() != 0)
    return 0;
  else if (psol.empty())
    return 0;
  else if (!myVtf)
    return -99;

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

    if (!pch->evalSolution(field,lovec,opt.nViz))
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
          field.fillColumn(j,pSol(Vec4(*cit,time)).ptr());
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

      if (!this->writeScalarFields(field,geomID,nBlock,xID))
        return -3;

      if (!myVtf->writeVres(field,++nBlock,geomID,nVcomp))
        return -2;

      vID[1].push_back(nBlock);
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

  if (idBlock <= (int)this->getNoSpaceDim())
    idBlock = this->getNoSpaceDim()+1; // Since we might have written BCs above

  std::vector<std::string> xname;
  if (haveXsol) xname.reserve(nf);
  if (nf > 1) pname += "_w";
  for (i = 0; i < sID.size() && !sID[i].empty() && ok; i++)
  {
    if (myProblem && (!pvecName || nf > nVcomp))
      pname = myProblem->getField1Name(i);
    else if (nf > 1)
      (*pname.rbegin()) ++;
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
  if (adm.dd.isPartitioned() && adm.getProcId() != 0)
    return true;
  else if (psol.empty())
    return true; // No primary solution
  else if (!myVtf || !myProblem)
    return false;
  else if (myProblem->getNoSolutions() < 1)
    return true; // No patch-level primary solution

  size_t nf = myProblem->getNoFields(2);
  if (nf < 1) return true; // No secondary solution

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
  if (psol.empty())
    return true;
  else if (!myProblem)
    return false;
  else if (myProblem->getNoFields(2) < 1)
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
  \brief Static helper to update (patch-wise) maximum result values.
  \param maxVal Array of maximum values
  \param[in] res Result value at current point
  \param[in] pidx Zero-based patch index of current point
  \param[in] Xpt Coordinates of current point
*/

static bool updateMaxVal (PointValues& maxVal, double res,
                          size_t pidx, const Vec3& Xpt)
{
  if (maxVal.size() != 1 && pidx >= maxVal.size())
    return false; // Patch index out of range

  PointValue& pv = maxVal.size() > 1 ? maxVal[pidx] : maxVal.front();
  if (fabs(res) > fabs(pv.second)) pv = std::make_pair(Xpt,res);
  return true;
}


bool SIMoutput::writeGlvP (const Vector& ssol, int iStep, int& nBlock,
                           int idBlock, const char* prefix,
                           std::vector<PointValues>* maxVal)
{
  if (adm.dd.isPartitioned() && adm.getProcId() != 0)
    return true;
  else if (ssol.empty())
    return true;
  else if (!myVtf)
    return false;

  std::vector<IntVec> sID;
  size_t nComp = myProblem->getNoFields(2);
  sID.reserve(nComp);

  Matrix field;
  Vector lovec;

  Vector::const_iterator ssolIt = ssol.begin();
  int geomID = myGeomID;
  for (const ASMbase* pch : myModel)
  {
    if (pch->empty()) continue; // skip empty patches

    if (msgLevel > 1)
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

    size_t j = 0; // Write out to VTF-file as scalar fields
    const ElementBlock* grid = myVtf->getBlock(++geomID);
    if (!this->writeScalarFields(field,geomID,nBlock,sID,&j))
      return false;

    if (maxVal) // Update extremal values
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


bool SIMoutput::evalProjSolution (const Vector& ssol,
                                  std::vector<PointValues>& maxVal)
{
  if (ssol.empty())
    return true;
  else if (!myVtf)
    return true;

  Matrix field;
  Vector lovec;

  bool empty = false;
  int geomID = myGeomID;
  for (const ASMbase* pch : myModel)
  {
    if (!this->extractNodeVec(ssol,lovec,pch,myProblem->getNoFields(2),empty))
      return false;
    else if (pch->empty())
      continue; // skip empty patches

    // Evaluate the solution variables at the visualization points
    if (!pch->evalSolution(field,lovec,opt.nViz))
      return false;

    const ElementBlock* grid = myVtf->getBlock(++geomID);
    for (size_t j = 0; j < field.rows() && j < maxVal.size(); j++)
    {
      // Update extremal values
      size_t indx = 0;
      double cmax = field.getRow(1+j).normInf(indx,1,true);
      updateMaxVal(maxVal[j],cmax,pch->idx,grid->getCoord(indx-1));
    }
  }

  return true;
}


bool SIMoutput::writeGlvF (const RealFunc& f, const char* fname,
                           int iStep, int& nBlock, int idBlock, double time)
{
  if (adm.dd.isPartitioned() && adm.getProcId() != 0)
    return true;
  else if (!myVtf)
    return false;

  IntVec sID;

  int geomID = myGeomID;
  for (const ASMbase* pch : myModel)
  {
    if (pch->empty()) continue; // skip empty patches

    if (msgLevel > 1)
      IFEM::cout <<"Writing function '"<< fname
                 <<"' for patch "<< pch->idx+1 << std::endl;

    if (!myVtf->writeNfunc(f,time,++nBlock,++geomID))
      return false;

    sID.push_back(nBlock);
  }

  return myVtf->writeSblk(sID,fname,idBlock,iStep);
}


bool SIMoutput::writeGlvStep (int iStep, double value, int itype)
{
  if (adm.dd.isPartitioned() && adm.getProcId() != 0)
    return true;
  else if (!myVtf)
    return false;

  myVtf->writeGeometryBlocks(iStep);

  if (itype == 0)
    return myVtf->writeState(iStep,"Time %g",value,itype);
  else
    return myVtf->writeState(iStep,"Step %g",value,itype);
}


bool SIMoutput::writeGlvM (const Mode& mode, bool freq, int& nBlock)
{
  if (adm.dd.isPartitioned() && adm.getProcId() != 0)
    return true;
  else if (mode.eigVec.empty())
    return true;
  else if (!myVtf)
    return false;

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
  if (adm.dd.isPartitioned() && adm.getProcId() != 0)
    return true;
  else if (norms.empty())
    return true;
  else if (!myVtf)
    return false;

  NormBase* norm = myProblem->getNormIntegrand(mySol);

  // Lambda function telling whether a norm quantity should be saved or not
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
                           const char* name)
{
  if (adm.dd.isPartitioned() && adm.getProcId() != 0)
    return true;
  else if (!myVtf)
    return false;

  Matrix infield(1,vec.size());
  infield.fillRow(1,vec.ptr());
  Matrix field;
  IntVec sID;

  int geomID = myGeomID;
  for (const ASMbase* pch : myModel)
  {
    if (pch->empty()) continue; // skip empty patches

    if (msgLevel > 1)
      IFEM::cout <<"Writing element field '"<< name
                 <<"' for patch "<< pch->idx+1 << std::endl;

    pch->extractElmRes(infield,field);

    if (!myVtf->writeEres(field.getRow(1),++nBlock,++geomID))
      return false;

    sID.push_back(nBlock);
  }

  int idBlock = 300;
  if (!myVtf->writeSblk(sID,name,++idBlock,iStep,true))
    return false;

  return true;
}


void SIMoutput::closeGlv ()
{
  if (myVtf) delete myVtf;
  myVtf = nullptr;
}


bool SIMoutput::dumpMatlabGrid (std::ostream& os, const std::string& name,
                                const std::vector<std::string>& sets,
                                double scale) const
{
  if (this->getNoParamDim() != 2)
  {
    std::cerr <<" *** SIMoutput::dumpMatlabGrid: For 2D only."<< std::endl;
    return false;
  }

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
                <<" Only linear standard elements are supported.";
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

  Matrix  sol1, sol2;
  Vector  sol3, reactionFS;
  Vec3    sol4;
  IntVec  points;
  Vec3Vec Xp;

  const Vector* reactionForces = this->getReactionForces();
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

      if (psolScl)
        sol4.x = (*psolScl)(Vec4(Xp[j],time));
      else if (psolVec)
        sol4 = (*psolVec)(Vec4(Xp[j],time));
      if (nxsol > 0)
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

        if (ssolScl || ssolVec || ssolStr)
          os <<"\n\t\texact2";
        if (ssolScl)
          sol3 = (*ssolScl)(Vec4(Xp[j],time)).vec(this->getNoSpaceDim());
        else if (ssolVec)
          sol3 = (*ssolVec)(Vec4(Xp[j],time));
        else if (ssolStr)
          sol3 = (*ssolStr)(Vec4(Xp[j],time));
        for (double s : sol3)
          os << std::setw(flWidth) << utl::trunc(s);
      }

      if (reactionForces && points[j] > 0)
        // Print nodal reaction forces for nodes with prescribed DOFs
        if (mySam->getNodalReactions(points[j],*reactionForces,reactionFS))
        {
          if (formatted)
            os <<"\n\t\treac =";
          for (double f : reactionFS)
            os << std::setw(flWidth) << utl::trunc(f);
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

  int jPoint = 1;
  for (const ResultPoint& pt : gPoints)
  {
    if (pt.patch == patch->idx+1)
    {
      if (opt.discretization >= ASM::Spline)
      {
        points.push_back(pt.inod > 0 ? pt.inod : -jPoint);
        for (unsigned short int k = 0; k < patch->getNoParamDim(); k++)
          params[k].push_back(pt.u[k]);
        if (mySol) Xp.push_back(pt.X);
      }
      else if (pt.inod > 0)
      {
        points.push_back(pt.inod);
        if (mySol) Xp.push_back(pt.X);
      }
    }
    ++jPoint;
  }

  if (points.empty())
    return true; // no points in this patch

  // Extract patch-level control/nodal point values of the primary solution
  patch->extractNodeVec(psol,myProblem->getSolution());

  if (opt.discretization < ASM::Spline)
    // Extract primary solution variables, for nodal points only
    return patch->getSolution(sol1,myProblem->getSolution(),points);

  // Evaluate the primary solution variables
  if (!patch->evalSolution(sol1,myProblem->getSolution(),params.data(),false))
    return false;

  if (myProblem->getNoFields(2) < 1)
    return true; // No secondary solution variables

  // Initialize patch for secondary solution evaluation
  if (!this->initPatchForEvaluation(patch->idx+1))
    return false;

  // Evaluate the secondary solution variables
  return patch->evalSolution(sol2,*myProblem,params.data(),false);
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
  if (nComp*ngNodes != vsol.size() || nComp == this->getNoFields())
    nComp = 0; // Using the number of primary field components

  Matrix sol1;
  Vector lsol;

  bool ok = true;
  std::ofstream* fs = nullptr;
  std::streamsize flWidth = 8 + precision;
  std::streamsize oldPrec;
  std::ios::fmtflags oldFlgs;

  for (const ResPtPair& rptp : myPoints)
  {
    if (fname)
    {
      // Formatted output, use scientific notation with fixed field width
      oldPrec = os.precision(precision);
      oldFlgs = os.flags(std::ios::scientific | std::ios::right);
    }
    else if (!rptp.first.empty())
    {
      // Output to a separate file for plotting
      fs = new std::ofstream(rptp.first.c_str(),std::ios::out);
      oldPrec = fs->precision(precision);
      oldFlgs = fs->flags(std::ios::scientific | std::ios::right);
    }
    else
      continue; // Skip output for this point group

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

  myProblem->initResultPoints(time);

  for (const ResPtPair& rptp : myPoints)
    if (!formatted || rptp.first.empty())
      if (!this->dumpResults(psol,time,os,rptp.second,formatted,precision))
        return false;

  return true;
}


bool SIMoutput::savePoints (const Vector& psol, double time, int step) const
{
  if (step < 1 || psol.empty())
    return true;

  myProblem->initResultPoints(time);

  for (const ResPtPair& rptp : myPoints)
    if (!rptp.first.empty())
      for (const ResultPoint& resPt : rptp.second)
        if (this->getLocalPatchIndex(resPt.patch) > 0)
        {
          std::ofstream fs(rptp.first.c_str(),
                           step == 1 ? std::ios::out : std::ios::app);
          utl::LogStream logs(fs);
          if (!this->dumpResults(psol,time,logs,rptp.second,false,myPrec))
            return false;
          break;
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
      return false; // Can't use MADOF since different number of nodal DOFs
    }
  }
  else if (emptyPatches)
    patch->extractNodalVec(glbVec,locVec,mySam->getMADOF());
  else
  {
    if (nodalCmps < 0) nodalCmps = 0;
    patch->extractNodeVec(glbVec,locVec,nodalCmps);
  }

  return true;
}


bool SIMoutput::hasPointResultFile () const
{
  for (const ResPtPair& rptp : myPoints)
    if (!rptp.first.empty())
      return true;

  return false;
}


size_t SIMoutput::getVCPindex (size_t idx) const
{
  if (extrFunc.empty() || idx == 0 || idx > extrFunc.size())
    return 0;

  return (this->haveAnaSol() ? 4 : 2) + idx;
}


void SIMoutput::printNorms (const Vectors& norms, size_t w) const
{
  if (norms.empty()) return;

  NormBase* norm = this->getNormIntegrand();
  const Vector& n = norms.front();

  IFEM::cout <<"Energy norm"
             << utl::adjustRight(w-11,norm->getName(1,1)) << n(1);
  if (n(2) != 0.0)
    IFEM::cout <<"\nExternal energy"
               << utl::adjustRight(w-15,norm->getName(1,2)) << n(2);

  if (this->haveAnaSol() && n.size() >= 4)
    IFEM::cout <<"\nExact norm"
               << utl::adjustRight(w-10,norm->getName(1,3)) << n(3)
               <<"\nExact error"
               << utl::adjustRight(w-11,norm->getName(1,4)) << n(4)
               <<"\nExact relative error (%) : "<< 100.0*n(4)/n(3);

  size_t i, k = 0;
  while ((i = this->getVCPindex(++k)) && i <= n.size())
    IFEM::cout <<"\nVCP quantity"
               << utl::adjustRight(w-12,norm->getName(1,i)) << n(i);

  size_t j = 0;
  for (const SIMoptions::ProjectionMap::value_type& prj : opt.project)
    if (++j < norms.size())
      this->printNormGroup(norms[j],n,prj.second);

  IFEM::cout << std::endl;
  delete norm;
}


double SIMoutput::getReferenceNorm (const Vectors& gNorm, size_t adaptor) const
{
  if (gNorm.empty() || gNorm.front().empty())
    return 0.0;

  const Vector& fNorm = gNorm.front();
  if (adaptor < 1 && fNorm.size() > 2)
    return fNorm(3); // Using the analytical solution, |u|_ref = |u|
  else if (adaptor >= gNorm.size() || gNorm[adaptor].size() < 2)
    return -(double)adaptor; // Norm group index is out of range

  // |u|_ref = sqrt( |u^h|^2 + |e^*|^2 )
  return hypot(fNorm(1),gNorm[adaptor](2));
}


double SIMoutput::getEffectivityIndex (const Vectors& gNorm,
                                       size_t idx, size_t inorm) const
{
  return gNorm[idx](inorm) / gNorm.front()(4);
}


bool SIMoutput::serialize (std::map<std::string,std::string>&) const
{
  std::cerr <<" *** SIMoutput::serialize: Must be implemented in sub-class.\n"
            <<"     Restart not supported for "<< this->getName() << std::endl;
  return false;
}


bool SIMoutput::writeAddFuncs (int iStep, int& nBlock, int idBlock, double time)
{
  for (const std::pair<std::string,RealFunc*>& func : myAddScalars)
    if (!this->writeGlvF(*func.second, func.first.c_str(), iStep,
                         nBlock, idBlock++, time))
      return false;

  return true;
}
