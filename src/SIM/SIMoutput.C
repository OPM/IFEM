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
#include "AlgEqSystem.h"
#include "AnaSol.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "ElementBlock.h"
#include "VTF.h"
#include "Utilities.h"
#include "tinyxml.h"
#include <fstream>
#include <iomanip>


SIMoutput::~SIMoutput ()
{
  if (myVtf) delete myVtf;
}


void SIMoutput::clearProperties ()
{
  myPoints.clear();
  this->SIMbase::clearProperties();
}


bool SIMoutput::parseOutputTag (const TiXmlElement* elem)
{
  if (myPid == 0)
    std::cout <<"  Parsing <"<< elem->Value() <<">"<< std::endl;

  if (strcasecmp(elem->Value(),"resultpoints"))
    return this->SIMbase::parseOutputTag(elem);

  const TiXmlElement* point = elem->FirstChildElement("point");
  for (int i = 1; point; i++, point = point->NextSiblingElement())
  {
    int patch = 0;
    ResultPoint thePoint;
    if (utl::getAttribute(point,"patch",patch) && patch > 0)
      thePoint.patch = patch;
    if (myPid == 0)
      std::cout <<"\tPoint "<< i <<": P"<< thePoint.patch <<" xi =";
    if (utl::getAttribute(point,"u",thePoint.par[0]) && myPid == 0)
      std::cout <<' '<< thePoint.par[0];
    if (utl::getAttribute(point,"v",thePoint.par[1]) && myPid == 0)
      std::cout <<' '<< thePoint.par[1];
    if (utl::getAttribute(point,"w",thePoint.par[2]) && myPid == 0)
      std::cout <<' '<< thePoint.par[2];
    if (myPid == 0)
      std::cout << std::endl;
    myPoints.push_back(thePoint);
  }
  if (myPid == 0)
    std::cout << std::endl;

  return true;
}


bool SIMoutput::parse (const TiXmlElement* elem)
{
  return this->SIMbase::parse(elem);
}


bool SIMoutput::parse (char* keyWord, std::istream& is)
{
  if (strncasecmp(keyWord,"RESULTPOINTS",12))
    return this->SIMbase::parse(keyWord,is);

  int nres = atoi(keyWord+12);
  if (myPid == 0 && nres > 0)
    std::cout <<"\nNumber of result points: "<< nres;

  char* cline = NULL;
  myPoints.resize(nres);
  for (int i = 0; i < nres && (cline = utl::readLine(is)); i++)
  {
    myPoints[i].patch = atoi(strtok(cline," "));
    if (myPid == 0)
      std::cout <<"\n\tPoint "<< i+1 <<": P"<< myPoints[i].patch <<" xi =";
    for (int j = 0; j < 3 && (cline = strtok(NULL," ")); j++)
    {
      myPoints[i].par[j] = atof(cline);
      if (myPid == 0)
        std::cout <<' '<< myPoints[i].par[j];
    }
  }
  if (myPid == 0 && nres > 0)
    std::cout << std::endl;

  return true;
}


void SIMoutput::preprocessResultPoints ()
{
  for (ResPointVec::iterator p = myPoints.begin(); p != myPoints.end();)
  {
    int pid = this->getLocalPatchIndex(p->patch);
    if (pid < 1 || myModel[pid-1]->empty())
      p = myPoints.erase(p);
    else if ((p->inod = myModel[pid-1]->evalPoint(p->par,p->par,p->X)) < 0)
      p = myPoints.erase(p);
    else
    {
      p->npar = myModel[pid-1]->getNoParamDim();
      int ipt = 1 + (int)(p-myPoints.begin());
      if (ipt == 1) std::cout <<'\n';
      std::cout <<"Result point #"<< ipt <<": patch #"<< p->patch;
      switch (p->npar) {
      case 1: std::cout <<" u="; break;
      case 2: std::cout <<" (u,v)=("; break;
      case 3: std::cout <<" (u,v,w)=("; break;
      }
      std::cout << p->par[0];
      for (unsigned char c = 1; c < p->npar; c++)
        std::cout <<','<< p->par[c];
      if (p->npar > 1) std::cout <<')';
      if (p->inod > 0) std::cout <<", node #"<< p->inod;
      std::cout <<", X = "<< p->X << std::endl;
      p++;
    }
  }
}


bool SIMoutput::writeGlvG (int& nBlock, const char* inpFile, bool doClear)
{
  if (inpFile)
  {
    if (myVtf) return false;

#if HAS_VTFAPI == 2
    const char* ext = ".vtfx";
#else
    const char* ext = ".vtf";
#endif

    // Open a new VTF-file
    char* vtfName = new char[strlen(inpFile)+10];
    strtok(strcpy(vtfName,inpFile),".");
    if (nProc > 1)
      sprintf(vtfName+strlen(vtfName),"_p%04d%s",myPid,ext);
    else
      strcat(vtfName,ext);

    if (myPid == 0)
      std::cout <<"\nWriting VTF-file "<< vtfName << std::endl;
    myVtf = new VTF(vtfName,opt.format);
    delete[] vtfName;
  }
  else if (doClear)
    myVtf->clearGeometryBlocks();

  ElementBlock* lvb;
  char pname[32];
  size_t i;

  // Convert and write model geometry
  for (i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    if (msgLevel > 0)
      std::cout <<"Writing geometry for patch "<< i+1 << std::endl;
    size_t nd = myModel[i]->getNoParamDim();
    lvb = new ElementBlock(nd == 3 ? 8 : (nd == 2 ? 4 : 2));
    if (!myModel[i]->tesselate(*lvb,opt.nViz))
      return false;

    sprintf(pname,"Patch %ld",myModel[i]->idx);
    if (!myVtf->writeGrid(lvb,pname,++nBlock))
      return false;
  }

  // Additional geometry for immersed boundaries
  for (i = 0; i < myModel.size(); i++)
    if ((lvb = myModel[i]->immersedGeometry()))
    {
      sprintf(pname,"Immersed boundary %ld",i+1);
      if (!myVtf->writeGrid(lvb,pname,++nBlock))
        return false;
    }

  // Do not write the geometry blocks to file yet, writeVectors might create
  // an additional block for the point vectors, see method VTF::writeVectors
  return true;
}


bool SIMoutput::writeGlvBC (int& nBlock, int iStep) const
{
  if (!myVtf) return false;

  Matrix field;
  IntVec dID[3];

  size_t i, j;
  int node, geomID = myGeomID;
  for (i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    geomID++;
    size_t nbc = myModel[i]->getNoFields(1);
    Matrix bc(nbc,myModel[i]->getNoNodes(-1));
    RealArray flag(3,0.0);
    ASMbase::BCVec::const_iterator bit;
    for (bit = myModel[i]->begin_BC(); bit != myModel[i]->end_BC(); bit++)
      if ((node = myModel[i]->getNodeIndex(bit->node,true)))
      {
        if (!bit->CX && nbc > 0) bc(1,node) = flag[0] = 1.0;
        if (!bit->CY && nbc > 1) bc(2,node) = flag[1] = 1.0;
        if (!bit->CZ && nbc > 2) bc(3,node) = flag[2] = 1.0;
      }

    if (flag[0]+flag[1]+flag[2] == 0.0)
      continue; // nothing on this patch

    if (msgLevel > 1)
      std::cout <<"Writing boundary conditions for patch "<< i+1 << std::endl;

    if (!myModel[i]->evalSolution(field,bc,opt.nViz))
      return false;

    // The BC fields should either be 0.0 or 1.0
    if (opt.nViz[0] > 2 || opt.nViz[1] > 2 || opt.nViz[2] > 2)
      for (j = 1; j <= 3; j++)
        if (flag[j-1] == 1.0)
          for (size_t n = 1; n <= field.cols(); n++)
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
  if (myProblem->hasTractionValues())
  {
    if (msgLevel > 1)
      std::cout <<"Writing boundary tractions"<< std::endl;
    return myProblem->writeGlvT(myVtf,iStep,geoBlk,nBlock);
  }

  return true;
}


bool SIMoutput::writeGlvV (const Vector& vec, const char* fieldName,
                           int iStep, int& nBlock, int idBlock) const
{
  if (vec.empty())
    return true;
  else if (!myVtf)
    return false;

  Matrix field;
  Vector lovec;
  IntVec vID;

  int geomID = myGeomID;
  for (size_t i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    if (msgLevel > 1)
      std::cout <<"Writing vector field for patch "<< i+1 << std::endl;

    myModel[i]->extractNodeVec(vec,lovec);
    if (!myModel[i]->evalSolution(field,lovec,opt.nViz))
      return false;

    if (!myVtf->writeVres(field,++nBlock,++geomID,this->getNoSpaceDim()))
      return false;
    else
      vID.push_back(nBlock);
  }

  return myVtf->writeVblk(vID,fieldName,idBlock,iStep);
}


/*!
  This method writes both the primary and all secondary solution fields to
  the VTF-file. The primary solution is written as a deformation plot (labelled
  "Solution") if \a pvecName is NULL. If the primary solution is a scalar field,
  the field value is then interpreted as a deformation along the global Z-axis.
  If the primary solution is a vector field and \a pvecName is not NULL, it is
  written as a named vector field instead (no deformation plot).

  If the primary solution is a vector field, each vector component is written
  as a scalar field in addition. If \a psolOnly is greater than 1, it is only
  written as scalar field components (no deformation or vector field output).

  If analytical solution fields are available, those fields are written as well.
*/

bool SIMoutput::writeGlvS (const Vector& psol, int iStep, int& nBlock,
                           double time, char psolOnly, const char* pvecName,
                           int idBlock, int psolComps)
{
  if (psol.empty())
    return true;
  else if (!myVtf)
    return false;

  if (!psolOnly && myProblem)
    myProblem->initResultPoints(time);

  Matrix field;
  Vector lovec;
  const size_t pMAX = 6;
  const size_t sMAX = 21;
  IntVec vID[3], dID[pMAX], sID[sMAX];
  bool scalarEq = this->getNoFields() == 1 || psolOnly == 3;
  size_t nVcomp = scalarEq ? 1 : this->getNoFields();
  if (nVcomp > this->getNoSpaceDim())
    nVcomp = this->getNoSpaceDim();

  bool haveAsol = false;
  bool haveXsol = false;
  if (mySol)
    if (scalarEq)
    {
      haveAsol = mySol->hasScalarSol() > 1;
      haveXsol = mySol->getScalarSol() != NULL;
    }
    else
    {
      haveAsol = mySol->hasVectorSol() > 1;
      haveXsol = mySol->getVectorSol() != NULL;
    }

  bool project = (opt.discretization == ASM::Spline ||
                  opt.discretization == ASM::SplineC1) &&
                 opt.project.find(SIMoptions::GLOBAL) != opt.project.end();

  size_t i, j, k;
  int geomID = myGeomID;
  for (i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    if (msgLevel > 1)
      std::cout <<"Writing FE solution for patch "<< i+1 << std::endl;

    // 1. Evaluate primary solution variables

    Vector& locvec = myProblem ? myProblem->getSolution() : lovec;
    myModel[i]->extractNodeVec(psol,locvec,psolComps,0);
    this->extractPatchDependencies(myProblem,myModel,i);
    if (!myModel[i]->evalSolution(field,locvec,opt.nViz))
      return false;

    myModel[i]->filterResults(field,myVtf->getBlock(geomID+1));

    if (psolOnly > 1 || (nVcomp == 1 && pvecName))
      geomID++; // skip output as vector field
    else if (!myVtf->writeVres(field,++nBlock,++geomID,nVcomp))
      return false;
    else
      vID[0].push_back(nBlock);

    for (j = 1, k = 0; j <= field.rows() && k < pMAX; j++)
      if (!myVtf->writeNres(field.getRow(j),++nBlock,geomID))
        return false;
      else
        dID[k++].push_back(nBlock);

    if (haveXsol)
    {
      if (msgLevel > 1)
        std::cout <<"Writing exact solution for patch "<< i+1 << std::endl;

      // 1.a Evaluate exact primary solution

      const ElementBlock* grid = myVtf->getBlock(geomID);
      Vec3Vec::const_iterator cit = grid->begin_XYZ();
      field.fill(0.0);
      if (scalarEq)
      {
        const RealFunc& pSol = *mySol->getScalarSol();
        for (j = 1; cit != grid->end_XYZ() && haveXsol; j++, cit++)
          field(1,j) = pSol(Vec4(*cit,time));
      }
      else
      {
        const VecFunc& pSol = *mySol->getVectorSol();
        for (j = 1; cit != grid->end_XYZ() && haveXsol; j++, cit++)
          field.fillColumn(j,pSol(Vec4(*cit,time)).ptr());
      }

      for (j = 1; j <= field.rows() && k < pMAX && haveXsol; j++)
        if (!myVtf->writeNres(field.getRow(j),++nBlock,geomID))
          return false;
        else
          dID[k++].push_back(nBlock);

      if (!myVtf->writeVres(field,++nBlock,geomID,nVcomp))
        return false;
      else
        vID[1].push_back(nBlock);
    }

    if (psolOnly || !myProblem) continue; // skip secondary solution
    if (myProblem->getNoFields(2) < 1) continue; // no secondary solution

    // 2. Direct evaluation of secondary solution variables

    LocalSystem::patch = i;
    if (!myModel[i]->evalSolution(field,*myProblem,opt.nViz))
      return false;

    myModel[i]->filterResults(field,myVtf->getBlock(geomID));

    if (scalarEq && field.rows() == this->getNoSpaceDim())
      if (!myVtf->writeVres(field,++nBlock,geomID))
        return false;
      else
        vID[2].push_back(nBlock);

    for (j = 1, k = 0; j <= field.rows() && k < sMAX; j++)
      if (!myVtf->writeNres(field.getRow(j),++nBlock,geomID))
        return false;
      else
        sID[k++].push_back(nBlock);

    if (project)
    {
      // 3. Projection of secondary solution variables (tensorial splines only)

      if (!myModel[i]->evalSolution(field,*myProblem,opt.nViz,'D'))
        return false;

      myModel[i]->filterResults(field,myVtf->getBlock(geomID));

      for (j = 1; j <= field.rows() && k < sMAX; j++)
        if (!myVtf->writeNres(field.getRow(j),++nBlock,geomID))
          return false;
        else
          sID[k++].push_back(nBlock);
    }

    if (haveAsol)
    {
      // 4. Evaluate analytical solution variables

      if (msgLevel > 1)
        std::cout <<"Writing exact solution for patch "<< i+1 << std::endl;

      const ElementBlock* grid = myVtf->getBlock(geomID);
      Vec3Vec::const_iterator cit = grid->begin_XYZ();
      field.fill(0.0);
      for (j = 1; cit != grid->end_XYZ() && haveAsol; j++, cit++)
      {
        Vec4 Xt(*cit,time);
        if (mySol->hasScalarSol() == 3 || mySol->hasVectorSol() == 3)
          haveAsol = myProblem->evalSol(lovec,*mySol->getStressSol(),Xt);
        else if (scalarEq)
          haveAsol = myProblem->evalSol(lovec,*mySol->getScalarSecSol(),Xt);
        else
          haveAsol = myProblem->evalSol(lovec,*mySol->getVectorSecSol(),Xt);
        if (haveAsol)
          field.fillColumn(j,lovec);
      }

      for (j = 1; j <= field.rows() && k < sMAX && haveAsol; j++)
        if (!myVtf->writeNres(field.getRow(j),++nBlock,geomID))
          return false;
        else
          sID[k++].push_back(nBlock);
    }
  }

  // Write result block identifications

  bool ok = true;
  std::string pname(pvecName ? pvecName : "Solution");

  if (!vID[0].empty())
    if (pvecName)
      ok = myVtf->writeVblk(vID[0],pname.c_str(),idBlock,iStep);
    else
      ok = myVtf->writeDblk(vID[0],pname.c_str(),idBlock,iStep);

  if (ok && !vID[1].empty())
  {
    std::string ename = "Exact " + pname;
    if (pvecName)
      ok = myVtf->writeVblk(vID[1],ename.c_str(),idBlock+1,iStep);
    else
      ok = myVtf->writeDblk(vID[1],ename.c_str(),idBlock+1,iStep);
  }

  if (ok && !vID[2].empty())
    ok = myVtf->writeVblk(vID[2],"Flux",idBlock+2,iStep);

  size_t nf = scalarEq ? 1 : this->getNoFields();
  std::vector<std::string> xname(haveXsol ? nf : 0);
  if (nf > 1) pname += "_w";
  for (i = j = 0; i < nf && j < pMAX && !dID[j].empty() && ok; i++, j++)
  {
    if ((!pvecName || nf > nVcomp) && myProblem)
      pname = myProblem->getField1Name(i);
    else if (nf > 1)
      (*pname.rbegin()) ++;
    ok = myVtf->writeSblk(dID[j],pname.c_str(),++idBlock,iStep);
    if (haveXsol) xname[i] = "Exact " + pname;
  }

  if (haveXsol)
    for (i = 0; i < nf && j < pMAX && !dID[j].empty() && ok; i++, j++)
      ok = myVtf->writeSblk(dID[j],xname[i].c_str(),++idBlock,iStep);

  if (psolOnly || !myProblem || !ok) return ok;

  nf = myProblem->getNoFields(2);
  for (i = j = 0; i < nf && j < sMAX && !sID[j].empty(); i++, j++)
    if (!myVtf->writeSblk(sID[j],myProblem->getField2Name(i,haveAsol?"FE":0),
                          ++idBlock,iStep)) return false;

  if (project)
    for (i = 0; i < nf && j < sMAX && !sID[j].empty(); i++, j++)
      if (!myVtf->writeSblk(sID[j],myProblem->getField2Name(i,"Projected"),
                            ++idBlock,iStep)) return false;

  if (haveAsol)
    for (i = 0; i < nf && j < sMAX && !sID[j].empty(); i++, j++)
      if (!myVtf->writeSblk(sID[j],myProblem->getField2Name(i,"Exact"),
                            ++idBlock,iStep)) return false;

  return true;
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
  for (size_t i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    LocalSystem::patch = i;
    myModel[i]->extractNodeVec(psol,myProblem->getSolution(),psolComps,0);
    this->extractPatchDependencies(myProblem,myModel,i);
    if (!myModel[i]->evalSolution(field,*myProblem,opt.nViz))
      return false;
  }

  return true;
}


bool SIMoutput::writeGlvP (const Vector& ssol, int iStep, int& nBlock,
                           int idBlock, const char* prefix,
                           std::vector<PointValue>* maxVal)
{
  if (ssol.empty())
    return true;
  else if (!myVtf)
    return false;

  Matrix field;
  Vector lovec;
  const size_t nf = myProblem->getNoFields(2);
  IntVec sID[nf];

  size_t i, j;
  int geomID = myGeomID;
  for (i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    if (msgLevel > 1)
      std::cout <<"Writing projected solution for patch "<< i+1 << std::endl;

    // Evaluate the solution variables at the visualization points
    myModel[i]->extractNodeVec(ssol,lovec,nf);
    if (!myModel[i]->evalSolution(field,lovec,opt.nViz))
      return false;

    // Write out to VTF-file as scalar fields
    const ElementBlock* grid = myVtf->getBlock(++geomID);
    for (j = 0; j < field.rows() && j < nf; j++)
    {
      Vector comp(field.getRow(1+j));
      if (!myVtf->writeNres(comp,++nBlock,geomID))
        return false;
      else
        sID[j].push_back(nBlock);

      if (maxVal && j < maxVal->size())
      {
        // Update extremal values
        size_t indx = 0;
        double cmax = comp.normInf(indx);
        PointValue& pv = (*maxVal)[j];
        if (fabs(cmax) > fabs(pv.second))
          pv = std::make_pair(grid->getCoord(indx-1),cmax);
      }
    }
  }

  // Write result block identifications
  for (j = 0; j < nf && !sID[j].empty(); j++)
    if (!myVtf->writeSblk(sID[j],myProblem->getField2Name(j,prefix),
                          ++idBlock,iStep)) return false;

  return true;
}


bool SIMoutput::evalProjSolution (const Vector& ssol,
                                  std::vector<PointValue>& maxVal)
{
  if (ssol.empty())
    return true;
  else if (!myVtf)
    return true;

  Matrix field;
  Vector lovec;

  size_t i, j;
  int geomID = myGeomID;
  for (i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    // Evaluate the solution variables at the visualization points
    myModel[i]->extractNodeVec(ssol,lovec,myProblem->getNoFields(2));
    if (!myModel[i]->evalSolution(field,lovec,opt.nViz))
      return false;

    const ElementBlock* grid = myVtf->getBlock(++geomID);
    for (j = 0; j < field.rows() && j < maxVal.size(); j++)
    {
      // Update extremal values
      size_t indx = 0;
      double cmax = field.getRow(1+j).normInf(indx);
      if (fabs(cmax) > fabs(maxVal[j].second))
        maxVal[j] = std::make_pair(grid->getCoord(indx-1),cmax);
    }
  }

  return true;
}


bool SIMoutput::writeGlvF (const RealFunc& f, const char* fname,
                           int iStep, int& nBlock, int idBlock, double time)
{
  if (!myVtf) return false;

  IntVec sID;

  int geomID = myGeomID;
  for (size_t i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    if (msgLevel > 1)
      std::cout <<"Writing function "<< fname
                <<" for patch "<< i+1 << std::endl;

    if (!myVtf->writeNfunc(f,time,++nBlock,++geomID))
      return false;
    else
      sID.push_back(nBlock);
  }

  return myVtf->writeSblk(sID,fname,idBlock,iStep);
}


bool SIMoutput::writeGlvStep (int iStep, double value, int itype)
{
  myVtf->writeGeometryBlocks(iStep);

  if (itype == 0)
    return myVtf->writeState(iStep,"Time %g",value,itype);
  else
    return myVtf->writeState(iStep,"Step %g",value,itype);
}


bool SIMoutput::writeGlvM (const Mode& mode, bool freq, int& nBlock)
{
  if (mode.eigVec.empty())
    return true;
  else if (!myVtf)
    return false;

  if (msgLevel > 1)
    std::cout <<"Writing eigenvector for Mode "<< mode.eigNo << std::endl;

  Vector displ;
  Matrix field;
  IntVec vID;

  int geomID = myGeomID;
  for (size_t i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches
    if (myModel.size() > 1 && msgLevel > 1) std::cout <<"."<< std::flush;

    geomID++;
    myModel[i]->extractNodeVec(mode.eigVec,displ);
    if (!myModel[i]->evalSolution(field,displ,opt.nViz))
      return false;

    if (!myVtf->writeVres(field,++nBlock,geomID))
      return false;
    else
      vID.push_back(nBlock);
  }
  if (myModel.size() > 1 && msgLevel > 1) std::cout << std::endl;

  int idBlock = 10;
  if (!myVtf->writeDblk(vID,"Mode Shape",idBlock,mode.eigNo))
    return false;

  return myVtf->writeState(mode.eigNo, freq ? "Frequency %g" : "Eigenvalue %g",
                           mode.eigVal, 1);
}


bool SIMoutput::writeGlvN (const Matrix& norms, int iStep, int& nBlock,
                           const char** prefix)
{
  if (norms.empty())
    return true;
  else if (!myVtf)
    return false;

  NormBase* norm = myProblem->getNormIntegrand(mySol);

  Matrix field;
  const size_t maxN = 20;
  IntVec sID[maxN];

  size_t i, j, k, l, m;
  int geomID = myGeomID;
  for (i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    if (msgLevel > 1)
      std::cout <<"Writing element norms for patch "<< i+1 << std::endl;

    const ElementBlock* grid = myVtf->getBlock(++geomID);
    myModel[i]->extractElmRes(norms,field);

    if (grid->getNoElms() > field.cols())
    {
      // Expand the element result array
      Matrix efield(field);
      field.resize(field.rows(),grid->getNoElms());
      for (j = 1; j <= field.cols(); j++)
        field.fillColumn(j,efield.getColumn(grid->getElmId(j)));
    }

    j = l = 1;
    for (k = m = 0; m < field.rows() && k < maxN; m++)
    {
      if (l > norm->getNoFields(j))
        l = 1, ++j;
      if (norm->hasElementContributions(j,l++))
        if (!myVtf->writeEres(field.getRow(1+m),++nBlock,geomID))
          return false;
        else
          sID[k++].push_back(nBlock);
    }
  }

  int idBlock = 200;
  j = l = 1;
  for (k = 0; k < maxN && !sID[k].empty(); l++)
  {
    if (l > norm->getNoFields(j))
      l = 1, ++j;
    if (!norm->hasElementContributions(j,l))
      continue;

    const char* normName = norm->getName(j,l,j>1?prefix[j-2]:0);
    if (!myVtf->writeSblk(sID[k++],normName,++idBlock,iStep,true))
      return false;
  }

  delete norm;
  return true;
}


void SIMoutput::closeGlv ()
{
  if (myVtf) delete myVtf;
  myVtf = 0;
}


bool SIMoutput::dumpGeometry (std::ostream& os) const
{
  for (size_t i = 0; i < myModel.size(); i++)
    if (!myModel[i]->empty())
      if (!myModel[i]->write(os))
        return false;

  return true;
}


bool SIMoutput::dumpBasis (std::ostream& os, int basis, size_t patch) const
{
  size_t begin = patch > 0 ? patch-1 : 0;
  size_t end = patch > 0 ? patch : myModel.size();
  for (size_t i = begin; i < end && i < myModel.size(); i++)
    if (!myModel[i]->empty())
      if (!myModel[i]->write(os,basis))
        return false;

  return true;
}


void SIMoutput::dumpPrimSol (const Vector& psol, std::ostream& os,
                             bool withID) const
{
  if (psol.empty()) return;

  size_t i, j, ip;
  unsigned char k, n;
  for (i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    Vector patchSol;
    myModel[i]->extractNodeVec(psol,patchSol,mySam->getMADOF());

    if (withID)
    {
      if (myModel.size() > 1)
        os <<"\n# Patch: "<< i+1;
      os <<"\n# inod/gnod\tNodal Coordinates\tSolution\n";
    }
    for (ip = 0, j = 1; j <= myModel[i]->getNoNodes(); j++)
    {
      if ((n = myModel[i]->getNodalDOFs(j)) == 0)
        continue;
      else if (withID)
        os << j <<' '<< myModel[i]->getNodeID(j)
           <<"\t\t"<< myModel[i]->getCoord(j) <<"\t\t";
      os << utl::trunc(patchSol[ip++]);
      for (k = 1; k < n; k++)
        os <<' '<< utl::trunc(patchSol[ip++]);
      os <<'\n';
    }
  }

  os.flush();
}


bool SIMoutput::dumpSolution (const Vector& psol, std::ostream& os) const
{
  if (psol.empty())
    return true;

  Matrix field;
  size_t i, j, k;

  for (i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    if (myModel.size() > 1)
      os <<"\n# Patch: "<< i+1;

    // Extract and write primary solution
    size_t nf = myModel[i]->getNoFields(1);
    Vector& patchSol = myProblem->getSolution();
    myModel[i]->extractNodeVec(psol,patchSol,mySam->getMADOF());
    for (k = 1; k <= nf; k++)
    {
      os << myProblem->getField1Name(k,"# FE");
      for (j = 1; j <= myModel[i]->getNoNodes(); j++)
      {
        std::pair<int,int> dofs = mySam->getNodeDOFs(j);
        int idof = dofs.first+k-1;
        if (idof <= dofs.second)
          os <<"\n"<< utl::trunc(patchSol[idof-1]);
      }
      os << std::endl;
    }

    // Evaluate secondary solution variables
    LocalSystem::patch = i;
    if (!myModel[i]->evalSolution(field,*myProblem))
      return false;

    // Write the secondary solution
    for (j = 1; j <= field.rows(); j++)
    {
      os << myProblem->getField2Name(j-1,"# FE");
      for (k = 1; k <= field.cols(); k++)
        os <<"\n"<< utl::trunc(field(j,k));
      os << std::endl;
    }
  }

  return true;
}


bool SIMoutput::dumpResults (const Vector& psol, double time, std::ostream& os,
                             bool formatted, std::streamsize precision) const
{
  if (psol.empty() || myPoints.empty())
    return true;

  myProblem->initResultPoints(time);

  size_t i, j, k;
  Matrix sol1, sol2;
  Vector reactionFS;
  const Vector* reactionForces = myEqSys->getReactions();

  for (i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    ResPointVec::const_iterator p;
    RealArray params[3];
    IntVec points;

    // Find all evaluation points within this patch, if any
    for (j = 0, p = myPoints.begin(); p != myPoints.end(); j++, p++)
      if (this->getLocalPatchIndex(p->patch) == (int)(i+1))
        if (opt.discretization >= ASM::Spline)
        {
          points.push_back(p->inod > 0 ? p->inod : -(j+1));
          for (k = 0; k < myModel[i]->getNoParamDim(); k++)
            params[k].push_back(p->par[k]);
        }
        else if (p->inod > 0)
          points.push_back(p->inod);

    if (points.empty()) continue; // no points in this patch

    myModel[i]->extractNodeVec(psol,myProblem->getSolution());
    if (opt.discretization >= ASM::Spline)
    {
      // Evaluate the primary solution variables
      if (!myModel[i]->evalSolution(sol1,myProblem->getSolution(),params,false))
        return false;

      // Evaluate the secondary solution variables
      LocalSystem::patch = i;
      if (myProblem->getNoFields(2) > 0)
        if (!myModel[i]->evalSolution(sol2,*myProblem,params,false))
          return false;
    }
    else
      // Extract nodal primary solution variables
      if (!myModel[i]->getSolution(sol1,myProblem->getSolution(),points))
        return false;

    // Formatted output, use scientific notation with fixed field width
    std::streamsize flWidth = 8 + precision;
    std::streamsize oldPrec = os.precision(precision);
    std::ios::fmtflags oldF = os.flags(std::ios::scientific | std::ios::right);
    for (j = 0; j < points.size(); j++)
    {
      if (!formatted)
        os << time <<" ";
      else if (points[j] < 0)
        os <<"  Point #"<< -points[j] <<":\tsol1 =";
      else
      {
        points[j] = myModel[i]->getNodeID(points[j]);
        os <<"  Node #"<< points[j] <<":\tsol1 =";
      }

      for (k = 1; k <= sol1.rows(); k++)
        os << std::setw(flWidth) << utl::trunc(sol1(k,j+1));

      if (opt.discretization >= ASM::Spline)
      {
        if (formatted && sol2.rows() > 0)
          os <<"\n\t\tsol2 =";
        for (k = 1; k <= sol2.rows(); k++)
          os << std::setw(flWidth) << utl::trunc(sol2(k,j+1));
      }

      if (reactionForces && points[j] > 0)
        // Print nodal reaction forces for nodes with prescribed DOFs
        if (mySam->getNodalReactions(points[j],*reactionForces,reactionFS))
        {
          if (formatted)
            os <<"\n\t\treac =";
          for (k = 0; k < reactionFS.size(); k++)
            os << std::setw(flWidth) << utl::trunc(reactionFS[k]);
        }

      os << std::endl;
    }
    os.precision(oldPrec);
    os.flags(oldF);
  }

  return true;
}


bool SIMoutput::savePoints (const std::string& fileName,
                            const Vector& psol, double time, int step,
                            std::streamsize precision) const
{
  if (step < 1 || myPoints.empty() || fileName.empty())
    return true;

  // Append extension .dat if nothing yet
  std::string datafile(fileName);
  size_t idot = datafile.find_last_of('.');
  if (idot == std::string::npos)
  {
    idot = datafile.size();
    datafile.append(".dat");
  }

  if (step == 1)
  {
    // Dump result point coordinates to file
    std::string coordfile(datafile);
    coordfile.insert(idot,"_coord");

    // Formatted output, use scientific notation with fixed field width
    std::streamsize flWidth = 8 + precision;
    std::ofstream f(coordfile.c_str(),std::ios::out);
    f.flags(std::ios::scientific | std::ios::right);
    f.precision(precision);

    for (size_t i = 0; i < myPoints.size(); i++)
    {
      f << time <<" ";
      for (unsigned char k = 0; k < myPoints[i].npar; k++)
        f << std::setw(flWidth) << myPoints[i].X[k];
      f << std::endl;
    }
  }

  std::ofstream f(datafile.c_str(), step == 1 ? std::ios::out : std::ios::app);
  return this->dumpResults(psol,time,f,false,precision);
}
