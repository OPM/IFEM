// $Id$
//==============================================================================
//!
//! \file VTF.C
//!
//! \date Dec 1 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Output of FE model and results to VTF file.
//!
//==============================================================================

#include "VTF.h"
#if HAS_VTFAPI == 1
#include "VTFAPI.h"
#include "VTOAPIPropertyIDs.h"
#elif HAS_VTFAPI == 2
#include "VTFXAPI.h"
#include "VTOAPIPropertyIDs.h"
#define VTFA_FAILURE VTFXA_FAILURE
#define VTFA_SUCCESS VTFXA_SUCCESS
#else
#define VTFA_FAILURE(x) x <= 0
#define VTFA_SUCCESS(x) x > 0
#endif
#include "ElementBlock.h"
#include "Tensor.h"
#include <iostream>
#include <cstdio>


Real VTF::vecOffset[3] = { 0.0, 0.0, 0.0 };


VTF::VTF (const char* filename, int type)
{
  myFile = 0;
  myState = 0;
  myGBlock = 0;
  pointGeoID = 0;
  if (!filename) return;

#if HAS_VTFAPI == 1
  // Create the VTF file object
  myFile = new VTFAFile();
  // Enable debug info to stderr/console
  myFile->SetOutputDebugError(1);

  if (!VTFA_FAILURE(myFile->CreateVTFFile(filename,type > 0)))
    return;

  delete myFile;
  showError("Error creating VTF file");
#elif HAS_VTFAPI == 2
  myFile = new VTFXAFile();
  VTFXAFileSettings settings;
  VTFXAInitFileSettings(&settings);
  settings.bBinary = type > 0 ? VTFXA_TRUE : VTFXA_FALSE;
  settings.pszApplicationName = "IFEM";
  settings.pszVendorName = "SINTEF ICT";
  settings.iVendorID = 1001;
  settings.pszVendorCode = "Test";
  if (!VTFA_FAILURE(myFile->CreateVTFxFile(filename,&settings))) {
    myDatabase = new VTFXADatabase(myFile,"Single",1);
    return;
  }
  delete myFile;
  showError("Error creating VTFx file");
#else
  showError("Not available in this version");
#endif
  myFile = 0;
}


VTF::~VTF ()
{
  if (!myFile) return;

  size_t i;
  for (i = 0; i < myBlocks.size(); i++)
    delete myBlocks[i].second;

#if HAS_VTFAPI == 1
  if (myGBlock)
  {
    if (VTFA_FAILURE(myFile->WriteBlock(myGBlock)))
      showError("Error writing Geometry Block");
    delete myGBlock;
  }

  for (i = 0; i < myTBlock.size(); i++)
    if (myTBlock[i])
    {
      if (VTFA_FAILURE(myFile->WriteBlock(myTBlock[i])))
	showError("Error writing Transformation Block");
      delete myTBlock[i];
    }

  for (i = 0; i < myDBlock.size(); i++)
    if (myDBlock[i])
    {
      if (VTFA_FAILURE(myFile->WriteBlock(myDBlock[i])))
	showError("Error writing Displacement Block");
      delete myDBlock[i];
    }

  for (i = 0; i < myVBlock.size(); i++)
    if (myVBlock[i])
    {
      if (VTFA_FAILURE(myFile->WriteBlock(myVBlock[i])))
	showError("Error writing Vector Block");
      delete myVBlock[i];
    }

  for (i = 0; i < mySBlock.size(); i++)
    if (mySBlock[i])
    {
      if (VTFA_FAILURE(myFile->WriteBlock(mySBlock[i])))
	showError("Error writing Scalar Block");
      delete mySBlock[i];
    }

  if (myState)
  {
    if (VTFA_FAILURE(myFile->WriteBlock(myState)))
      showError("Error writing state info block");
    delete myState;
  }

  if (VTFA_FAILURE(myFile->CloseFile()))
    showError("Error closing VTF file");

  delete myFile;
#elif HAS_VTFAPI == 2
  if (myGBlock)
  {
    if (VTFA_FAILURE(myDatabase->WriteBlock(myGBlock)))
      showError("Error writing Geometry Block");
    delete myGBlock;
  }

  for (i = 0; i < myTBlock.size(); i++)
    if (myTBlock[i])
    {
      if (VTFA_FAILURE(myDatabase->WriteBlock(myTBlock[i])))
        showError("Error writing Transformation Block");
      delete myTBlock[i];
    }

  for (i = 0; i < myDBlock.size(); i++)
    if (myDBlock[i])
    {
      if (VTFA_FAILURE(myDatabase->WriteBlock(myDBlock[i])))
        showError("Error writing Displacement Block");
      delete myDBlock[i];
    }

  for (i = 0; i < myVBlock.size(); i++)
    if (myVBlock[i])
    {
      if (VTFA_FAILURE(myDatabase->WriteBlock(myVBlock[i])))
        showError("Error writing Vector Block");
      delete myVBlock[i];
    }

  for (i = 0; i < mySBlock.size(); i++)
    if (mySBlock[i])
    {
      if (VTFA_FAILURE(myDatabase->WriteBlock(mySBlock[i])))
	showError("Error writing Scalar Block");
      delete mySBlock[i];
    }

  if (myState)
  {
    if (VTFA_FAILURE(myDatabase->WriteBlock(myState)))
      showError("Error writing state info block");
    delete myState;
  }

  VTFXACase                singleCase(myFile,"Single case",1,1);
  VTFXACasePropertiesBlock frameGeneratorProps(VT_CT_FRAME_GENERATOR_SETTINGS);

  frameGeneratorProps.AddInt(VT_PI_FG_FEM_MODEL_IDS,1);
  singleCase.WritePropertiesBlock(&frameGeneratorProps);
  for (i = 0; i < myBlocks.size(); i++) {
    VTFXACasePropertiesBlock partAttr(VT_CT_PART_ATTRIBUTES);
    partAttr.SetPartID(i+1);
    partAttr.AddBool(VT_PB_PA_MESH, VTFXA_FALSE);
    partAttr.AddBool(VT_PB_PA_DISPLACEMENTS, VTFXA_FALSE);
    singleCase.WritePropertiesBlock(&partAttr);
  }
  if (VTFA_FAILURE(myFile->CloseFile()))
    showError("Error closing VTF file");

  delete myDatabase;
  delete myFile;
#endif
}


void VTF::writeGeometryBlocks (int iStep)
{
  if (myBlocks.empty())
    return;

  std::vector<int> geomID(myBlocks.size());
  for (size_t i = 0; i < myBlocks.size(); i++)
    geomID[i] = myBlocks[i].first;

#ifdef HAS_VTFAPI
  if (!myGBlock) myGBlock = new VTFAGeometryBlock();
#if HAS_VTFAPI == 1
  myGBlock->SetGeometryElementBlocks(&geomID.front(),geomID.size(),iStep);
#elif HAS_VTFAPI == 2
  myGBlock->SetElementBlocksForState(&geomID.front(),geomID.size(),iStep);
#endif
#endif
}


void VTF::clearGeometryBlocks ()
{
  for (size_t i = 0; i < myBlocks.size(); i++)
    delete myBlocks[i].second;

  myBlocks.clear();
  pointGeoID = 0;
}


const ElementBlock* VTF::getBlock (int geomID) const
{
  if (geomID > 0 && (size_t)geomID <= myBlocks.size())
    return myBlocks[geomID-1].second;

  // Return a pointer to an empty block if index is out of range
  static ElementBlock emptyBlock;
  return &emptyBlock;
}


bool VTF::writeGrid (const ElementBlock* block, const char* partname, int gID)
{
  if (!myFile) return true;
  if (!block) return false;

  myBlocks.push_back(std::make_pair(gID,block));

  if (!this->writeNodes(gID))
    return showError("Error writing node block",gID);

  if (!this->writeElements(partname,gID,gID))
    return showError("Error writing element block",gID);

  return true;
}


bool VTF::writeTransformation (const Vec3& X, const Tensor& T,
			       int idBlock, int geomID)
{
  if (!myFile) return true;
  if (!this->getBlock(geomID)) return false;

#if HAS_VTFAPI == 1
  // Cast to float
  float trans[12];
  size_t i, j, k = 0;
  for (j = 1; j <= 3; j++)
    for (i = 1; i <= 3; i++)
      trans[k++] = T(i,j);
  for (j = 0; j < 3; j++)
    trans[k++] = X[j];

  VTFAMatrixResultBlock tBlock(idBlock);
  if (VTFA_FAILURE(tBlock.SetMatrix(trans)))
    return showError("Error defining result block",idBlock);

  tBlock.SetMapToElementBlockID(myBlocks[geomID-1].first);
  if (VTFA_FAILURE(myFile->WriteBlock(&tBlock)))
    return showError("Error writing result block",idBlock);
#elif HAS_VTFAPI == 2
  std::cerr <<"VTF: Transformation not yet implemented for VTFx"<< std::endl;
#endif

  return true;
}


bool VTF::writeVres (const std::vector<Real>& nodeResult,
		     int idBlock, int geomID, size_t nvc)
{
  if (!myFile) return true;

  const ElementBlock* grid = this->getBlock(geomID);
  if (!grid) return false;

  const size_t nnod = grid->getNoNodes();
  const size_t nres = nodeResult.size();
  const size_t ncmp = nres/(nnod > 0 ? nnod : 1);
  if (nres != ncmp*nnod)
    return showError("Invalid size of result array",nres);
  else if (nvc < 1 || nvc > ncmp)
    nvc = ncmp;

  // Cast to float
  float* resVec = new float[3*nnod];
  if (nres == 3*nnod && nvc == 3)
    for (size_t i = 0; i < nres; i++)
      resVec[i] = nodeResult[i];
  else if (nres == nnod && nvc == 1)
    for (size_t i = 0; i < nnod; i++)
    {
      // Writing a scalar solution as Z-deflection
      resVec[3*i] = resVec[3*i+1] = 0.0f;
      resVec[3*i+2] = nodeResult[i];
    }
  else
    for (size_t i = 0; i < nnod; i++)
      for (size_t j = 0; j < 3; j++)
	resVec[3*i+j] = j < nvc ? nodeResult[ncmp*i+j] : 0.0f;

#if HAS_VTFAPI == 1
  VTFAResultBlock dBlock(idBlock,VTFA_DIM_VECTOR,VTFA_RESMAP_NODE,0);

  if (VTFA_FAILURE(dBlock.SetResults3D(resVec,nnod)))
    return showError("Error defining result block",idBlock);

  dBlock.SetMapToBlockID(myBlocks[geomID-1].first);
  if (VTFA_FAILURE(myFile->WriteBlock(&dBlock)))
    return showError("Error writing result block",idBlock);
#elif HAS_VTFAPI == 2
  VTFXAResultValuesBlock dBlock(idBlock,VTFXA_DIM_VECTOR,VTFXA_FALSE);
  dBlock.SetMapToBlockID(myBlocks[geomID-1].first,VTFXA_NODES);
  dBlock.SetResultValues3D(resVec,nnod);
  if (VTFA_FAILURE(myDatabase->WriteBlock(&dBlock)))
    return showError("Error writing result block",idBlock);
#endif
  delete[] resVec;

  return true;
}


bool VTF::writeEres (const std::vector<Real>& elementResult,
		     int idBlock, int geomID)
{
  if (!myFile) return true;

  const ElementBlock* grid = this->getBlock(geomID);
  if (!grid) return false;

  const size_t nres = elementResult.size();
  if (nres > grid->getNoElms())
    return showError("Invalid size of result array",nres);
  else if (nres < grid->getNoElms())
    showError("Warning: Fewer element results that anticipated",nres);

  // Cast to float
  float* resVec = new float[nres];
  for (size_t i = 0; i < nres; i++)
    resVec[i] = elementResult[i];

#if HAS_VTFAPI == 1
  VTFAResultBlock dBlock(idBlock,VTFA_DIM_SCALAR,VTFA_RESMAP_ELEMENT,0);

  if (VTFA_FAILURE(dBlock.SetResults1D(resVec,nres)))
    return showError("Error defining result block",idBlock);

  dBlock.SetMapToBlockID(myBlocks[geomID-1].first);
  if (VTFA_FAILURE(myFile->WriteBlock(&dBlock)))
    return showError("Error writing result block",idBlock);
#elif HAS_VTFAPI == 2
  VTFXAResultValuesBlock dBlock(idBlock,VTFXA_DIM_SCALAR,VTFXA_FALSE);
  dBlock.SetMapToBlockID(myBlocks[geomID-1].first,VTFXA_ELEMENTS);
  dBlock.SetResultValues1D(resVec,nres);
  if (VTFA_FAILURE(myDatabase->WriteBlock(&dBlock)))
    return showError("Error writing result block",idBlock);
#endif

  delete[] resVec;
  return true;
}


bool VTF::writeNres (const std::vector<Real>& nodalResult,
		     int idBlock, int geomID)
{
  if (!myFile) return true;

  const ElementBlock* grid = this->getBlock(geomID);
  if (!grid) return false;

  const size_t nres = nodalResult.size();
  if (nres != grid->getNoNodes())
    return showError("Invalid size of result array",nres);

  // Cast to float
  float* resVec = new float[nres];
  for (size_t i = 0; i < nres; i++)
    resVec[i] = nodalResult[i];

#if HAS_VTFAPI == 1
  VTFAResultBlock dBlock(idBlock,VTFA_DIM_SCALAR,VTFA_RESMAP_NODE,0);

  if (VTFA_FAILURE(dBlock.SetResults1D(resVec,nres)))
    return showError("Error defining result block",idBlock);

  dBlock.SetMapToBlockID(myBlocks[geomID-1].first);
  if (VTFA_FAILURE(myFile->WriteBlock(&dBlock)))
    return showError("Error writing result block",idBlock);
#elif HAS_VTFAPI == 2
  VTFXAResultValuesBlock dBlock(idBlock,VTFXA_DIM_SCALAR,VTFXA_FALSE);
  dBlock.SetMapToBlockID(myBlocks[geomID-1].first,VTFXA_NODES);
  dBlock.SetResultValues1D(resVec,nres);
  if (VTFA_FAILURE(myDatabase->WriteBlock(&dBlock)))
    return showError("Error writing result block",idBlock);
#endif
  delete[] resVec;

  return true;
}


bool VTF::writeNfunc (const RealFunc& f, Real time, int idBlock, int geomID)
{
  if (!myFile) return true;

  const ElementBlock* grid = this->getBlock(geomID);
  if (!grid) return false;

  const size_t nres = grid->getNoNodes();

  // Evaluate the function at the grid points
  float* resVec = new float[nres];
  std::vector<Vec3>::const_iterator cit = grid->begin_XYZ();
  for (size_t i = 0; i < nres; i++, cit++)
    resVec[i] = f(Vec4(*cit,time));

#if HAS_VTFAPI == 1
  VTFAResultBlock dBlock(idBlock,VTFA_DIM_SCALAR,VTFA_RESMAP_NODE,0);

  if (VTFA_FAILURE(dBlock.SetResults1D(resVec,nres)))
    return showError("Error defining result block",idBlock);

  dBlock.SetMapToBlockID(myBlocks[geomID-1].first);
  if (VTFA_FAILURE(myFile->WriteBlock(&dBlock)))
    return showError("Error writing result block",idBlock);
#elif HAS_VTFAPI == 2
  VTFXAResultValuesBlock dBlock(idBlock,VTFXA_DIM_SCALAR,VTFXA_FALSE);
  dBlock.SetMapToBlockID(myBlocks[geomID-1].first,VTFXA_NODES);
  dBlock.SetResultValues1D(resVec,nres);
  if (VTFA_FAILURE(myDatabase->WriteBlock(&dBlock)))
    return showError("Error writing result block",idBlock);
#endif
  delete[] resVec;

  return true;
}


bool VTF::writeVectors (const std::vector<Vec3Pair>& pntResult, int& gID,
                        int idBlock, const char* resultName,
                        int iStep, int iBlock)
{
#if HAS_VTFAPI == 1
  bool writePoints = false;
  if (pointGeoID == 0)
  {
    // The big assumption here is that we have only one call to writeVectors
    // per time step, and that all subsequent calls are with the same points
    pointGeoID = gID < 0 && !myBlocks.empty() ? myBlocks.back().first+1 : ++gID;
    myBlocks.push_back(std::make_pair(pointGeoID,new ElementBlock()));
    writePoints = true;
  }

  VTFANodeBlock   nBlock(pointGeoID,0);
  VTFAResultBlock rBlock(idBlock,VTFA_DIM_VECTOR,VTFA_RESMAP_NODE,0);

  size_t i = 0, np = pntResult.size();
  if (writePoints && VTFA_FAILURE(nBlock.SetNumNodes(np)))
    return showError("Error defining node block",pointGeoID);
  else if (VTFA_FAILURE(rBlock.SetNumResults(np)))
    return showError("Error defining result block",idBlock);
  else
    rBlock.SetMapToBlockID(pointGeoID);

  int* mnpc = writePoints ? new int[np] : 0;
  std::vector<Vec3Pair>::const_iterator cit;
  for (cit = pntResult.begin(); cit != pntResult.end(); cit++, i++)
    if (writePoints && VTFA_FAILURE(nBlock.AddNode(vecOffset[0]+cit->first.x,
						   vecOffset[1]+cit->first.y,
						   vecOffset[2]+cit->first.z)))
      return showError("Error adding node to block",pointGeoID);
    else if (VTFA_FAILURE(rBlock.AddResult(cit->second.x,
					   cit->second.y,
					   cit->second.z)))
      return showError("Error adding result to block",idBlock);
    else if (writePoints)
      mnpc[i] = i;

  if (writePoints)
  {
    // We must define an element block (with point elements) also,
    // otherwise GLview does not visualize the vectors
    VTFAElementBlock eBlock(pointGeoID,0,0);
    eBlock.SetPartID(pointGeoID);
    eBlock.SetNodeBlockID(pointGeoID);
    if (VTFA_FAILURE(eBlock.AddElements(VTFA_POINTS,mnpc,np)))
      return showError("Error defining element block",pointGeoID);
    delete[] mnpc;

    if (VTFA_FAILURE(myFile->WriteBlock(&nBlock)))
      return showError("Error writing node block",pointGeoID);
    else if (VTFA_FAILURE(myFile->WriteBlock(&eBlock)))
      return showError("Error writing element block",pointGeoID);
  }
  if (VTFA_FAILURE(myFile->WriteBlock(&rBlock)))
    return showError("Error writing result block",idBlock);
#elif HAS_VTFAPI == 2
  std::cerr <<"VTF: Vector points are not yet implemented for VTFx"<< std::endl;
  return true;
#endif

  return iStep > 0 ? this->writeVblk(idBlock,resultName,iBlock,iStep) : true;
}


bool VTF::writePoints (const Vec3Vec& points, int& gID)
{
#if HAS_VTFAPI == 1
  myBlocks.push_back(std::make_pair(++gID,new ElementBlock()));

  VTFANodeBlock nBlock(gID,0);

  size_t i, np = points.size();
  if (VTFA_FAILURE(nBlock.SetNumNodes(np)))
    return showError("Error defining node block",gID);

  int* mnpc = new int[np];
  Vec3Vec::const_iterator cit;
  for (cit = points.begin(), i = 0; cit != points.end(); cit++, i++)
    if (VTFA_FAILURE(nBlock.AddNode(vecOffset[0]+cit->x,
				    vecOffset[1]+cit->y,
				    vecOffset[2]+cit->z)))
      return showError("Error adding node to block",gID);
    else
      mnpc[i] = i;

  // We must define an element block (with point elements) also,
  // otherwise GLview does not visualize the points
  VTFAElementBlock eBlock(gID,0,0);
  eBlock.SetPartID(gID);
  eBlock.SetNodeBlockID(gID);
  if (VTFA_FAILURE(eBlock.AddElements(VTFA_POINTS,mnpc,np)))
    return showError("Error defining element block",gID);
  delete[] mnpc;

  if (VTFA_FAILURE(myFile->WriteBlock(&nBlock)))
    return showError("Error writing node block",gID);
  else if (VTFA_FAILURE(myFile->WriteBlock(&eBlock)))
    return showError("Error writing element block",gID);
#elif HAS_VTFAPI == 2
  std::cerr <<"VTF: Points are not yet implemented for VTFx"<< std::endl;
  return true;
#endif

  return true;
}


bool VTF::writeTblk (const std::vector<int>& tBlockIDs, const char* resultName,
		     int idBlock, int iStep)
{
  if ((int)myTBlock.size() < idBlock) myTBlock.resize(idBlock,0);

  int status = 1;
#if HAS_VTFAPI == 1
  if (!myTBlock[--idBlock])
  {
    myTBlock[idBlock] = new VTFATransformationBlock(idBlock+1);
    if (resultName) myTBlock[idBlock]->SetName(resultName);
    status = myTBlock[idBlock]->SetResultBlocks(&tBlockIDs.front(),
                                                tBlockIDs.size(),iStep);
  }
  else for (size_t i = 0; i < tBlockIDs.size() && VTFA_SUCCESS(status); i++)
    status = myTBlock[idBlock]->AddResultBlock(tBlockIDs[i],iStep);
#elif HAS_VTFAPI == 2
  std::cerr <<"VTF: Transformation not yet implemented for VTFx"<< std::endl;
#endif
  if (VTFA_FAILURE(status))
    return showError("Error defining transformation block",idBlock);

  return true;
}


bool VTF::writeDblk (const std::vector<int>& dBlockIDs, const char* resultName,
		     int idBlock, int iStep)
{
  if ((int)myDBlock.size() < idBlock) myDBlock.resize(idBlock,0);

  int status = 1;
#if HAS_VTFAPI == 1
  if (!myDBlock[--idBlock])
  {
    myDBlock[idBlock] = new VTFADisplacementBlock(idBlock+1);
    if (resultName) myDBlock[idBlock]->SetName(resultName);
    myDBlock[idBlock]->SetRelativeDisplacementResults(1);
    status = myDBlock[idBlock]->SetResultBlocks(&dBlockIDs.front(),
                                                dBlockIDs.size(),iStep);
  }
  else for (size_t i = 0; i < dBlockIDs.size() && VTFA_SUCCESS(status); i++)
    status = myDBlock[idBlock]->AddResultBlock(dBlockIDs[i],iStep);
#elif HAS_VTFAPI == 2
  if (!myDBlock[--idBlock])
  {
    myDBlock[idBlock] = new VTFXAResultBlock(idBlock+1,
					     VTFXA_RESTYPE_DISPLACEMENT,
					     VTFXA_RESMAP_NODE);
    if (resultName) myDBlock[idBlock]->SetName(resultName);
  }
  myDBlock[idBlock]->SetResultID(idBlock);
  status = myDBlock[idBlock]->SetResultValuesBlocks(&dBlockIDs.front(),
                                                    dBlockIDs.size(),iStep);
#endif
  if (VTFA_FAILURE(status))
    return showError("Error defining displacement block",idBlock);

  return true;
}


bool VTF::writeVblk (int vBlockID, const char* resultName,
		     int idBlock, int iStep)
{
  if ((int)myVBlock.size() < idBlock) myVBlock.resize(idBlock,0);

  int status = 1;
#if HAS_VTFAPI == 1
  if (!myVBlock[--idBlock])
  {
    myVBlock[idBlock] = new VTFAVectorBlock(idBlock+1);
    if (resultName) myVBlock[idBlock]->SetName(resultName);
    status = myVBlock[idBlock]->SetResultBlocks(&vBlockID,1,iStep);
  }
  else
    status = myVBlock[idBlock]->AddResultBlock(vBlockID,iStep);
#elif HAS_VTFAPI == 2
  if (!myVBlock[--idBlock])
  {
    myVBlock[idBlock] = new VTFXAResultBlock(idBlock+1,
					     VTFXA_RESTYPE_VECTOR,
					     VTFXA_RESMAP_NODE);
    if (resultName) myVBlock[idBlock]->SetName(resultName);
  }
  myVBlock[idBlock]->SetResultID(idBlock);
  status = myVBlock[idBlock]->SetResultValuesBlocks(&vBlockID,1,iStep);
#endif
  if (VTFA_FAILURE(status))
    return showError("Error defining vector block",idBlock);

  return true;
}


bool VTF::writeVblk (const std::vector<int>& vBlockIDs, const char* resultName,
		     int idBlock, int iStep)
{
  if ((int)myVBlock.size() < idBlock) myVBlock.resize(idBlock,0);

  int status = 1;
#if HAS_VTFAPI == 1
  if (!myVBlock[--idBlock])
  {
    myVBlock[idBlock] = new VTFAVectorBlock(idBlock+1);
    if (resultName) myVBlock[idBlock]->SetName(resultName);
    status = myVBlock[idBlock]->SetResultBlocks(&vBlockIDs.front(),
                                                vBlockIDs.size(),iStep);
  }
  else for (size_t i = 0; i < vBlockIDs.size() && VTFA_SUCCESS(status); i++)
    status = myVBlock[idBlock]->AddResultBlock(vBlockIDs[i],iStep);
#elif HAS_VTFAPI == 2
  if (!myVBlock[--idBlock])
  {
    myVBlock[idBlock] = new VTFXAResultBlock(idBlock+1,
					     VTFXA_RESTYPE_VECTOR,
					     VTFXA_RESMAP_NODE);
    if (resultName) myVBlock[idBlock]->SetName(resultName);
  }
  myVBlock[idBlock]->SetResultID(idBlock);
  status = myVBlock[idBlock]->SetResultValuesBlocks(&vBlockIDs.front(),
                                                    vBlockIDs.size(),iStep);
#endif
  if (VTFA_FAILURE(status))
    return showError("Error defining vector block",idBlock);

  return true;
}


bool VTF::writeSblk (int sBlockID, const char* resultName,
		     int idBlock, int iStep, bool elementData)
{
  if ((int)mySBlock.size() < idBlock) mySBlock.resize(idBlock,0);

  int status = 1;
#if HAS_VTFAPI == 1
  if (!mySBlock[--idBlock])
  {
    mySBlock[idBlock] = new VTFAScalarBlock(idBlock+1);
    if (resultName) mySBlock[idBlock]->SetName(resultName);
    status = mySBlock[idBlock]->SetResultBlocks(&sBlockID,1,iStep);
  }
  else
    status = mySBlock[idBlock]->AddResultBlock(sBlockID,iStep);
#elif HAS_VTFAPI == 2
  if (!mySBlock[--idBlock])
  {
    int type = elementData?VTFXA_RESMAP_ELEMENT:VTFXA_RESMAP_NODE;
    mySBlock[idBlock] = new VTFXAResultBlock(idBlock+1,
					     VTFXA_RESTYPE_SCALAR,type);
    if (resultName) mySBlock[idBlock]->SetName(resultName);
  }
  mySBlock[idBlock]->SetResultID(idBlock);
  status = mySBlock[idBlock]->SetResultValuesBlocks(&sBlockID,1,iStep);
#endif
  if (VTFA_FAILURE(status))
    return showError("Error defining scalar block",idBlock);

  return true;
}


bool VTF::writeSblk (const std::vector<int>& sBlockIDs, const char* resultName,
		     int idBlock, int iStep, bool elementData)
{
  if ((int)mySBlock.size() < idBlock) mySBlock.resize(idBlock,0);

  int status = 1;
#if HAS_VTFAPI == 1
  if (!mySBlock[--idBlock])
  {
    mySBlock[idBlock] = new VTFAScalarBlock(idBlock+1);
    if (resultName) mySBlock[idBlock]->SetName(resultName);
    status = mySBlock[idBlock]->SetResultBlocks(&sBlockIDs.front(),
                                                sBlockIDs.size(),iStep);
  }
  else for (size_t i = 0; i < sBlockIDs.size() && VTFA_SUCCESS(status); i++)
    status = mySBlock[idBlock]->AddResultBlock(sBlockIDs[i],iStep);
#elif HAS_VTFAPI == 2
  if (!mySBlock[--idBlock])
  {
    int type = elementData?VTFXA_RESMAP_ELEMENT:VTFXA_RESMAP_NODE;
    mySBlock[idBlock] = new VTFXAResultBlock(idBlock+1,
					     VTFXA_RESTYPE_SCALAR,type);
    if (resultName) mySBlock[idBlock]->SetName(resultName);
  }
  mySBlock[idBlock]->SetResultID(idBlock);
  status = mySBlock[idBlock]->SetResultValuesBlocks(&sBlockIDs.front(),
                                                    sBlockIDs.size(),iStep);
#endif
  if (VTFA_FAILURE(status))
    return showError("Error defining scalar block",idBlock);

  return true;
}


bool VTF::writeState (int iStep, const char* fmt, Real refValue, int refType)
{

  char stepName[32];
  sprintf(stepName,fmt,refValue);
#if HAS_VTFAPI == 1
  if (!myState) myState = new VTFAStateInfoBlock();
  if (VTFA_FAILURE(myState->SetStepData(iStep,stepName,refValue,refType)))
    return showError("Error defining state info block");
#elif HAS_VTFAPI == 2
  if (!myState) myState = new VTFXAStateInfoBlock();
  if (VTFA_FAILURE(myState->AddStateInfo(iStep,stepName,refValue,
					 VTFXA_REFVALUETYPE_TIME)))
    return showError("Error defining state info block");
#endif

  return true;
}


bool VTF::writeNodes (int iBlockID)
{
  bool ok = true;

#if HAS_VTFAPI == 1
  VTFANodeBlock nBlock(iBlockID,0);
#elif HAS_VTFAPI == 2
  VTFXANodeBlock nBlock(iBlockID,false);
#endif

#ifdef HAS_VTFAPI
  const ElementBlock* grid = myBlocks.back().second;
  if (VTFA_FAILURE(nBlock.SetNumNodes(grid->getNoNodes())))
    ok = false;

  std::vector<Vec3>::const_iterator cit;
  for (cit = grid->begin_XYZ(); cit != grid->end_XYZ() && ok; cit++)
    if (VTFA_FAILURE(nBlock.AddNode(cit->x, cit->y, cit->z))) ok = false;
#endif

#if HAS_VTFAPI == 1
  if (VTFA_FAILURE(myFile->WriteBlock(&nBlock))) ok = false;
#elif HAS_VTFAPI == 2
  if (VTFA_FAILURE(myDatabase->WriteBlock(&nBlock))) ok = false;
#endif

  return ok;
}


bool VTF::writeElements (const char* partName, int iBlockID, int iNodeBlockID)
{
  bool ok = true;

#if HAS_VTFAPI == 1
  VTFAElementBlock eBlock(iBlockID,0,0);

  const ElementBlock* grid = myBlocks.back().second;
  const int* mnpc = grid->getElements();
  int nel = grid->getNoElms();
  switch (grid->getNoElmNodes()) {
  case 2:
    ok = VTFA_SUCCESS(eBlock.AddElements(VTFA_BEAMS,mnpc,nel));
    break;
  case 4:
    ok = VTFA_SUCCESS(eBlock.AddElements(VTFA_QUADS,mnpc,nel));
    break;
  case 8:
    ok = VTFA_SUCCESS(eBlock.AddElements(VTFA_HEXAHEDRONS,mnpc,nel));
    break;
  default:
    ok = false;
  }

  eBlock.SetPartID(iBlockID);
  eBlock.SetPartName(partName);
  eBlock.SetNodeBlockID(iNodeBlockID);

  if (VTFA_FAILURE(myFile->WriteBlock(&eBlock)))
    ok = false;

#elif HAS_VTFAPI == 2
  VTFXAElementBlock eBlock(iBlockID,0,0);

  const ElementBlock* grid = myBlocks.back().second;
  const int* mnpc = grid->getElements();
  int nel = grid->getNoElms();
  switch (grid->getNoElmNodes()) {
  case 2:
    ok = VTFA_SUCCESS(eBlock.AddElements(VTFXA_BEAMS,mnpc,nel));
    break;
  case 4:
    ok = VTFA_SUCCESS(eBlock.AddElements(VTFXA_QUADS,mnpc,nel));
    break;
  case 8:
    ok = VTFA_SUCCESS(eBlock.AddElements(VTFXA_HEXAHEDRONS,mnpc,nel));
    break;
  default:
    ok = false;
  }

  eBlock.SetNodeBlockID(iNodeBlockID);
  if (VTFA_FAILURE(myDatabase->WriteBlock(&eBlock)))
    ok = false;
#endif

  return ok;
}


bool VTF::showError (const char* msg, int ID)
{
  std::cerr <<"VTF: "<< msg;
  if (ID) std::cerr <<" "<< ID;
  std::cerr << std::endl;
  return false;
}
