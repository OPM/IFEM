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
#endif
#include "ElementBlock.h"
#include <iostream>
#include <stdio.h>


real VTF::vecOffset[3] = { 0.0, 0.0, 0.0 };


VTF::VTF (const char* filename, int type)
{
  myFile = 0;
  geoBlock = new VTFAGeometryBlock;
  myState = 0;
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
  showError("VTF export is not available in this version");
#endif
  myFile = 0;
}


VTF::~VTF ()
{
  if (!myFile) return;

  size_t i;
#if HAS_VTFAPI == 1
  if (VTFA_FAILURE(myFile->WriteBlock(geoBlock)))
    showError("Error writing Geometry Block");
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
  if (VTFA_FAILURE(myDatabase->WriteBlock(geoBlock)))
    showError("Error writing Geometry Block");
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

  VTFXACase* singleCase = new VTFXACase(myFile,"Single case",1,1);
  VTFXACasePropertiesBlock frameGeneratorProps(VT_CT_FRAME_GENERATOR_SETTINGS);

  frameGeneratorProps.AddInt(VT_PI_FG_FEM_MODEL_IDS, 1); // for VTFx just use always "1" here
  singleCase->WritePropertiesBlock(&frameGeneratorProps);
  for (i = 0; i < myBlocks.size(); i++) {
    VTFXACasePropertiesBlock partAttr(VT_CT_PART_ATTRIBUTES);
    partAttr.SetPartID(i+1);
    // Turn on mesh
    partAttr.AddBool(VT_PB_PA_MESH, VTFXA_FALSE);
    partAttr.AddBool(VT_PB_PA_DISPLACEMENTS,      VTFXA_FALSE);

    singleCase->WritePropertiesBlock(&partAttr);
  }
  if (VTFA_FAILURE(myFile->CloseFile()))
    showError("Error closing VTF file");

  delete singleCase;
  delete myDatabase;
  delete myFile;
#endif
  delete geoBlock;
}

void VTF::writeGeometryBlocks(int iStep)
{
  if (!myBlocks.size())
    return;
  size_t i;
  std::vector<int> geomID(myBlocks.size());
  for (i = 0; i < myBlocks.size(); i++)
  {
    geomID[i] = myBlocks[i].first;
  }

#if HAS_VTFAPI == 1
  geoBlock->SetGeometryElementBlocks(&geomID.front(),geomID.size(),iStep);
#elif HAS_VTFAPI == 2
  std::cout << "yo " << iStep << std::endl;
  geoBlock->SetElementBlocksForState(&geomID.front(),geomID.size(),iStep);
#endif
}


void VTF::clearGeometryBlocks()
{ 
  for (size_t i=0;i<myBlocks.size();++i) 
    delete myBlocks[i].second;
  myBlocks.clear();
}


bool VTF::writeGrid (const ElementBlock* block, const char* partname, int nBlock)
{
  if (!myFile) return true;

  myBlocks.push_back(std::make_pair(nBlock,block));

  if (!writeNodes(nBlock))
    return showError("Error writing node block",nBlock);

  if (!writeElements(partname,nBlock,nBlock))
    return showError("Error writing element block",nBlock);

  return true;
}


bool VTF::writeVres (const std::vector<real>& nodeResult,
		     int idBlock, int geomID, size_t nvc)
{
  if (!myFile) return true;

  const size_t nnod = myBlocks[geomID-1].second->getNoNodes();
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
  dBlock.SetMapToBlockID(geomID,VTFXA_NODES);
  dBlock.SetResultValues3D(resVec,nnod);
  if (VTFA_FAILURE(myDatabase->WriteBlock(&dBlock)))
    return showError("Error writing result block",idBlock);
#endif
  delete[] resVec;

  return true;
}


bool VTF::writeEres (const std::vector<real>& elementResult,
		     int idBlock, int geomID)
{
  if (!myFile) return true;
  if (geomID < 1 || (size_t)geomID > myBlocks.size()) return false;

  const size_t nres = elementResult.size();
  if (nres > myBlocks[geomID-1].second->getNoElms())
    return showError("Invalid size of result array",nres);
  else if (nres < myBlocks[geomID-1].second->getNoElms())
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
  dBlock.SetMapToBlockID(geomID,VTFXA_ELEMENTS);
  dBlock.SetResultValues1D(resVec,nres);
  if (VTFA_FAILURE(myDatabase->WriteBlock(&dBlock)))
    return showError("Error writing result block",idBlock);
#endif

  delete[] resVec;
  return true;
}


bool VTF::writeNres (const std::vector<real>& nodalResult,
		     int idBlock, int geomID)
{
  if (!myFile) return true;
  if (geomID < 1 || (size_t)geomID > myBlocks.size()) return false;

  const size_t nres = nodalResult.size();
  if (nres != myBlocks[geomID-1].second->getNoNodes())
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
  dBlock.SetMapToBlockID(geomID,VTFXA_NODES);
  dBlock.SetResultValues1D(resVec,nres);
  if (VTFA_FAILURE(myDatabase->WriteBlock(&dBlock)))
    return showError("Error writing result block",idBlock);
#endif
  delete[] resVec;

  return true;
}


bool VTF::writeVectors (const std::map<Vec3,Vec3>& pntResult, int idBlock)
{
#if HAS_VTFAPI == 1
  bool writePoints = false;
  static int geomID = 0;
  if (geomID == 0)
  {
    // The big assumption here is that we have only one call to writeVectors
    // per time step, and that all subsequent calls are with the same points
    myBlocks.push_back(std::make_pair(myBlocks.size()+1,new ElementBlock()));
    geomID = myBlocks.size();
    writePoints = true;
  }

  VTFANodeBlock   nBlock(geomID,0);
  VTFAResultBlock rBlock(idBlock,VTFA_DIM_VECTOR,VTFA_RESMAP_NODE,0);

  size_t i = 0, np = pntResult.size();
  if (writePoints && VTFA_FAILURE(nBlock.SetNumNodes(np)))
    return showError("Error defining node block",geomID);
  else if (VTFA_FAILURE(rBlock.SetNumResults(np)))
    return showError("Error defining result block",idBlock);
  else
    rBlock.SetMapToBlockID(geomID);

  int* mnpc = writePoints ? new int[np] : 0;
  std::map<Vec3,Vec3>::const_iterator cit;
  for (cit = pntResult.begin(); cit != pntResult.end(); cit++, i++)
    if (writePoints && VTFA_FAILURE(nBlock.AddNode(vecOffset[0]+cit->first.x,
						   vecOffset[1]+cit->first.y,
						   vecOffset[2]+cit->first.z)))
      return showError("Error adding node to block",geomID);
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
    VTFAElementBlock eBlock(geomID,0,0);
    eBlock.SetPartID(geomID);
    eBlock.SetNodeBlockID(geomID);
    if (VTFA_FAILURE(eBlock.AddElements(VTFA_POINTS,mnpc,np)))
      return showError("Error defining element block",geomID);
    delete[] mnpc;

    if (VTFA_FAILURE(myFile->WriteBlock(&nBlock)))
      return showError("Error writing node block",geomID);
    else if (VTFA_FAILURE(myFile->WriteBlock(&eBlock)))
      return showError("Error writing element block",geomID);
  }
  if (VTFA_FAILURE(myFile->WriteBlock(&rBlock)))
    return showError("Error writing result block",idBlock);
#elif HAS_VTFAPI == 2
  showError("Note: Vector points are not yet implemented for VTFx");
#endif

  return true;
}


bool VTF::writeDblk (const std::vector<int>& dBlockIDs, const char* resultName,
		     int idBlock, int iStep)
{
  if ((int)myDBlock.size() < idBlock) myDBlock.resize(idBlock,0);

#if HAS_VTFAPI == 1
  if (!myDBlock[--idBlock])
  {
    myDBlock[idBlock] = new VTFADisplacementBlock(idBlock+1);
    if (resultName) myDBlock[idBlock]->SetName(resultName);
    myDBlock[idBlock]->SetRelativeDisplacementResults(1);
  }
  if (VTFA_FAILURE(myDBlock[idBlock]->SetResultBlocks(&dBlockIDs.front(),
						      dBlockIDs.size(),iStep)))
    return showError("Error defining displacement block",idBlock);
#elif HAS_VTFAPI == 2
  if (!myDBlock[--idBlock])
  {
    myDBlock[idBlock] = new VTFXAResultBlock(idBlock+1,
					     VTFXA_RESTYPE_DISPLACEMENT,
					     VTFXA_RESMAP_NODE);
    if (resultName) myDBlock[idBlock]->SetName(resultName);
  }
  myDBlock[idBlock]->SetResultID(idBlock);
  if (VTFA_FAILURE(myDBlock[idBlock]->SetResultValuesBlocks(&dBlockIDs.front(),
						      dBlockIDs.size(),iStep)))
    return showError("Error defining displacement block",idBlock);
#endif

  return true;
}


bool VTF::writeVblk (int vBlockID, const char* resultName,
		     int idBlock, int iStep)
{
  if ((int)myVBlock.size() < idBlock) myVBlock.resize(idBlock,0);

#if HAS_VTFAPI == 1
  if (!myVBlock[--idBlock])
  {
    myVBlock[idBlock] = new VTFAVectorBlock(idBlock+1);
    if (resultName) myVBlock[idBlock]->SetName(resultName);
  }
  if (VTFA_FAILURE(myVBlock[idBlock]->SetResultBlocks(&vBlockID,1,iStep)))
    return showError("Error defining vector block",idBlock);
#elif HAS_VTFAPI == 2
  if (!myVBlock[--idBlock])
  {
    myVBlock[idBlock] = new VTFXAResultBlock(idBlock+1,
					     VTFXA_RESTYPE_VECTOR,
					     VTFXA_RESMAP_NODE);
    if (resultName) myVBlock[idBlock]->SetName(resultName);
  }
  myVBlock[idBlock]->SetResultID(idBlock);
  if (VTFA_FAILURE(myVBlock[idBlock]->SetResultValuesBlocks(&vBlockID,1,iStep)))
    return showError("Error defining vector block",idBlock);
#endif

  return true;
}


bool VTF::writeVblk (const std::vector<int>& vBlockIDs, const char* resultName,
		     int idBlock, int iStep)
{
  if ((int)myVBlock.size() < idBlock) myVBlock.resize(idBlock,0);
#if HAS_VTFAPI == 1
  if (!myVBlock[--idBlock])
  {
    myVBlock[idBlock] = new VTFAVectorBlock(idBlock+1);
    if (resultName) myVBlock[idBlock]->SetName(resultName);
  }
  if (VTFA_FAILURE(myVBlock[idBlock]->SetResultBlocks(&vBlockIDs.front(),
						      vBlockIDs.size(),iStep)))
    return showError("Error defining vector block",idBlock);
#elif HAS_VTFAPI == 2
  if (!myVBlock[--idBlock])
  {
    myVBlock[idBlock] = new VTFXAResultBlock(idBlock+1,
					     VTFXA_RESTYPE_VECTOR,
					     VTFXA_RESMAP_NODE);
    if (resultName) myVBlock[idBlock]->SetName(resultName);
  }
  myVBlock[idBlock]->SetResultID(idBlock);
  if (VTFA_FAILURE(myVBlock[idBlock]->SetResultValuesBlocks(&vBlockIDs.front(),
                                                            vBlockIDs.size(),
							    iStep)))
    return showError("Error defining vector block",idBlock);
#endif

  return true;
}


bool VTF::writeSblk (int sBlockID, const char* resultName,
		     int idBlock, int iStep, bool elementData)
{
#if HAS_VTFAPI == 1
  if ((int)mySBlock.size() < idBlock) mySBlock.resize(idBlock,0);

  if (!mySBlock[--idBlock])
  {
    mySBlock[idBlock] = new VTFAScalarBlock(idBlock+1);
    if (resultName) mySBlock[idBlock]->SetName(resultName);
  }
  if (VTFA_FAILURE(mySBlock[idBlock]->SetResultBlocks(&sBlockID,1,iStep)))
    return showError("Error defining scalar block",idBlock);
#elif HAS_VTFAPI == 2
  if ((int)mySBlock.size() < idBlock) mySBlock.resize(idBlock,0);

  if (!mySBlock[--idBlock])
  {
    int type = elementData?VTFXA_RESMAP_ELEMENT:VTFXA_RESMAP_NODE;
    mySBlock[idBlock] = new VTFXAResultBlock(idBlock+1,
					     VTFXA_RESTYPE_SCALAR,type);
    if (resultName) mySBlock[idBlock]->SetName(resultName);
  }
  mySBlock[idBlock]->SetResultID(idBlock);
  if (VTFA_FAILURE(mySBlock[idBlock]->SetResultValuesBlocks(&sBlockID,1,iStep)))
    return showError("Error defining scalar block",idBlock);
#endif

  return true;
}


bool VTF::writeSblk (const std::vector<int>& sBlockIDs, const char* resultName,
                     int idBlock, int iStep, bool elementData)
{
  if ((int)mySBlock.size() < idBlock) mySBlock.resize(idBlock,0);

#if HAS_VTFAPI == 1
  if (!mySBlock[--idBlock])
  {
    mySBlock[idBlock] = new VTFAScalarBlock(idBlock+1);
    if (resultName) mySBlock[idBlock]->SetName(resultName);
  }
  if (VTFA_FAILURE(mySBlock[idBlock]->SetResultBlocks(&sBlockIDs.front(),
						      sBlockIDs.size(),iStep)))
    return showError("Error defining scalar block",idBlock);
#elif HAS_VTFAPI == 2
  if (!mySBlock[--idBlock])
  {
    int type = elementData?VTFXA_RESMAP_ELEMENT:VTFXA_RESMAP_NODE;
    mySBlock[idBlock] = new VTFXAResultBlock(idBlock+1,
                                             VTFXA_RESTYPE_SCALAR,type);
    if (resultName) mySBlock[idBlock]->SetName(resultName);
  }
  mySBlock[idBlock]->SetResultID(idBlock);
  if (VTFA_FAILURE(mySBlock[idBlock]->SetResultValuesBlocks(&sBlockIDs.front(),
							    sBlockIDs.size(),
							    iStep)))
    return showError("Error defining scalar block",idBlock);
#endif

  return true;
}


bool VTF::writeState (int iStep, const char* fmt, real refValue, int refType)
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
