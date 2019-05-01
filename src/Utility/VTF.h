// $Id$
//==============================================================================
//!
//! \file VTF.h
//!
//! \date Dec 1 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Output of FE model and results to VTF file.
//!
//==============================================================================

#ifndef _VTF_H
#define _VTF_H

#include <vector>
#include <cstddef>

class ElementBlock;
class RealFunc;
class Tensor;
class Vec3;

typedef std::pair<int,const ElementBlock*> GridBlock; //!< Convenience type
typedef std::pair<Vec3,Vec3>               Vec3Pair;  //!< Convenience type

class VTFAFile;
class VTFAStateInfoBlock;
class VTFATransformationBlock;
class VTFADisplacementBlock;
class VTFAVectorBlock;
class VTFAScalarBlock;
class VTFAGeometryBlock;

#if HAS_VTFAPI == 2
class VTFXAFile;
class VTFXADatabase;
class VTFXAStateInfoBlock;
class VTFXAResultBlock;
class VTFXAGeometryBlock;
#define VTFAFile                VTFXAFile
#define VTFAStateInfoBlock      VTFXAStateInfoBlock
#define VTFATransformationBlock VTFXAResultBlock
#define VTFADisplacementBlock   VTFXAResultBlock
#define VTFAVectorBlock         VTFXAResultBlock
#define VTFAScalarBlock         VTFXAResultBlock
#define VTFAGeometryBlock       VTFXAGeometryBlock
#endif


/*!
  \brief Class for output of FE model and results to VTF file.
*/

class VTF
{
public:
  //! \brief The constructor opens a new VTF-file.
  //! \param[in] filename Name of the VTF-file
  //! \param[in] type     File type (0 = ASCII, 1 = BINARY)
  VTF(const char* filename, int type);
  //! \brief No copying of this class.
  VTF(const VTF&) = delete;
  //! \brief The destructor finalizes and closes the VTF-file.
  virtual ~VTF();

  //! \brief Writes the FE geometry to the VTF-file.
  //! \param[in] g The FE grid that all results written are referred to
  //! \param[in] partname Name of the geometry being written
  //! \param[in] partID Part identifier
  //! \param[in] geomID Geometry block identifier
  virtual bool writeGrid(const ElementBlock* g, const char* partname,
                         int partID, int geomID);

  //! \brief Writes a transformation matrix to the VTF-file.
  //! \param[in] X Position part of the transformation
  //! \param[in] T Orientation part of the transformation
  //! \param[in] idBlock Result block identifier
  //! \param[in] gID Geometry block identifier
  bool writeTransformation(const Vec3& X, const Tensor& T,
                           int idBlock = 1, int gID = 1);
  //! \brief Writes a block of scalar nodal results to the VTF-file.
  //! \param[in] nodeResult Vector of nodal results,
  //!            length must equal the number of nodes in the geometry block
  //! \param[in] idBlock Result block identifier
  //! \param[in] gID Geometry block identifier
  virtual bool writeNres(const std::vector<Real>& nodeResult,
                         int idBlock = 1, int gID = 1);
  //! \brief Writes a block of scalar element results to the VTF-file.
  //! \param[in] elmResult Vector of element results (one per element)
  //!            length must equal the number of elements in the geometry block
  //! \param[in] idBlock Result block identifier
  //! \param[in] gID Geometry block identifier
  virtual bool writeEres(const std::vector<Real>& elmResult,
                         int idBlock = 1, int gID = 1);
  //! \brief Writes a block of vector nodal results to the VTF-file.
  //! \param[in] nodeResult Vector of nodal results, length must be equal the
  //!            number of nodes in the geometry block times 1...5
  //! \param[in] idBlock Result block identifier
  //! \param[in] gID Geometry block identifier
  //! \param[in] nvc Number of components per node in \a nodeResult
  virtual bool writeVres(const std::vector<Real>& nodeResult,
                         int idBlock = 1, int gID = 1, size_t nvc = 0);
  //! \brief Writes a block of scalar nodal function values to the VTF-file.
  //! \param[in] f The scalar function to evaluate at the grid points
  //! \param[in] time Current time
  //! \param[in] idBlock Result block identifier
  //! \param[in] gID Geometry block identifier
  bool writeNfunc(const RealFunc& f, Real time = Real(0),
                  int idBlock = 1, int gID = 1);
  //! \brief Writes a block of point vector results to the VTF-file.
  //! \param[in] pntResult A set of result vectors with associated attack points
  //! \param[in] partID Part identifier
  //! \param gID Running geometry block identifier
  //! \param[in] idBlock Result block identifier
  //! \param[in] resultName Name of the result quantity
  //! \param[in] iStep Load/Time step identifier
  //! \param[in] iBlock Vector block identifier
  //!
  //! \details This method creates a separate geometry block consisting of the
  //! attack points of the result vectors, since they are independent of the
  //! FE geometry created by the writeGrid() method.
  virtual bool writeVectors(const std::vector<Vec3Pair>& pntResult,
                            int partID, int& gID, int idBlock = 1,
                            const char* resultName = nullptr,
                            int iStep = 0, int iBlock = 1);
  //! \brief Writes a block of points (no results) to the VTF-file.
  //! \param[in] points Vector of point coordinates
  //! \param[in] partID Part identifier
  //! \param gID Running geometry block identifier
  bool writePoints(const std::vector<Vec3>& points, int partID, int& gID);

  //! \brief Writes a scalar block definition to the VTF-file.
  //! \param[in] sBlockID The result block that makes up this scalar block
  //! \param[in] resultName Name of the result quantity
  //! \param[in] idBlock Scalar block identifier
  //! \param[in] iStep Load/Time step identifier
  //! \param[in] elementData false -> data per node, true -> data per element
  bool writeSblk(int sBlockID, const char* resultName = nullptr,
                 int idBlock = 1, int iStep = 1, bool elementData = false);
  //! \brief Writes a scalar block definition to the VTF-file.
  //! \param[in] sBlockIDs All result blocks that make up this scalar block
  //! \param[in] resultName Name of the result quantity
  //! \param[in] idBlock Scalar block identifier
  //! \param[in] iStep Load/Time step identifier
  //! \param[in] elementData false -> data per node, true -> data per element
  virtual bool writeSblk(const std::vector<int>& sBlockIDs,
                         const char* resultName = nullptr, int idBlock = 1,
                         int iStep = 1, bool elementData = false);
  //! \brief Writes a vector block definition to the VTF-file.
  //! \param[in] vBlockID The result block that makes up this vector block
  //! \param[in] resultName Name of the result quantity
  //! \param[in] idBlock Vector block identifier
  //! \param[in] iStep Load/Time step identifier
  bool writeVblk(int vBlockID, const char* resultName = nullptr,
                 int idBlock = 1, int iStep = 1);
  //! \brief Writes a vector block definition to the VTF-file.
  //! \param[in] vBlockIDs All result blocks that make up this vector block
  //! \param[in] resultName Name of the result quantity
  //! \param[in] idBlock Vector block identifier
  //! \param[in] iStep Load/Time step identifier
  virtual bool writeVblk(const std::vector<int>& vBlockIDs,
                         const char* resultName = nullptr,
                         int idBlock = 1, int iStep = 1);
  //! \brief Writes a displacement block definition to the VTF-file.
  //! \param[in] dBlockIDs All result blocks that make up the displacement block
  //! \param[in] resultName Name of the result quantity
  //! \param[in] idBlock Displacement block identifier
  //! \param[in] iStep Load/Time step identifier
  virtual bool writeDblk(const std::vector<int>& dBlockIDs,
                         const char* resultName = nullptr,
                         int idBlock = 1, int iStep = 1);
  //! \brief Writes a transformation block definition to the VTF-file.
  //! \param[in] tBlockIDs All result blocks that make the transformation block
  //! \param[in] resultName Name of the result quantity
  //! \param[in] idBlock Transformation block identifier
  //! \param[in] iStep Load/Time step identifier
  bool writeTblk(const std::vector<int>& tBlockIDs,
                 const char* resultName = nullptr,
                 int idBlock = 1, int iStep = 1);

  //! \brief Writes a state info block to the VTF-file.
  //! \param[in] iStep Load/Time step identifier
  //! \param[in] fmt Format string for step name
  //! \param[in] refValue Reference value for the step (time, frequency, etc.)
  //! \param[in] refType Reference value type (0=Time, 1=Frequency, 2=Load case)
  virtual bool writeState(int iStep, const char* fmt,
                          Real refValue, int refType = 0);

  //! \brief Returns the pointer to a geometry block.
  virtual const ElementBlock* getBlock(int geomID) const;

  //! \brief Adds the current FE geometry blocks to the description block.
  void writeGeometryBlocks(int iStep);
  //! \brief Drops the current FE geometry blocks.
  virtual void clearGeometryBlocks();

private:
  //! \brief Writes a node block to the VTF-file.
  //! \details The coordinates of the last added element block are written.
  //! \param[in] iBlockID Node block identifier
  bool writeNodes(int iBlockID);
  //! \brief Writes an element block to the VTF-file.
  //! \details The element topology of the last added element block is written.
  //! \param[in] partName Name used to identify the part in GLview
  //! \param[in] partID Part identifier
  //! \param[in] iBlockID Element block identifier
  //! \param[in] iNodeBlockID Node block which the element block refers to
  bool writeElements(const char* partName, int partID,
                     int iBlockID, int iNodeBlockID);
  //! \brief Prints an error message to \a std::cerr.
  //! \param[in] msg The message to print
  //! \param[in] ID If non-negative, the value is appended to the message
  //! \return \e false (always)
  static bool showError(const char* msg, int ID = -1);

public:
  static Real vecOffset[3]; //!< Optional offset for vector attack points

private:
  VTFAFile* myFile; //!< Pointer to the actual VTF-file being written
#if HAS_VTFAPI == 2
  VTFXADatabase* myDatabase; //!< Pointer to VTFx database object for this file
#endif
  VTFAStateInfoBlock* myState; //!< The state info block for this file
  VTFAGeometryBlock* myGBlock; //!< The geometry description block for this file
  std::vector<VTFATransformationBlock*> myTBlock; //!< Transformation blocks
  std::vector<VTFADisplacementBlock*>   myDBlock; //!< Displacement blocks
  std::vector<VTFAVectorBlock*>         myVBlock; //!< Vector field blocks
  std::vector<VTFAScalarBlock*>         mySBlock; //!< Scalar field blocks

  int pointGeoID;                  //!< ID of point vector geometry block
  int lastStep;                    //!< ID of the last state written to file
  std::vector<GridBlock> myBlocks; //!< The FE geometry of the whole model
};

#endif
