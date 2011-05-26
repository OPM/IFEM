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

#include <map>
#include <vector>

struct ElementBlock;
class Vec3;
class VTFAFile;
class VTFAStateInfoBlock;
class VTFADisplacementBlock;
class VTFAVectorBlock;
class VTFAScalarBlock;


/*!
  \brief Class for output of FE model and results to VTF file.
*/

class VTF
{
public:
  //! \brief Constructor which opens a new VTF-file.
  //! \param[in] filename Name of the VTF-file
  //! \param[in] type     File type (0 = ASCII, 1 = BINARY)
  VTF(const char* filename, int type = 0);
  //! \brief The destructor finalizes and closes the VTF-file.
  ~VTF();

  //! \brief Writes the FE geometry to the VTF-file.
  //! \param[in] g The FE grid that all results written are referred to
  //! \param[in] partname Name of the geometry being written
  bool writeGrid(const ElementBlock* g, const char* partname);

  //! \brief Writes a block of scalar nodal results to the VTF-file.
  //! \param[in] nodeResult Vector of nodal results,
  //!            length must equal the number of nodes in the geometry block
  //! \param[in] idBlock Result block identifier
  //! \param[in] gID Geometry block identifier
  bool writeNres(const std::vector<real>& nodeResult,
		 int idBlock = 1, int gID = 1);
  //! \brief Writes a block of scalar element results to the VTF-file.
  //! \param[in] elmResult Vector of element results (one per element)
  //!            length must equal the number of elements in the geometry block
  //! \param[in] idBlock Result block identifier
  //! \param[in] gID Geometry block identifier
  bool writeEres(const std::vector<real>& elmResult,
		 int idBlock = 1, int gID = 1);
  //! \brief Writes a block of vector nodal results to the VTF-file.
  //! \param[in] nodeResult Vector of nodal results, length must be equal the
  //!            number of nodes in the geometry block times 1...5
  //! \param[in] idBlock Result block identifier
  //! \param[in] gID Geometry block identifier
  //! \param[in] nvc Number of components per node in \a nodeResult
  bool writeVres(const std::vector<real>& nodeResult,
		 int idBlock = 1, int gID = 1, size_t nvc = 0);
  //! \brief Writes a block of point vector results to the VTF-file.
  //! \details This method creates a separate geometry block consisting of the
  //! attack points of the result vectors, since they are independent of the
  //! FE geometry created by the \a writeGrid method.
  //! \param[in] pntResult A set of result vectors with associated attack points
  //! \param[in] idBlock Vector block identifier
  bool writeVectors(const std::map<Vec3,Vec3>& pntResult, int idBlock = 1);

  //! \brief Writes a scalar block definition to the VTF-file.
  //! \param[in] sBlockID The result block that makes up this scalar block
  //! \param[in] resultName Name of the result quantity
  //! \param[in] idBlock Scalar block identifier
  //! \param[in] iStep Load/Time step identifier
  bool writeSblk(int sBlockID,
                 const char* resultName = 0, int idBlock = 1, int iStep = 1);
  //! \brief Writes a scalar block definition to the VTF-file.
  //! \param[in] sBlockIDs All result blocks that make up this scalar block
  //! \param[in] resultName Name of the result quantity
  //! \param[in] idBlock Scalar block identifier
  //! \param[in] iStep Load/Time step identifier
  bool writeSblk(const std::vector<int>& sBlockIDs,
                 const char* resultName = 0, int idBlock = 1, int iStep = 1);
  //! \brief Writes a vector block definition to the VTF-file.
  //! \param[in] vBlockID The result block that makes up this vector block
  //! \param[in] resultName Name of the result quantity
  //! \param[in] idBlock Vector block identifier
  //! \param[in] iStep Load/Time step identifier
  bool writeVblk(int vBlockID,
                 const char* resultName = 0, int idBlock = 1, int iStep = 1);
  //! \brief Writes a vector block definition to the VTF-file.
  //! \param[in] vBlockIDs All result blocks that make up this vector block
  //! \param[in] resultName Name of the result quantity
  //! \param[in] idBlock Vector block identifier
  //! \param[in] iStep Load/Time step identifier
  bool writeVblk(const std::vector<int>& vBlockIDs,
                 const char* resultName = 0, int idBlock = 1, int iStep = 1);
  //! \brief Writes a displacement block definition to the VTF-file.
  //! \param[in] dBlockIDs All result blocks that make up the displacement block
  //! \param[in] resultName Name of the result quantity
  //! \param[in] idBlock Displacement block identifier
  //! \param[in] iStep Load/Time step identifier
  bool writeDblk(const std::vector<int>& dBlockIDs,
		 const char* resultName = 0, int idBlock = 1, int iStep = 1);

  //! \brief Writes a state info block to the VTF-file.
  //! \param[in] iStep Load/Time step identifier
  //! \param[in] fmt Format string for step name
  //! \param[in] refValue Reference value for the step (time, frequency, etc.)
  //! \param[in] refType Reference value type (0=Time, 1=Frequency, 2=Load case)
  bool writeState(int iStep, const char* fmt, real refValue, int refType = 0);

  //! \brief Returns the pointer to a geometry block.
  const ElementBlock* getBlock(int geomID) const { return myBlocks[geomID-1]; }

private:
  //! \brief Writes a geometry block to the VTF-file.
  //! \param[in] parts Array of geometry part indices
  //! \param[in] nPart Number of geometry parts
  bool writeGeometry(const int* parts, int nPart);
  //! \brief Writes a node block to the VTF-file.
  //! \details The coordinates of the last added element block are written.
  //! \param[in] blockID Node block identifier
  bool writeNodes(int blockID);
  //! \brief Writes an element block to the VTF-file.
  //! \details The element topology of the last added element block is written.
  //! \param[in] partName Name used to identify the part in GLview
  //! \param[in] blockID Element block identifier
  //! \param[in] nodesID Node block which the written element block refers to
  bool writeElements(const char* partName, int blockID, int nodesID);
  //! \brief Prints an error message to \a std::cerr.
  //! \param[in] msg The message to print
  //! \param[in] ID If non-zero, the value is appended to the printed message
  //! \return \e false (always)
  static bool showError(const char* msg, int ID = 0);

public:
  static real vecOffset[3]; //!< Optional offset for vector attack points

private:
  VTFAFile* myFile; //!< Pointer to the actual VTF-file being written
  VTFAStateInfoBlock* myState; //!< The state info block for this file
  std::vector<VTFADisplacementBlock*> myDBlock; //!< Displacement blocks
  std::vector<VTFAVectorBlock*>       myVBlock; //!< Vector field blocks
  std::vector<VTFAScalarBlock*>       mySBlock; //!< Scalar field blocks
  std::vector<const ElementBlock*>    myBlocks; //!< The FE geometry
};

#endif
