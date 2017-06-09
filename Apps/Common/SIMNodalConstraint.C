// $Id$
//==============================================================================
//!
//! \file SIMNodalConstraint.C
//!
//! \date Nov 4 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Base class for simulators constraining a topologyset to a given node.
//!
//==============================================================================

#include "SIMNodalConstraint.h"
#include "ASMs1D.h"
#include "ASMs2D.h"
#include "ASMs3D.h"
#ifdef HAS_LRSPLINE
#include "ASMu2D.h"
#include "LRSpline/LRSpline.h"
#endif


//! \brief Helper for applying operation to different ASM types
class NodalConstraintASMHelper {
public:
  //! \brief Constructor
  //! \param pch The associated ASM class
  NodalConstraintASMHelper(ASMbase* pch) : bpch(pch) {}

  //! \brief Empty destructor
  virtual ~NodalConstraintASMHelper() {}

  //! \brief Obtain the global node number of a given corner node on patch.
  //! \param vertex Vertex to obtain
  //! \param basis Basis for vertex
  virtual int getCorner(int vertex, int basis) = 0;

  //! \brief Constrain a given edge to a given node.
  //! \param item Edge index on patch
  //! \param comp Component to constrain
  //! \param basis Basis to constrain edge for
  //! \param idx Global node to constrain edge to.
  virtual void constrainEdge(int item, int comp, int basis, int idx) = 0;

  //! \brief Constrain a given vertex to a given node.
  //! \param item item Vertex index on patch.
  //! \param comp Component to constrain
  //! \param basis Basis to constrain vertex for.
  //! \param idx Global node to constrain edge to.
  void constrainVertex(int item, int comp, int basis, int idx)
  {
    int gn = bpch->getNodeID(getCorner(item, basis));
    if (gn != idx)
      bpch->add2PC(gn, comp, idx);
  }

  //! \brief Constrain the patch to a given node.
  //! \param[in] comp Component to constrain
  //! \param[in] basis Basis to constrain vertex for.
  //! \param[in] idx Global node to constrain patch to.
  void constrainPatch(int comp, int basis, int idx)
  {
    size_t ofs = getStartNode(basis);
    for (size_t i = 1; i <= bpch->getNoNodes(basis); ++i) {
      int gn = bpch->getNodeID(ofs+i);
      if (gn != idx)
        bpch->add2PC(bpch->getNodeID(ofs+i), comp, idx);
    }
  }

  //! \brief Obtain the starting node for a given basis.
  //! \param[in] basis Basis to obtain the starting node for
  size_t getStartNode(size_t basis)
  {
    size_t ofs = 0;
    for (size_t i=1;i<basis;++i)
      ofs += bpch->getNoNodes(i);

    return ofs;
  }

protected:
  ASMbase* bpch; //!< ASMbase pointer to associated patch
};


//! \brief Helper for apply constraints to a structured 1D model.
class NodalConstraintASMs1DHelper : public NodalConstraintASMHelper {
public:
  //! \brief Constructor
  //! \param pch The associated ASM class
  NodalConstraintASMs1DHelper(ASMs1D* spch) :
    NodalConstraintASMHelper(spch), pch(spch) {}

  //! \brief Obtain the global node number of a given corner node on patch.
  //! \param vertex Vertex to obtain
  //! \param basis Basis for vertex
  int getCorner(int vertex, int basis)
  {
    size_t ofs = getStartNode(basis);
    return pch->getNodeID(ofs+(vertex==1?1:pch->getSize(basis)));
  }

  //! \brief Constrain a given edge to a given node.
  //! \param item Edge index on patch
  //! \param comp Component to constrain
  //! \param basis Basis to constrain edge for
  //! \param idx Global node to constrain edge to.
  void constrainEdge(int item, int comp, int basis, int idx) {}
protected:
  ASMs1D* pch; //!< The associated patch.
};


//! \brief Helper for apply constraints to a structured 2D model.
class NodalConstraintASMs2DHelper : public NodalConstraintASMHelper {
public:
  //! \brief Constructor
  //! \param pch The associated ASM class
  NodalConstraintASMs2DHelper(ASMs2D* spch) :
    NodalConstraintASMHelper(spch), pch(spch) {}

  //! \copydoc NodalConstrainASM2DHelper::getCorner
  int getCorner(int vertex, int basis)
  {
    int n1, n2;
    pch->getSize(n1, n2, basis);
    size_t ofs = getStartNode(basis);
    const std::vector<int> idxs = {1, n1, n1*(n2-1)+1, n1*n2};
    return pch->getNodeID(idxs[vertex-1]+ofs);
  }

  //! \copydoc NodalConstrainASM2DHelper::constrainEdge
  void constrainEdge(int item, int comp, int basis, int idx)
  {
    size_t ofs = getStartNode(basis);

    int n1, n2, node = 1;
    pch->getSize(n1,n2,basis);

    switch (item) {
      case  2: // Right edge (positive I-direction)
        node += n1-1;
      case 1: // Left edge (negative I-direction)
        for (int i2 = 1; i2 <= n2; i2++, node += n1) {
          int gn = pch->getNodeID(ofs+node);
          if (gn != idx)
            pch->add2PC(pch->getNodeID(ofs+node), comp, idx);
        }
        break;

      case  4: // Back edge (positive J-direction)
        node += n1*(n2-1);
      case 3: // Front edge (negative J-direction)
        for (int i1 = 1; i1 <= n1; i1++, node++) {
          int gn = pch->getNodeID(ofs+node);
          if (gn != idx)
            pch->add2PC(gn, comp, idx);
        }
      default:
        break;
    }
  }
protected:
  ASMs2D* pch; //!< The associated patch.
};


//! \brief Helper for apply constraints to a structured 3D model.
class NodalConstraintASMs3DHelper : public NodalConstraintASMHelper {
public:
  //! \brief Constructor
  //! \param spch The associated patch
  NodalConstraintASMs3DHelper(ASMs3D* spch) :
    NodalConstraintASMHelper(spch), pch(spch) {}

  //! \copydoc NodalConstraintASMHelper::getCorner
  int getCorner(int vertex, int basis)
  {
    size_t ofs = getStartNode(basis);
    int n1, n2, n3;
    pch->getSize(n1, n2, n3, basis);
    int ofs_j = n1*(n2-1);
    int ofs_k = n1*n2*(n3-1);
    const std::vector<int> idxs = {             1,             n1,
                                          ofs_j+1,       ofs_j+n1,
                                          ofs_k+1,       ofs_k+n1,
                                    ofs_k+ofs_j+1, ofs_k+ofs_j+n1};
    return pch->getNodeID(idxs[vertex-1]+ofs);
  }

  //! \copydoc NodalConstrainASMHelper::constrainEdge
  void constrainEdge(int item, int comp, int basis, int idx)
  {
    size_t node = getStartNode(basis)+1;

    int n1, n2, n3;
    pch->getSize(n1,n2,n3,basis);

    size_t inc = 1;
    int n;
    if (item > 8) {
      inc = n1*n2;
      n = n3;
    } else if (item > 4) {
      inc = n1;
      n = n2;
    } else
      n = n1;

    switch (item)
    {
      case  6:
      case 10:
        node += n1 - 1;
        break;
      case  2:
      case 11:
        node += n1*(n2-1);
        break;
      case 12:
        node += n1*n2 - 1;
        break;
      case  3:
      case  7:
        node += n1*n2*(n3-1);
        break;
      case  8:
        node += n1*(n2*(n3-1) + 1) - 1;
        break;
      case  4:
        node += n1*(n2*n3-1);
        break;
    }

    for (int i = 1; i <= n; i++, node += inc) {
      int gn = pch->getNodeID(node);
      if (gn != idx)
        pch->add2PC(gn, comp, idx);
    }
  }

  //! \brief Constrain a given face to a given node.
  //! \param item Face index on patch
  //! \param comp Component to constrain
  //! \param basis Basis to constrain edge for
  //! \param idx Global node to constrain edge to.
  void constrainFace(int item, int comp, int basis, int idx)
  {
    int node = getStartNode(basis)+1;
    int n1, n2, n3;
    pch->getSize(n1,n2,n3,basis);

    const std::vector<int> faceDir = {-1, 1, -2, 2, -3, 3};
    switch (faceDir[item-1])
    {
      case  1: // Right face (positive I-direction)
        node += n1-1;
      case -1: // Left face (negative I-direction)
        for (int i3 = 1; i3 <= n3; i3++)
          for (int i2 = 1; i2 <= n2; i2++, node += n1) {
            int gn = pch->getNodeID(node);
            if (gn != idx)
              pch->add2PC(gn, comp, idx);
          }
        break;

      case  2: // Back face (positive J-direction)
        node += n1*(n2-1);
      case -2: // Front face (negative J-direction)
        for (int i3 = 1; i3 <= n3; i3++, node += n1*(n2-1))
          for (int i1 = 1; i1 <= n1; i1++, node++) {
            int gn = pch->getNodeID(node);
            if (gn != idx)
              pch->add2PC(gn, comp, idx);
          }
        break;

      case  3: // Top face (positive K-direction)
        node += n1*n2*(n3-1);
      case -3: // Bottom face (negative K-direction)
        for (int i2 = 1; i2 <= n2; i2++)
          for (int i1 = 1; i1 <= n1; i1++, node++) {
            int gn = pch->getNodeID(node);
            if (gn != idx)
              pch->add2PC(gn, comp, idx);
          }
        break;
    }
  }
protected:
  ASMs3D* pch; //!< The associated patch.
};


#ifdef HAS_LRSPLINE
class NodalConstraintASMu2DHelper : public NodalConstraintASMHelper {
public:
  //! \brief Constructor
  //! \param upch Associated patch
  NodalConstraintASMu2DHelper(ASMu2D* upch) :
    NodalConstraintASMHelper(upch), pch(upch) {}

  //! \copydoc NodalConstraintASM2DHelper::getCorner
  int getCorner(int vertex, int basis)
  {
    static const int indices[4][2] = {{-1,-1}, {1, -1}, {-1, 1}, {1,1}};
    return pch->getCorner(indices[vertex-1][0], indices[vertex-1][1], basis);
  }

  //! \copydoc NodalConstraintASM2DHelper::constrainEdge
  void constrainEdge(int item, int comp, int basis, int idx)
  {
    std::vector<int> map = { LR::WEST, LR::EAST, LR::SOUTH, LR::NORTH };
    std::vector<int> nodes = pch->getEdgeNodes(map[item-1], basis, 0);
    for (auto& it : nodes) {
      int gn = pch->getNodeID(it);
      if (gn != idx)
        pch->add2PC(gn, comp, idx);
    }
  }
protected:
  ASMu2D* pch; //!< The associated patch.
};
#endif


/*!
  \brief Helper function to create a NodalConstraintASMHelper instance.
*/

static NodalConstraintASMHelper* get2DHelper(ASMbase* pch)
{
  ASMs2D* spch = dynamic_cast<ASMs2D*>(pch);
  if (spch)
    return new NodalConstraintASMs2DHelper(spch);

#ifdef HAS_LRSPLINE
  ASMu2D* upch = dynamic_cast<ASMu2D*>(pch);
  if (upch)
    return new NodalConstraintASMu2DHelper(upch);
#endif

  return nullptr;
}


template<> bool SIMNodalConstraint<SIM1D>::applyConstraint()
{
  for (const auto& it3 : vertConstraints) {
    TopologySet::const_iterator it = SIM1D::myEntitys.find(it3.topset);
    if (it != SIM1D::myEntitys.end()) {
      ASMs1D* pch = static_cast<ASMs1D*>(this->getPatch(it3.patch));
      if (!pch)
        continue;
      NodalConstraintASMs1DHelper helper(pch);
      int idx = pch->getNodeID(helper.getCorner(it3.vertex, it3.basis));
      for (const auto& it2 : it->second) {
        ASMs1D* pch2 = static_cast<ASMs1D*>(this->getPatch(it2.patch));
        if (!pch2)
          continue;
        NodalConstraintASMs1DHelper helper2(pch2);
        if (it2.idim == 1)
          helper2.constrainPatch(it3.comp, it3.basis, idx);
        else if (it2.idim == 0) // vertex constraints
          helper.constrainVertex(it2.item, it3.comp, it3.basis, idx);
      }
    }
  }
  return true;
}


template<> bool SIMNodalConstraint<SIM2D>::applyConstraint()
{
  for (const auto& it3 : vertConstraints) {
    TopologySet::const_iterator it = SIM2D::myEntitys.find(it3.topset);
    if (it != SIM2D::myEntitys.end()) {
      std::unique_ptr<NodalConstraintASMHelper> helper(get2DHelper(this->getPatch(it3.patch)));
      int idx = helper->getCorner(it3.vertex, it3.basis);
      for (const auto& it2 : it->second) {
        std::unique_ptr<NodalConstraintASMHelper> helper2(get2DHelper(this->getPatch(it2.patch)));
        if (it2.idim == 2)
          helper2->constrainPatch(it3.comp, it3.basis, idx);
        else if (it2.idim == 1) // Edge constraints
          helper2->constrainEdge(it2.item, it3.comp, it3.basis, idx);
        else if (it2.idim == 0) // Vertex constraint
          helper2->constrainVertex(it2.item, it3.comp, it3.basis, idx);
      }
    }
  }
  return true;
}


template<> bool SIMNodalConstraint<SIM3D>::applyConstraint()
{
  for (const auto& it3 : vertConstraints) {
    TopologySet::const_iterator it = SIM3D::myEntitys.find(it3.topset);
    if (it != SIM3D::myEntitys.end()) {
      ASMs3D* pch = static_cast<ASMs3D*>(this->getPatch(it3.patch));
      if (!pch)
        continue;
      NodalConstraintASMs3DHelper helper(pch);
      int idx = helper.getCorner(it3.vertex, it3.basis);
      for (const auto& it2 : it->second) {
        ASMs3D* pch2 = static_cast<ASMs3D*>(this->getPatch(it2.patch));
        if (!pch2)
          continue;
        NodalConstraintASMs3DHelper helper2(pch2);
        if (it2.idim == 3)
          helper2.constrainPatch(it3.comp, it3.basis, idx);
        else if (it2.idim == 2) // Face constraints
          helper2.constrainFace(it2.item, it3.comp, it3.basis, idx);
        else if (it2.idim == 1) // Edge constraints
          helper2.constrainEdge(it2.item, it3.comp, it3.basis, idx);
        else if (it2.idim == 0) // Vertex constraint
          helper2.constrainVertex(it2.item, it3.comp, it3.basis, idx);
      }
    }
  }
  return true;
}
