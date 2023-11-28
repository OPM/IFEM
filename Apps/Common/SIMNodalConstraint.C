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
#include "IFEM.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "Utilities.h"

#ifdef HAS_LRSPLINE
#include "ASMu2D.h"
#endif

#include <tinyxml.h>


//! \brief Base class for helpers applying nodal constraints.
class NodalConstraintASMHelper {
public:
  //! \brief The constructor initializes the patch pointer.
  explicit NodalConstraintASMHelper(ASMbase* pch) : bpch(pch) {}

  //! \brief Empty destructor.
  virtual ~NodalConstraintASMHelper() {}

  //! \brief Returns the local node number of a given corner of the patch.
  //! \param[in] vertex Vertex index to return the node number for
  //! \param[in] basis Basis for vertex
  virtual int getCorner(int vertex, int basis) const = 0;

  //! \brief Constrains a given edge to a given node.
  //! \param[in] item Edge index on patch
  //! \param[in] comp Component to constrain
  //! \param[in] basis Basis to constrain edge for
  //! \param[in] idx Global node to constrain edge to
  virtual void constrainEdge(int item, int comp, int basis, int idx) {}

  //! \brief Constrains a given face to a given node.
  //! \param[in] item Face index on patch
  //! \param[in] comp Component to constrain
  //! \param[in] basis Basis to constrain face for
  //! \param[in] idx Global node to constrain face to
  virtual void constrainFace(int item, int comp, int basis, int idx) {}

  //! \brief Constrains a given vertex to a given node.
  //! \param[in] item Vertex index on patch
  //! \param[in] comp Component to constrain
  //! \param[in] basis Basis to constrain vertex for
  //! \param[in] idx Global node to constrain edge to
  void constrainVertex(int item, int comp, int basis, int idx)
  {
    this->constrainNode(this->getCorner(item,basis),comp,idx);
  }

  //! \brief Constrains the whole patch to a given node.
  //! \param[in] comp Component to constrain
  //! \param[in] basis Basis to constrain patch for
  //! \param[in] idx Global node to constrain patch to
  void constrainPatch(int comp, int basis, int idx)
  {
    int nnod = bpch->getNoNodes(basis);
    int node = this->getStartNode(basis);
    for (int i = 0; i < nnod; i++)
      this->constrainNode(++node,comp,idx);
  }

protected:
  //! \brief Constrains a given node to another node.
  //! \param[in] item Local index of node to constrain
  //! \param[in] comp Component to constrain
  //! \param[in] idx Global index of the node to constrain to
  void constrainNode(int item, int comp, int idx)
  {
    int node = bpch->getNodeID(item);
    if (node != idx)
      bpch->add2PC(node,comp,idx);
  }

  //! \brief Returns the starting node for the given \a basis.
  int getStartNode(int basis) const
  {
    int ofs = 0;
    for (int i = 1; i < basis; i++)
      ofs += bpch->getNoNodes(i);
    return ofs;
  }

protected:
  ASMbase* bpch; //!< Pointer to the associated patch
};


//! \brief Helper for apply constraints to a structured 1D model.
class NodalConstraintASMs1DHelper : public NodalConstraintASMHelper {
public:
  //! \brief The constructor forwards to the parent class constructor.
  explicit NodalConstraintASMs1DHelper(ASMs1D* pch)
    : NodalConstraintASMHelper(pch) {}

  //! \copydoc NodalConstraintASMHelper::getCorner
  virtual int getCorner(int vertex, int basis) const
  {
    int n1 = static_cast<ASMs1D*>(bpch)->getSize(basis);
    int ofs = this->getStartNode(basis);
    return bpch->getNodeID(ofs + (vertex == 1 ? 1 : n1));
  }
};


//! \brief Helper for apply constraints to a structured 2D model.
class NodalConstraintASMs2DHelper : public NodalConstraintASMHelper {
public:
  //! \brief The constructor forwards to the parent class constructor.
  explicit NodalConstraintASMs2DHelper(ASMs2D* pch)
    : NodalConstraintASMHelper(pch) {}

  //! \copydoc NodalConstraintASMHelper::getCorner
  virtual int getCorner(int vertex, int basis) const
  {
    int n1, n2, ofs = this->getStartNode(basis);
    static_cast<ASMs2D*>(bpch)->getSize(n1,n2,basis);
    const IntVec idxs = { 1, n1, n1*(n2-1)+1, n1*n2 };
    return bpch->getNodeID(idxs[vertex-1]+ofs);
  }

  //! \copydoc NodalConstraintASMHelper::constrainEdge
  virtual void constrainEdge(int item, int comp, int basis, int idx)
  {
    int n1, n2, node = 1 + this->getStartNode(basis);
    static_cast<ASMs2D*>(bpch)->getSize(n1,n2,basis);

    switch (item) {
      case  2: // Right edge (positive I-direction)
        node += n1-1;
      case 1: // Left edge (negative I-direction)
        for (int i2 = 1; i2 <= n2; i2++, node += n1)
          this->constrainNode(node,comp,idx);
        break;

      case  4: // Back edge (positive J-direction)
        node += n1*(n2-1);
      case 3: // Front edge (negative J-direction)
        for (int i1 = 1; i1 <= n1; i1++, node++)
          this->constrainNode(node,comp,idx);
      default:
        break;
    }
  }
};


//! \brief Helper for apply constraints to a structured 3D model.
class NodalConstraintASMs3DHelper : public NodalConstraintASMHelper {
public:
  //! \brief The constructor forwards to the parent class constructor.
  explicit NodalConstraintASMs3DHelper(ASMs3D* pch)
    : NodalConstraintASMHelper(pch) {}

  //! \copydoc NodalConstraintASMHelper::getCorner
  virtual int getCorner(int vertex, int basis) const
  {
    int n1, n2, n3, ofs = this->getStartNode(basis);
    static_cast<ASMs3D*>(bpch)->getSize(n1,n2,n3,basis);
    int ofs_j = n1*(n2-1);
    int ofs_k = n1*n2*(n3-1);

    const IntVec idxs = {             1,             n1,
                                ofs_j+1,       ofs_j+n1,
                                ofs_k+1,       ofs_k+n1,
                          ofs_k+ofs_j+1, ofs_k+ofs_j+n1};
    return bpch->getNodeID(idxs[vertex-1]+ofs);
  }

  //! \copydoc NodalConstraintASMHelper::constrainEdge
  virtual void constrainEdge(int item, int comp, int basis, int idx)
  {
    int n1, n2, n3, node = 1 + this->getStartNode(basis);
    static_cast<ASMs3D*>(bpch)->getSize(n1,n2,n3,basis);

    int inc = 1;
    int n = n1;
    if (item > 8) {
      inc = n1*n2;
      n = n3;
    } else if (item > 4) {
      inc = n1;
      n = n2;
    }

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

    for (int i = 1; i <= n; i++, node += inc)
      this->constrainNode(node,comp,idx);
  }

  //! \copydoc NodalConstraintASMHelper::constrainFace
  virtual void constrainFace(int item, int comp, int basis, int idx)
  {
    int n1, n2, n3, node = 1 + this->getStartNode(basis);
    static_cast<ASMs3D*>(bpch)->getSize(n1,n2,n3,basis);

    const std::vector<int> faceDir = {-1, 1, -2, 2, -3, 3};
    switch (faceDir[item-1])
    {
      case  1: // Right face (positive I-direction)
        node += n1-1;
      case -1: // Left face (negative I-direction)
        for (int i3 = 1; i3 <= n3; i3++)
          for (int i2 = 1; i2 <= n2; i2++, node += n1)
            this->constrainNode(node,comp,idx);
        break;

      case  2: // Back face (positive J-direction)
        node += n1*(n2-1);
      case -2: // Front face (negative J-direction)
        for (int i3 = 1; i3 <= n3; i3++, node += n1*(n2-1))
          for (int i1 = 1; i1 <= n1; i1++, node++)
            this->constrainNode(node,comp,idx);
        break;

      case  3: // Top face (positive K-direction)
        node += n1*n2*(n3-1);
      case -3: // Bottom face (negative K-direction)
        for (int i2 = 1; i2 <= n2; i2++)
          for (int i1 = 1; i1 <= n1; i1++, node++)
            this->constrainNode(node,comp,idx);
        break;
    }
  }
};


#ifdef HAS_LRSPLINE
//! \brief Helper for apply constraints to an unstructured 2D model.
class NodalConstraintASMu2DHelper : public NodalConstraintASMHelper {
public:
  //! \brief Constructor
  //! \param upch Associated patch
  explicit NodalConstraintASMu2DHelper(ASMu2D* pch)
    : NodalConstraintASMHelper(pch) {}

  //! \copydoc NodalConstraintASMHelper::getCorner
  virtual int getCorner(int vertex, int basis) const
  {
    static const int indices[4][2] = {{-1,-1}, {1,-1}, {-1,1}, {1,1}};
    return static_cast<ASMu2D*>(bpch)->getCorner(indices[vertex-1][0],
                                                 indices[vertex-1][1],
                                                 basis);
  }

  //! \copydoc NodalConstraintASMHelper::constrainEdge
  virtual void constrainEdge(int item, int comp, int basis, int idx)
  {
    IntVec nodes;
    bpch->getBoundaryNodes(item,nodes,basis,1,0,true);
    for (int node : nodes)
      this->constrainNode(node,comp,idx);
  }
};
#endif


//! \brief Template specialization for 1D.
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


//! \brief Template specialization for 2D.
template<> bool SIMNodalConstraint<SIM2D>::applyConstraint()
{
  // Lambda function to create the right NodalConstraintASMHelper instance
  auto&& getHelper = [this](int pidx) -> NodalConstraintASMHelper*
  {
    ASMbase* pch = this->getPatch(pidx);
    ASMs2D* spch = dynamic_cast<ASMs2D*>(pch);
    if (spch) return new NodalConstraintASMs2DHelper(spch);
#ifdef HAS_LRSPLINE
    ASMu2D* upch = dynamic_cast<ASMu2D*>(pch);
    if (upch) return new NodalConstraintASMu2DHelper(upch);
#endif
    return nullptr;
  };

  for (const auto& it3 : vertConstraints) {
    TopologySet::const_iterator it = SIM2D::myEntitys.find(it3.topset);
    if (it != SIM2D::myEntitys.end()) {
      std::unique_ptr<NodalConstraintASMHelper> helper(getHelper(it3.patch));
      int idx = helper->getCorner(it3.vertex, it3.basis);
      for (const auto& it2 : it->second) {
        std::unique_ptr<NodalConstraintASMHelper> helper2(getHelper(it2.patch));
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


//! \brief Template specialization for 3D.
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


template <class Dim>
bool SIMNodalConstraint<Dim>::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"constraintovertex"))
    return this->Dim::parse(elem);

  TopSetToVertex topset;
  utl::getAttribute(elem,"set",topset.topset);
  utl::getAttribute(elem,"patch",topset.patch);
  utl::getAttribute(elem,"vertex",topset.vertex);
  utl::getAttribute(elem,"comp",topset.comp);
  utl::getAttribute(elem,"basis",topset.basis);
  vertConstraints.push_back(topset);
  IFEM::cout <<"\tConstraining set \""<< topset.topset
             <<"\" to P"<< topset.patch <<" V"<< topset.vertex
             <<" in direction "<< topset.comp;
  if (topset.basis > 1)
    IFEM::cout <<" (basis "<< topset.basis <<")";
  IFEM::cout << std::endl;
  return true;
}


template class SIMNodalConstraint<SIM1D>;
template class SIMNodalConstraint<SIM2D>;
template class SIMNodalConstraint<SIM3D>;
