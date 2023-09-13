// $Id$
//==============================================================================
//!
//! \file ASMsupel.h
//!
//! \date Mar 30 2021
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of general superelements.
//!
//==============================================================================

#ifndef _ASM_SUPEL_H
#define _ASM_SUPEL_H

#include "ASMbase.h"
#include "ASMutils.h"
#include "ElmMats.h"
#include "Vec3.h"


/*!
  \brief Driver for assembly of general superelements.
  \details This class contains methods for assembly of superelements
  resulting from static condensation or reduced order modeling.
*/

class ASMsupel : public ASMbase
{
public:
  //! \brief Default constructor.
  explicit ASMsupel(unsigned char n_f = 6) : ASMbase(0,3,n_f) {}
  //! \brief Special copy constructor for sharing of FE data.
  ASMsupel(const ASMsupel& patch, unsigned char n_f) : ASMbase(patch,n_f) {}
  //! \brief Default copy constructor copying everything.
  ASMsupel(const ASMsupel& patch) : ASMbase(patch) {}
  //! \brief Empty destructor.
  virtual ~ASMsupel() {}

  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream& is);
  //! \brief Writes the geometry of the patch to the given stream.
  virtual bool write(std::ostream&, int) const;
  //! \brief Generates the finite element topology data for this patch.
  virtual bool generateFEMTopology();
  //! \brief Checks if this patch is empty.
  virtual bool empty() const { return myElmMat.empty(); }

  //! \brief Returns (1-based) index of a predefined node set in the patch.
  virtual int getNodeSetIdx(const std::string& setName) const;
  //! \brief Returns an indexed pre-defined node set.
  virtual const IntVec& getNodeSet(int idx) const;
  //! \brief Returns a named node set for update.
  virtual IntVec& getNodeSet(const std::string& setName, int& idx);

  //! \brief Returns the global coordinates for the given node.
  //! \param[in] inod 1-based node index local to current patch.
  //! If \a inod is 0, the centroid of the patch is returned.
  virtual Vec3 getCoord(size_t inod) const;
  //! \brief Returns a matrix with all nodal coordinates within the patch.
  //! \param[out] X nsd\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in the patch
  virtual void getNodalCoordinates(Matrix& X, bool = false) const;
  //! \brief Returns a matrix with nodal coordinates for an element.
  //! \param[in] iel 1-based element index local to current patch
  //! \param[out] X 3\f$\times\f$n-matrix, where \a n is the number of nodes
  //! in one element
  virtual bool getElementCoordinates(Matrix& X, int iel, bool = false) const;

  //! \brief Finds the global (or patch-local) node numbers on a patch boundary.
  //! \param[in] lIndex Local index of the boundary face/edge
  //! \param nodes Array of node numbers
  //! \param[in] local If \e true, return patch-local numbers
  virtual void getBoundaryNodes(int lIndex, IntVec& nodes,
                                int, int, int, bool local) const;

  //! \brief Dummy method doing nothing.
  virtual void getBoundaryElms(int, int, IntVec&) const {}
  //! \brief Dummy method doing nothing.
  virtual bool getParameterDomain(Real2DMat&, IntVec*) const { return false; }
  //! \brief Dummy method doing nothing.
  virtual void getElmConnectivities(IntMat&) const {}
  //! \brief Dummy method doing nothing.
  virtual bool updateCoords(const Vector&) { return false; }

  //! \brief Evaluates an integral over the interior patch domain.
  //! \param integrand Object with problem-specific data and methods
  //! \param glbInt The integrated quantity
  virtual bool integrate(Integrand& integrand, GlobalIntegral& glbInt,
                         const TimeDomain&);

  //! \brief Dummy method (patch patch boundaries are not defined).
  virtual bool integrate(Integrand&, int, GlobalIntegral&,
                         const TimeDomain&) { return false; }

  //! \brief Dummy method doing nothing.
  virtual int evalPoint(const double*, double*, Vec3&) const { return 0; }
  //! \brief Dummy method doing nothing.
  virtual int findElementContaining(const double*) const { return 0; }

  //! \brief Creates a standard FE model of this patch for visualization.
  //! \param[out] grid The generated finite element grid
  virtual bool tesselate(ElementBlock& grid, const int*) const;

  //! \brief Evaluates the primary solution field at all visualization points.
  //! \param[out] sField Solution field
  //! \param[in] locSol Solution vector local to current patch
  virtual bool evalSolution(Matrix& sField, const Vector& locSol,
                            const int*, int) const;

  //! \brief Applies a transformation matrix from local to global system.
  virtual bool transform(const Matrix& Tlg);

private:
  Vec3Vec myNodes;  //!< Supernode coordinates
  ElmMats myElmMat; //!< Duperelement matrices

  std::vector<ASM::NodeSet> nodeSets; //!< Node sets for Dirichlet BCs
};

#endif
