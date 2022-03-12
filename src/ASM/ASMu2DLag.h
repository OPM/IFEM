// $Id$
//==============================================================================
//!
//! \file ASMu2DLag.h
//!
//! \date Aug 12 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of unstructured 2D Lagrange FE models.
//!
//==============================================================================

#ifndef _ASM_U2D_LAG_H
#define _ASM_U2D_LAG_H

#include "ASMs2DLag.h"
#include "ASMutils.h"


/*!
  \brief Driver for assembly of unstructured 2D Lagrange FE models.
  \details This class overrides the methods of its parent class such that
  it does not depend on a surface spline object for geometry discretization.
  It can therefore be used for any unstructured grid read from mesh files.
*/

class ASMu2DLag : public ASMs2DLag
{
public:
  //! \brief Default constructor.
  ASMu2DLag(unsigned char n = 2, unsigned char n_f = 2, char fType = 'm');
  //! \brief Special copy constructor for sharing of FE data.
  ASMu2DLag(const ASMu2DLag& pch, unsigned char n_f);
  //! \brief Default copy constructor copying everything.
  ASMu2DLag(const ASMu2DLag& pch);
  //! \brief Empty destructor.
  virtual ~ASMu2DLag() {}

  // Methods for model generation
  // ============================

  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream& is, int = 0);
  //! \brief Generates the finite element topology data for the patch.
  virtual bool generateFEMTopology();
  //! \brief Checks if this patch is empty.
  virtual bool empty() const { return nel == 0; }

  //! \brief Returns (1-based) index of a predefined node set in the patch.
  virtual int getNodeSetIdx(const std::string& setName) const;
  //! \brief Returns an indexed pre-defined node set.
  virtual const IntVec& getNodeSet(int idx) const;
  //! \brief Returns a named node set for update.
  virtual IntVec& getNodeSet(const std::string& setName, int& idx);

  //! \brief Finds the global (or patch-local) node numbers on a patch boundary.
  //! \param[in] lIndex Local index of the boundary node set
  //! \param nodes Array of node numbers
  //! \param[in] local If \e true, return patch-local numbers
  virtual void getBoundaryNodes(int lIndex, IntVec& nodes,
                                int, int, int, bool local) const;

  using ASMs2DLag::generateThreadGroups;
  //! \brief Generates element groups for multi-threading of interior integrals.
  virtual void generateThreadGroups(const Integrand&, bool, bool);

  // Post-processing methods
  // =======================

  //! \brief Creates a quad element model of this patch for visualization.
  //! \param[out] grid The generated quadrilateral grid
  //! \note The number of element nodes must be set in \a grid on input.
  virtual bool tesselate(ElementBlock& grid, const int*) const;

private:
  char                      fileType; //!< Mesh file format
  std::vector<ASM::NodeSet> nodeSets; //!< Node sets for Dirichlet BCs
};

#endif
