// $Id$
//==============================================================================
//!
//! \file ASMu1DLag.h
//!
//! \date May 26 2021
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of unstructured 1D Lagrange FE models.
//!
//==============================================================================

#ifndef _ASM_U1D_LAG_H
#define _ASM_U1D_LAG_H

#include "ASMs1DLag.h"
#include "ASMutils.h"


/*!
  \brief Driver for assembly of unstructured 1D Lagrange FE models.
  \details This class overrides the methods of its parent class such that
  it does not depend on a curve spline object for geometry discretization.
  It can therefore be used for any unstructured grid read from mesh files.
*/

class ASMu1DLag : public ASMs1DLag
{
public:
  //! \brief Default constructor.
  ASMu1DLag(unsigned char n = 1, unsigned char n_f = 1, char fType = 'm');
  //! \brief Special copy constructor for sharing of FE data.
  ASMu1DLag(const ASMu1DLag& pch, unsigned char n_f);
  //! \brief Default copy constructor copying everything.
  ASMu1DLag(const ASMu1DLag& pch);
  //! \brief Empty destructor.
  virtual ~ASMu1DLag() {}

  // Methods for model generation
  // ============================

  //! \brief Creates an instance by reading the given input stream.
  virtual bool read(std::istream& is);
  //! \brief Generates a beam finite element model for the patch.
  //! \param[in] Zaxis Vector defining a point in the local XZ-plane
  virtual bool generateOrientedFEModel(const Vec3& Zaxis);
  //! \brief Checks if this patch is empty.
  virtual bool empty() const { return nel == 0; }

  //! \brief Returns (1-based) index of a predefined node set in the patch.
  virtual int getNodeSetIdx(const std::string& setName) const;
  //! \brief Returns an indexed pre-defined node set.
  virtual const IntVec& getNodeSet(int idx) const;
  //! \brief Returns a named node set for update.
  virtual IntVec& getNodeSet(const std::string& setName);

  // Post-processing methods
  // =======================

  //! \brief Creates a line element model of this patch for visualization.
  //! \param[out] grid The generated line grid
  //! \note The number of element nodes must be set in \a grid on input.
  virtual bool tesselate(ElementBlock& grid, const int*) const;

private:
  char                      fileType; //!< Mesh file format
  std::vector<ASM::NodeSet> nodeSets; //!< Node sets for Dirichlet BCs
  std::vector<ASM::NodeSet> elemSets; //!< Element sets for properties
};

#endif
