// $Id$
//==============================================================================
//!
//! \file SAMpatch.h
//!
//! \date Dec 10 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of FE matrices into system matrices for multi-patch models.
//!
//==============================================================================

#ifndef _SAM_PATCH_H
#define _SAM_PATCH_H

#include "SAM.h"

class ASMbase;


/*!
  \brief This is a sub-class of SAM with added functionality for spline models.
  \details It contains some additional methods for initializing the SAM data
  for an FE model based on a set of ASMbase patches.
*/

class SAMpatch : public SAM
{
public:
  //! \brief Empty default constructor.
  SAMpatch() {}
  //! \brief Empty destructor.
  virtual ~SAMpatch() {}

  //! \brief Allocates the dynamic arrays and populates them with data.
  //! \param[in] patches All spline patches in the model
  //! \param[in] numNod Total number of unique nodes in the model
  //! \param[in] dTypes Nodal DOF type flags
  bool init(const std::vector<ASMbase*>& patches, int numNod,
            const std::vector<char>& dTypes);

  //! \brief Updates the multi-point constraint array \a TTCC.
  //! \param[in] prevSol Previous primary solution vector in DOF-order
  bool updateConstraintEqs(const Vector* prevSol = nullptr);

  //! \brief Returns the start iterator of the patch container.
  std::vector<ASMbase*>::const_iterator begin() const { return model.begin(); }
  //! \brief Returns the end iterator of the patch container.
  std::vector<ASMbase*>::const_iterator end() const { return model.end(); }
  //! \brief Returns the number of patches in model.
  size_t getNoPatches() const { return model.size(); }

  //! \brief Merges the assembly data from another %SIM with this.
  virtual bool merge(const SAM* other, const std::map<int,int>* old2new);

protected:
  //! \brief Initializes the nodal arrays \a MINEX, \a MADOF and \a MSC.
  bool initNodeDofs(const std::vector<char>& dTypes);
  //! \brief Initializes the element topology arrays \a MPMNPC and \a MMNPC.
  bool initElementConn();
  //! \brief Initializes the multi-point constraint arrays
  //! \a MPMCEQ, \a MMCEQ and \a TTCC.
  bool initConstraintEqs();

private:
  std::vector<ASMbase*> model; //!< The spline patches of the model
};

#endif
