// $Id$
//==============================================================================
//!
//! \file PETScSolParams.h
//!
//! \date Mar 10 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Linear solver parameters for PETSc matrices.
//! \details Includes linear solver method, preconditioner
//! and convergence criteria.
//!
//==============================================================================

#ifndef _PETSC_SOL_PARAMS_H
#define _PETSC_SOL_PARAMS_H

#include "LinSolParams.h"
#include "PETScSupport.h"

#include <set>
#include <string>
#include <vector>

class ProcessAdm;
class SettingMap;

typedef std::vector<int>         IntVec;       //!< Integer vector
typedef std::vector<IntVec>      IntMat;       //!< Integer matrix
typedef std::vector<std::string> StringVec;    //!< String vector
typedef std::vector<StringVec>   StringMat;    //!< String matrix
typedef std::vector<IS>          ISVec;        //!< Index set vector
typedef std::vector<ISVec>       ISMat;        //!< Index set matrix

//! \brief Schur preconditioner methods
enum SchurPrec { SIMPLE, MSIMPLER, PCD };


/*!
  \brief Class for PETSc solver parameters.
  \details It contains information about solver method, preconditioner
  and convergence criteria.
*/

class PETScSolParams
{
public:
  //! \brief Default constructor.
  //! \param[in] spar The base linear solver parameters
  //! \param[in] padm The process administrator
  PETScSolParams(const LinSolParams& spar, const ProcessAdm& padm) :
    params(spar), adm(padm)
  {}

  //! \brief Destructor.
  ~PETScSolParams()
  {
  }

  //! \brief Set up preconditioner parameters for PC object
  //! \param pc Preconditioner to configure
  //! \param block Block this preconditioner applies to
  //! \param prefix PETsc param prefix for block
  //! \param blockEqs The local equations belonging to block
  void setupPC(PC& pc, size_t block,
               const std::string& prefix,
               const std::set<int>& blockEqs);

  //! \brief Obtain number of blocks
  size_t getNoBlocks() const { return params.getNoBlocks(); }

  //! \brief Obtain settings for a given block
  const LinSolParams::BlockParams& getBlock(size_t i) const { return params.getBlock(i); }

  //! \brief Get integer setting
  int getIntValue(const std::string& key) const { return params.getIntValue(key); }

  //! \brief Get integer setting
  double getDoubleValue(const std::string& key) const { return params.getDoubleValue(key); }

  //! \brief Get string setting
  std::string getStringValue(const std::string& key) const { return params.getStringValue(key); }

  //! \brief Returns the linear system type.
  LinAlg::LinearSystemType getLinSysType() const { return params.getLinSysType(); }

protected:
  //! \brief Set directional smoother
  //! \param[in] pc The preconditioner to add smoother for
  //! \param[in] P The preconditioner matrix
  //! \param[in] iBlock The index of the block to add smoother to
  //! \param[in] dirIndexSet The index set for the smoother
  bool addDirSmoother(PC pc, const Mat& P,
                      int iBlock, const ISMat& dirIndexSet);

  //! \brief Set ML options
  //! \param[in] prefix The prefix of the block to set parameters for
  //! \param[in] map The map of settings to use
  void setMLOptions(const std::string& prefix, const SettingMap& map);

  //! \brief Set GAMG options
  //! \param[in] prefix The prefix of the block to set parameters for
  //! \param[in] map The map of settings to use
  void setGAMGOptions(const std::string& prefix, const SettingMap& map);

  //! \brief Set Hypre options
  //! \param[in] prefix The prefix of the block to set parameters for
  //! \param[in] map The settings to apply
  void setHypreOptions(const std::string& prefix, const SettingMap& map);

  //! \brief Setup the coarse solver in a multigrid
  //! \param[in] pc The preconditioner to set coarse solver for
  //! \param[in] prefix The prefix of the block to set parameters for
  //! \param[in] map The settings to apply
  void setupCoarseSolver(PC& pc, const std::string& prefix, const SettingMap& map);

  //! \brief Setup the smoothers in a multigrid
  //! \param[in] pc The preconditioner to set coarse solver for
  //! \param[in] iBlock The index of the  block to set parameters for
  //! \param[in] dirIndexSet The index set for direction smoothers
  //! \param blockEqs The local equations belonging to block
  void setupSmoothers(PC& pc, size_t iBlock,
                      const ISMat& dirIndexSet,
                      const std::set<int>& blockEqs);

  //! \brief Setup an additive Schwarz preconditioner
  //! \param pc The preconditioner to set coarse solver for
  //! \param[in] block The block the preconditioner belongs to
  //! \param[in] asmlu True to use LU subdomain solvers
  //! \param[in] smoother True if this is a smoother in multigrid
  //! \param blockEqs The local equations belonging to block
  void setupAdditiveSchwarz(PC& pc, size_t block,
                            bool asmlu, bool smoother,
                            const std::set<int>& blockEqs);

  const LinSolParams& params; //!< Reference to linear solver parameters.
  const ProcessAdm& adm;      //!< Reference to process administrator.
};

#endif
