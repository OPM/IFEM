// $Id$
//==============================================================================
//!
//! \file LinSolParams.h
//!
//! \date Jan 29 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Linear solver parameters for PETSc matrices.
//! \details Includes linear solver method, preconditioner
//! and convergence criteria.
//!
//==============================================================================

#ifndef _LIN_SOL_PARAMS_H
#define _LIN_SOL_PARAMS_H

#include "LinAlgenums.h"
#include <map>
#include <string>
#include <vector>

class TiXmlElement;


/*!
  \brief A key-value store for settings.
*/

class SettingMap
{
public:
  //! \brief Add a value to the store.
  //! \param[in] key The key
  //! \param[in] value The value
  void addValue(const std::string& key, const std::string& value);

  //! \brief Obtain a value as a string.
  //! \param[in] key The key
  std::string getStringValue(const std::string& key) const;

  //! \brief Obtain a value as an integer.
  //! \param[in] key The key
  int getIntValue(const std::string& key) const;

  //! \brief Obtain a value as an double.
  //! \param[in] key The key
  double getDoubleValue(const std::string& key) const;

  //! \brief Checks if the store holds a value for a key.
  //! \param[in] key The key
  bool hasValue(const std::string& key) const;

protected:
  std::map<std::string,std::string> values; //!< Map of key-value pairs
};


/*!
  \brief Class for linear solver parameters.
  \details It contains information about solver method, preconditioner
  and convergence criteria.
*/

class LinSolParams : public SettingMap
{
public:
  //! \brief Default constructor.
  LinSolParams(LinAlg::LinearSystemType ls = LinAlg::GENERAL_MATRIX);
  //! \brief Copy constructor.
  LinSolParams(const LinSolParams& par,
               LinAlg::LinearSystemType ls = LinAlg::GENERAL_MATRIX);

  //! \brief Read linear solver parameters from XML document.
  bool read(const TiXmlElement* elem);

  //! \brief Linear solver settings for a block of the linear system
  class BlockParams : public SettingMap
  {
  public:
    //! \brief Default constructor
    BlockParams();

    //! \brief Settings for a directional smoother
    struct DirSmoother
    {
      int order;        //!< Ordering of DOFs
      std::string type; //!< Directional smoother types
    };

    //! \brief Read settings from XML block
    //! \param[in] child XML block
    //! \param[in] prefix Prefix to add to read data
    bool read(const TiXmlElement* child, const std::string& prefix="");

    size_t basis; //!< Basis for block
    size_t comps; //!< Components from basis (1, 2, 3, 12, 13, 23, 123, ..., 0 = all)
    std::vector<DirSmoother> dirSmoother; //!< Directional smoother data
  };

  //! \brief Number of blocks in matrix system
  size_t getNoBlocks() const { return blocks.size(); }

  //! \brief Obtain settings for a given block
  const BlockParams& getBlock(size_t i) const { return blocks[i]; }

  //! \brief Returns the linear system type.
  LinAlg::LinearSystemType getLinSysType() const { return linSys; }

private:
  std::vector<BlockParams> blocks; //!< Parameters for each block
  LinAlg::LinearSystemType linSys; //!< Type of linear system matrix
};

#endif
