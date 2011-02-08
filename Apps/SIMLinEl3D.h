// $Id: SIMLinEl3D.h,v 1.11 2010-12-18 16:23:52 kmo Exp $
//==============================================================================
//!
//! \file SIMLinEl3D.h
//!
//! \date Dec 08 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Solution driver for 3D NURBS-based linear elastic FEM analysis.
//!
//==============================================================================

#ifndef _SIM_LIN_EL_3D_H
#define _SIM_LIN_EL_3D_H

#include "SIM3D.h"
#include "SIMenums.h"


/*!
  \brief Driver class for 3D isogeometric FEM analysis of elasticity problems.
  \details The class incapsulates data and methods for solving linear elasticity
  problems using NURBS-based finite elements. It reimplements the parse method
  of the parent class.
*/

class SIMLinEl3D : public SIM3D
{
  //! \brief Struct for storage of physical material property parameters.
  struct IsoMat
  {
    double E;   //!< Young's modulus
    double nu;  //!< Poisson's ratio
    double rho; //!< Mass density
    //! \brief Default constructor.
    IsoMat() { E = nu = rho = 0.0; }
    //! \brief Constructor initializing a material instance.
    IsoMat(double Emod, double Poiss, double D) : E(Emod), nu(Poiss), rho(D) {}
  };

  typedef std::vector<IsoMat> MaterialVec; //!< Vector of material data sets

public:
  //! \brief Default constructor.
  //! \param[in] checkRHS If \e true, ensure the model is in a right-hand system
  //! \param[in] form Problem formulation option
  SIMLinEl3D(bool checkRHS = false, int form = SIM::LINEAR);
  //! \brief The destructor frees the dynamically allocated data.
  virtual ~SIMLinEl3D();

protected:
  //! \brief Parses a data section from the input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);

  //! \brief Initializes material properties for integration of interior terms.
  //! \param[in] propInd Physical property index
  virtual bool initMaterial(size_t propInd);

  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  virtual bool initNeumann(size_t propInd);

  //! \brief Returns the analytical solution function, if any.
  virtual TensorFunc* getAnaSol() const { return asol; }

private:
  MaterialVec mVec; //!< Material data
  TensorFunc* asol; //!< Analytical stress field
};

#endif
