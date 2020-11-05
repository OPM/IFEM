// $Id$
//==============================================================================
//!
//! \file Interface.h
//!
//! \date Aug 28 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Representation of domain interfaces.
//!
//==============================================================================

#ifndef _INTERFACE_H
#define _INTERFACE_H


namespace ASM
{
  /*!
    \brief Struct for representing a domain interface.
  */
  struct Interface
  {
    int master = 0; //!< Master patch (global number)
    int slave  = 0; //!< Slave patch (global number)
    int midx   = 0; //!< Index of boundary on master
    int sidx   = 0; //!< Index of boundary on slave
    int orient = 0; //!< Orientation flag (for 3D patches)
    int dim    = 0; //!< Dimension of boundary
    int basis  = 0; //!< Basis of boundary
    int thick  = 1; //!< Thickness of connection
  };


  /*!
    \brief Base class to check for internal boundary integrand contributions.
  */
  class InterfaceChecker
  {
  public:
    //! \brief Default destructor.
    virtual ~InterfaceChecker() = default;
    //! \brief Returns non-zero if the specified element has contributions.
    //! \param[in] iel Element number
    //! \param[in] I Index in first parameter direction of the element
    //! \param[in] J Index in second parameter direction of the element
    //! \param[in] K Index in third parameter direction of the element
    virtual short int hasContribution(int iel,
                                      int I, int J, int K = -1) const = 0;
    //! \brief Returns a status mask based on the element boundary parameters.
    //! \param[in] u0 Parameter value of the west element boundary
    //! \param[in] u1 Parameter value of the east element boundary
    //! \param[in] v0 Parameter value of the south element boundary
    //! \param[in] v1 Parameter value of the north element boundary
    //! \param[in] w0 Parameter value of the back element boundary
    //! \param[in] w1 Parameter value of the front element boundary
    virtual short int elmBorderMask(double u0, double u1,
                                    double v0, double v1,
                                    double w0 = 0.0, double w1 = 1.0) const
    {
      return 63; // = 1+2+4+8+16+32 (i.e., don't mask off any boundary)
    }
  };
}

#endif
