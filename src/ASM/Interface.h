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

  struct Interface {
    int master; //!< Master patch (global number).
    int slave;  //!< Slave patch (global number).
    int midx;   //!< Index of boundary on master.
    int sidx;   //!< Index of boundary on slave.
    int orient; //!< Orientation.
    int dim;    //!< Dimension of boundary.
    int basis;  //!< Basis of boundary.
    int thick;  //!< Thickness of connection.
  };


  /*!
    \brief Base class for checking for internal boundary integrand contributions.
  */

  class InterfaceChecker
  {
  public:
    virtual ~InterfaceChecker () = default;
    //! \brief Returns non-zero if the specified element have contributions.
    //! \param[in] iel Element number
    //! \param[in] I Index in first parameter direction of the element
    //! \param[in] J Index in second parameter direction of the element
    //! \param[in] K Index in third parameter direction of the element
    virtual short int hasContribution(int iel, int I,
                                      int J, int K = -1) const = 0;
  };
}

#endif
