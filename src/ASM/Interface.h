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

}

#endif
