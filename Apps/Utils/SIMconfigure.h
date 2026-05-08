// $Id$
//==============================================================================
//!
//! \file SIMconfigure.h
//!
//! \date Oct 12 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief SIM solver configurator.
//!
//==============================================================================

#ifndef _SIM_CONFIGURE_H_
#define _SIM_CONFIGURE_H_


/*!
  \brief Struct for configuring a given simulator.
  \details Your SIM needs to specialize this for its type.
*/

template<class T> struct SolverConfigurator
{
  //! \brief Configures a simulator.
  //! \param sim The simulator to configure
  //! \param[in] props The setup properties for the simulator
  //! \param[in] infile The input file to parse
  int setup(T& sim, const typename T::SetupProps& props, char* infile);
};


//! \brief Configuration template.
//! \param sim The simulator to configure
//! \param[in] infile The input file to parse
//! \param[in] props The setup properties for the simulator
template<class T>
int ConfigureSIM(T& sim, char* infile,
                 const typename T::SetupProps& props = typename T::SetupProps())
{
  SolverConfigurator<T> setup;
  return setup.setup(sim,props,infile);
}

#endif
