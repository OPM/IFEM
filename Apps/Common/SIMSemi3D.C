// $Id$
//==============================================================================
//!
//! \file SIMSemi3D.C
//!
//! \date Jun 5 2014
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Solution driver for plane-decoupled solvers.
//!
//==============================================================================

#include <memory>
#include <vector>
#include "DataExporter.h"

std::vector<DataExporter*> plane_exporters; //!< Data exporters for planes
std::vector<std::shared_ptr<std::ostream>> plane_log_files;
