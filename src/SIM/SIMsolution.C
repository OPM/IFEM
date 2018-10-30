// $Id$
//==============================================================================
//!
//! \file SIMsolution.C
//!
//! \date Nov 11 2017
//!
//! \author Knut Morten Okstad
//!
//! \brief General solution vector container for simulator drivers.
//!
//==============================================================================

#include "SIMsolution.h"
#ifdef HAS_CEREAL
#include <cereal/cereal.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#endif


bool SIMsolution::initSolution (size_t ndof, size_t nsol)
{
  // Note: Always at least one vector at this point, even if nsol is zero
  solution.resize(nsol > 1 ? nsol : 1);
  for (Vector& sol : solution)
    sol.resize(ndof,true);
  return true;
}


void SIMsolution::pushSolution (size_t nsol)
{
  if (solution.empty())
    return;
  else if (nsol == 0 || nsol > solution.size())
    nsol = solution.size();

  for (size_t n = nsol-1; n > 0; n--)
    std::copy(solution[n-1].begin(),solution[n-1].end(),solution[n].begin());
}


bool SIMsolution::saveSolution (SerializeMap& data,
                                const std::string& name) const
{
#ifdef HAS_CEREAL
  std::ostringstream str;
  {
    cereal::BinaryOutputArchive archive(str);
    for (const Vector& sol : this->getSolutions())
      archive(sol);
  }
  data.insert(std::make_pair(name,str.str()));
  return true;
#else
  return false;
#endif
}


bool SIMsolution::restoreSolution (const SerializeMap& data,
                                   const std::string& name)
{
#ifdef HAS_CEREAL
  SerializeMap::const_iterator sit = data.find(name);
  if (sit != data.end()) {
    std::stringstream str(sit->second);
    cereal::BinaryInputArchive archive(str);
    for (Vector& sol : this->theSolutions())
      archive(sol);
    return true;
  }
#endif
  return false;
}


std::string SIMsolution::serialize (const double* v, size_t n)
{
#ifdef HAS_CEREAL
  std::ostringstream str;
  cereal::BinaryOutputArchive archive(str);
  archive.saveBinary(v,n*sizeof(double));
  return str.str();
#else
  return "";
#endif
}


void SIMsolution::deSerialize (const std::string& data, double* v, size_t n)
{
#ifdef HAS_CEREAL
  std::stringstream str(data);
  cereal::BinaryInputArchive archive(str);
  archive.loadBinary(v,n*sizeof(double));
#endif
}
