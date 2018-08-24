// $Id$
//==============================================================================
//!
//! \file SIMdependency.C
//!
//! \date May 22 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Administration of simulators with dependencies to other simulators.
//!
//==============================================================================

#include "SIMdependency.h"
#include "IntegrandBase.h"
#include "ASMbase.h"
#include "Fields.h"
#include "Field.h"


void SIMdependency::registerDependency (SIMdependency* sim,
                                        const std::string& name, short int nvc,
                                        const PatchVec& patches,
                                        char diffBasis, int component)
{
  this->registerDependency(sim,name,nvc);
  depFields.back().patches = patches;
  depFields.back().differentBasis = diffBasis;
  depFields.back().comp_use = component;
}


void SIMdependency::registerDependency (SIMdependency* sim,
                                        const std::string& name, short int nvc,
                                        const PatchVec& patches,
                                        const int* MADOF)
{
  this->registerDependency(sim,name,nvc);
  depFields.back().patches = patches;
  depFields.back().MADOF = MADOF;
}


void SIMdependency::registerDependency (SIMdependency* sim,
                                        const std::string& name, short int nvc)
{
#ifdef SP_DEBUG
  std::cout <<"SIMdependency: Registering \""<< name
            <<"\" from "<< sim->getName()
            <<" as dependent field in "<< this->getName() << std::endl;
#endif
  depFields.push_back(Dependency(sim,name,nvc));
}


bool SIMdependency::fillField (const std::string& name,
                               const std::vector<double>& values)
{
  FieldMap::iterator it = myFields.find(name);
  if (it == myFields.end()) return false;

  *const_cast<utl::vector<double>*>(it->second) = values;
  return true;
}


utl::vector<double>*
SIMdependency::getField (const std::string& name)
{
  FieldMap::iterator it = myFields.find(name);
  if (it == myFields.end()) return nullptr;

  return const_cast<utl::vector<double>*>(it->second);
}


const utl::vector<double>*
SIMdependency::getField (const std::string& name) const
{
  FieldMap::const_iterator it = myFields.find(name);
  return it == myFields.end() ? nullptr : it->second;
}


SIMdependency::DepVector::const_iterator
SIMdependency::getDependency (const std::string& name) const
{
  std::vector<Dependency>::const_iterator it;
  for (it = depFields.begin(); it != depFields.end(); ++it)
    if (it->name == name) break;

  return it;
}


const utl::vector<double>*
SIMdependency::getDependentField (const std::string& name) const
{
  DepVector::const_iterator it = this->getDependency(name);
  if (it == depFields.end())
    return nullptr;

  return it->sim->getField(name);
}


ASMbase* SIMdependency::getDependentPatch (const std::string& name,
                                           int pindx) const
{
  DepVector::const_iterator it = this->getDependency(name);
  if (it == depFields.end())
    return nullptr;

  if (pindx < 0 || (size_t)pindx >= it->patches.size())
    return nullptr;

  return it->patches[pindx];
}


void SIMdependency::registerField (const std::string& name,
                                   const utl::vector<double>& vec)
{
  myFields[name] = &vec;
}


bool SIMdependency::extractPatchDependencies (IntegrandBase* problem,
                                              const PatchVec& model,
                                              size_t pindx) const
{
  for (const Dependency& dp : depFields)
  {
    Vector* lvec = problem->getNamedVector(dp.name);
    if (!lvec) continue; // Ignore fields without corresponding integrand vector

    const Vector* gvec = dp.sim->getField(dp.name);
    if (!gvec)
    {
      std::cerr <<" *** SIMdependency::extractPatchDependencies: \""
                << dp.name <<"\" is not a registered field in the simulator \""
                << dp.sim->getName() <<"\""<< std::endl;
      return false;
    }
    else if (gvec->empty())
    {
      lvec->clear();
      continue; // No error, silently ignore empty fields (treated as zero)
    }

    // See ASMbase::extractNodeVec for interpretation of negative value on basis
    int basis = dp.components < 0 ? dp.components : dp.differentBasis;
    ASMbase* pch = pindx < dp.patches.size() ? dp.patches[pindx] : model[pindx];
    if (dp.MADOF)
      pch->extractNodalVec(*gvec,*lvec,dp.MADOF);
    else if (dp.differentBasis && dp.components != pch->getNoFields(basis))
      pch->extractNodeVec(*gvec,*lvec);
    else
      pch->extractNodeVec(*gvec,*lvec,abs(dp.components),basis);

#if SP_DEBUG > 2
    std::cout <<"SIMdependency: Dependent field \""<< dp.name
              <<"\" for patch "<< pindx+1 << *lvec;
#endif
    if (dp.differentBasis > 0)
    {
      // Create a field object to handle different interpolation basis
      if (dp.components == 1)
        problem->setNamedField(dp.name,Field::create(pch,*lvec,
                                                     dp.differentBasis,
                                                     dp.comp_use));
      else
        problem->setNamedFields(dp.name,Fields::create(pch,*lvec,
                                                       dp.differentBasis));
    }
  }

  return true;
}
