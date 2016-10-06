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
  this->SIMdependency::registerDependency(sim,name,nvc);
  depFields.back().patches = patches;
  depFields.back().differentBasis = diffBasis;
  depFields.back().comp_use = component;
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
  ASMbase* patch;
  DepVector::const_iterator it;
  for (it = depFields.begin(); it != depFields.end(); ++it)
  {
    Vector* lvec = problem->getNamedVector(it->name);
    if (!lvec) continue; // Ignore fields without corresponding integrand vector

    const Vector* gvec = it->sim->getField(it->name);
    if (!gvec)
    {
      std::cerr <<" *** SIMdependency::extractPatchDependencies: \""
                << it->name <<"\" is not a registered field in the simulator \""
                << it->sim->getName() <<"\""<< std::endl;
      return false;
    }
    else if (gvec->empty()) {
      lvec->clear();
      continue; // No error, silently ignore empty fields (treated as zero)
    }

    patch = pindx < it->patches.size() ? it->patches[pindx] : model[pindx];
    // See ASMbase::extractNodeVec for interpretation of negative value on basis
    int basis = it->components < 0 ? it->components : it->differentBasis;
    patch->extractNodeVec(*gvec,*lvec,abs(it->components),basis);
    if (it->differentBasis > 0) {
      if (it->components == 1)
        problem->setNamedField(it->name,Field::create(patch,*lvec,
                                                      it->differentBasis,
                                                      it->comp_use));
      else
        problem->setNamedFields(it->name,Fields::create(patch,*lvec,
                                                        it->differentBasis));
    }
#if SP_DEBUG > 2
    std::cout <<"SIMdependency: Dependent field \""<< it->name
              <<"\" for patch "<< pindx+1 << *lvec;
#endif
  }

  return true;
}
