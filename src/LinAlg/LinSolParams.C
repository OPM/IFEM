// $Id$
//==============================================================================
//!
//! \file LinSolParams.C
//!
//! \date Jan 29 2010
//!
//! \author Runar Holdahl / SINTEF
//!
//! \brief Linear solver parameters for PETSc matrices.
//!
//==============================================================================

#include "LinSolParams.h"
#include "Utilities.h"
#include "tinyxml2.h"
#include <cstdlib>


void SettingMap::addValue(const std::string& key, const std::string& value)
{
  values[key] = value;
}


std::string SettingMap::getStringValue(const std::string& key) const
{
  auto it = values.find(key);
  if (it != values.end())
    return it->second;

  return "";
}


int SettingMap::getIntValue(const std::string& key) const
{
  auto it = values.find(key);
  if (it != values.end())
    return atoi(it->second.c_str());

  return 0;
}


double SettingMap::getDoubleValue(const std::string& key) const
{
  auto it = values.find(key);
  if (it != values.end())
    return atof(it->second.c_str());

  return 0.0;
}


bool SettingMap::hasValue(const std::string& key) const
{
  return values.find(key) != values.end();
}


LinSolParams::BlockParams::BlockParams () : basis(1), comps(0)
{
  this->addValue("pc", "default");
  this->addValue("multigrid_ksp", "defrichardson");
}

LinSolParams::LinSolParams (LinAlg::LinearSystemType ls) : blocks(1), linSys(ls)
{
  this->addValue("type", "gmres");
  this->addValue("rtol", "1e-6");
  this->addValue("atol", "1e-20");
  this->addValue("dtol", "1e6");
  this->addValue("maxits", "1000");
  this->addValue("gmres_restart_iterations", "100");
  this->addValue("verbosity", "1");
}

LinSolParams::LinSolParams (const LinSolParams& p, LinAlg::LinearSystemType ls)
  : blocks(p.blocks), linSys(ls)
{
  values = p.values;
}


bool LinSolParams::BlockParams::read (const tinyxml2::XMLElement* elem,
                                      const std::string& prefix)
{
  utl::getAttribute(elem, "basis", basis);
  utl::getAttribute(elem, "components", comps);

  const char* value;
  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(child->Value(), "multigrid")) {
      std::string v;
      if (utl::getAttribute(child, "smoother", v))
        this->addValue("multigrid_smoother", v);
      if (utl::getAttribute(child, "levels", v))
        this->addValue("multigrid_levels", v);
      if (utl::getAttribute(child, "no_smooth", v)) {
        this->addValue("multigrid_no_smooth", v);
        this->addValue("multigrid_no_fine_smooth", v);
      }
      if (utl::getAttribute(child, "finesmoother", v))
        this->addValue("multigrid_finesmoother", v);
      if (utl::getAttribute(child, "multigrid_no_fine_smooth", v))
        this->addValue("multigrid_no_fine_smooth", v);
      if (utl::getAttribute(child, "ksp", v))
        this->addValue("multigrid_ksp", v);
      if (utl::getAttribute(child, "coarse_solver", v))
        this->addValue("multigrid_coarse_solver", v);
      if (utl::getAttribute(child, "max_coarse_size", v))
        this->addValue("multigrid_max_coarse_size", v);
    } else if (!strcasecmp(child->Value(),"dirsmoother")) {
      int order;
      std::string type;

      if (!utl::getAttribute(child,"type",type))
        return false;
      if (!utl::getAttribute(child,"order",order))
        return false;

      dirSmoother.push_back(DirSmoother{order, type});
    } else if (!strcasecmp(child->Value(), "asm")) {
      std::string v;
      if (utl::getAttribute(child, "nx", v))
        this->addValue("asm_nx", v);
      else if (utl::getAttribute(child, "ny", v))
        this->addValue("asm_ny", v);
      else if (utl::getAttribute(child, "nz", v))
        this->addValue("asm_nz", v);
      else if (utl::getAttribute(child, "overlap", v))
        this->addValue("asm_overlap", v);
    } else { // generic value or container tag - add with tag name as key
      std::string key = child->Value();
      if (child->FirstChildElement()) {
        if (!read(child, key+"_"))
          return false;
      } else
        if (value = utl::getValue(child, key.c_str()))
          this->addValue(prefix+key, value);
    }

  if (!this->hasValue("multigrid_finesmoother") && this->hasValue("multigrid_smoother"))
    this->addValue("multigrid_finesmoother", this->getStringValue("multigrid_smoother"));

  if (!this->hasValue("multigrid_finesmoother") && this->hasValue("multigrid_smoother"))
    this->addValue("multigrid_finesmoother", this->getStringValue("multigrid_smoother"));

  if (!this->hasValue("multigrid_no_smooth") || this->getIntValue("multigrid_no_smooth") < 1)
    this->addValue("multigrid_no_smooth", "1");

  return true;
}


bool LinSolParams::read (const tinyxml2::XMLElement* elem)
{
  if (elem->Attribute("verbosity"))
    this->addValue("verbosity", elem->Attribute("verbosity"));

  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  int parseblock = 0;
  for (; child; child = child->NextSiblingElement()) {
    const char* value;
    if ((value = utl::getValue(child,"type")))
      this->addValue("type", value);
    else if ((value = utl::getValue(child,"gmres_restart_iterations")))
      this->addValue("gmres_restart_iterations", value);
    else if ((value = utl::getValue(child,"reset_pc")))
      this->addValue("reset_pc", value);
    else if ((value = utl::getValue(child,"pc")))
      this->addValue("pc", value);
    else if ((value = utl::getValue(child,"schur")))
      this->addValue("schur", value);
    else if (!strcasecmp(child->Value(),"block")) {
      blocks.resize(++parseblock);
      blocks.back().read(child);
    }
    else if ((value = utl::getValue(child,"atol")))
      this->addValue("atol", value);
    else if ((value = utl::getValue(child,"rtol")))
      this->addValue("rtol", value);
    else if ((value = utl::getValue(child,"dtol")))
      this->addValue("dtol", value);
    else if ((value = utl::getValue(child,"maxits")))
      this->addValue("maxits", value);
  }

  if (parseblock == 0)
    blocks.back().read(elem);

  return true;
}
