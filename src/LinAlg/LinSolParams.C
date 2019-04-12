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
#include "tinyxml.h"
#include <fstream>
#include <sstream>
#include <utility>
#include <iterator>


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


LinSolParams::BlockParams::BlockParams() :
  basis(1), comps(0)
{
  addValue("pc", "default");
  addValue("multigrid_ksp", "defrichardson");
}

LinSolParams::LinSolParams()
  : blocks(1)
{
  addValue("type", "gmres");
  addValue("rtol", "1e-6");
  addValue("atol", "1e-20");
  addValue("dtol", "1e6");
  addValue("maxits", "1000");
  addValue("gmres_restart_iterations", "100");
  addValue("verbosity", "1");
}


bool LinSolParams::BlockParams::read(const TiXmlElement* elem, const std::string& prefix)
{
  utl::getAttribute(elem, "basis", basis);
  utl::getAttribute(elem, "components", comps);

  const char* value;
  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(child->Value(), "multigrid")) {
      std::string v;
      if (utl::getAttribute(child, "smoother", v))
        addValue("multigrid_smoother", v);
      if (utl::getAttribute(child, "levels", v))
        addValue("multigrid_levels", v);
      if (utl::getAttribute(child, "no_smooth", v)) {
        addValue("multigrid_no_smooth", v);
        addValue("multigrid_no_fine_smooth", v);
      }
      if (utl::getAttribute(child, "finesmoother", v))
        addValue("multigrid_finesmoother", v);
      if (utl::getAttribute(child, "multigrid_no_fine_smooth", v))
        addValue("multigrid_no_fine_smooth", v);
      if (utl::getAttribute(child, "ksp", v))
        addValue("multigrid_ksp", v);
      if (utl::getAttribute(child, "coarse_solver", v))
        addValue("multigrid_coarse_solver", v);
      if (utl::getAttribute(child, "max_coarse_size", v))
        addValue("multigrid_max_coarse_size", v);
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
        addValue("asm_nx", v);
      else if (utl::getAttribute(child, "ny", v))
        addValue("asm_ny", v);
      else if (utl::getAttribute(child, "nz", v))
        addValue("asm_nz", v);
      else if (utl::getAttribute(child, "overlap", v))
        addValue("asm_overlap", v);
    } else { // generic value or container tag - add with tag name as key
      std::string key = child->Value();
      if (child->FirstChildElement()) {
        if (!read(child, key+"_"))
          return false;
      } else
        if (value = utl::getValue(child, key.c_str()))
          addValue(prefix+key, value);
    }

  if (!hasValue("multigrid_finesmoother") && hasValue("multigrid_smoother"))
    addValue("multigrid_finesmoother", getStringValue("multigrid_smoother"));

  if (!hasValue("multigrid_finesmoother") && hasValue("multigrid_smoother"))
    addValue("multigrid_finesmoother", getStringValue("multigrid_smoother"));

  if (!hasValue("multigrid_no_smooth") || getIntValue("multgrid_no_smooth") < 1)
    addValue("multigrid_no_smooth", "1");

  return true;
}


bool LinSolParams::read (const TiXmlElement* elem)
{
  if (elem->Attribute("verbosity"))
    addValue("verbosity", elem->Attribute("verbosity"));

  const TiXmlElement* child = elem->FirstChildElement();
  int parseblock = 0;
  for (; child; child = child->NextSiblingElement()) {
    const char* value;
    if ((value = utl::getValue(child,"type")))
      addValue("type", value);
    else if ((value = utl::getValue(child,"gmres_restart_iterations")))
      addValue("gmres_restart_iterations", value);
    else if ((value = utl::getValue(child,"pc")))
      addValue("pc", value);
    else if ((value = utl::getValue(child,"schur")))
      addValue("schur", value);
    else if (!strcasecmp(child->Value(),"block")) {
      blocks.resize(++parseblock);
      blocks.back().read(child);
    }
    else if ((value = utl::getValue(child,"atol")))
      addValue("atol", value);
    else if ((value = utl::getValue(child,"rtol")))
      addValue("rtol", value);
    else if ((value = utl::getValue(child,"dtol")))
      addValue("dtol", value);
    else if ((value = utl::getValue(child,"maxits")))
      addValue("maxits", value);
  }

  if (parseblock == 0)
    blocks.back().read(elem);

  return true;
}
