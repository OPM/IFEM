// $Id$
//==============================================================================
//!
//! \file TextureProperties.h
//!
//! \date Jan 3 2020
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Properties defined through a texture map.
//!
//==============================================================================

#ifndef TEXTURE_PROPERTIES_H_
#define TEXTURE_PROPERTIES_H_

#include "Function.h"
#include "MatVec.h"

#include <map>
#include <memory>

class TiXmlElement;
class Vec3;


//! \brief Class containing a set of properties defined through a texture map.
class TextureProperties {
public:
  //! \brief Parse an XML definition.
  //! param elem XML element to parse
  void parse (const TiXmlElement* elem);

  //! \brief Print property information to log.
  void printLog() const;

  //! \brief Get value for a property
  //! \param[in] name Name of property
  //! \param[in] X Position (including parameter values) to evaluate property for
  //! \param[out] val Property value
  bool getProperty(const std::string& name, const Vec3& X, double& val) const;

  //! \brief Get value for a vector property
  //! \param[in] name Name of property
  //! \param[in] X Position (including parameter values) to evaluate property for
  //! \param[out] val Property value
  bool getProperty(const std::string& name, const Vec3& X, Vec3& val) const;

  //! \brief Check if a property is available.
  //! \param name Name of property
  bool hasProperty(const std::string& name) const;

protected:
  //! \brief Struct holding information about a property.
  struct Property {
    double min = 0.0; //!< Minimum value
    double max = 1.0; //!< Maximum value
    Matrix3D textureData; //!< Texture data
    bool prescaled = false; //!< True if data is already scaled
    std::string func_definition; //!< Non-empty if we have a function of the texture value
    std::unique_ptr<FunctionBase> function; //!< Function definition for property (if any)
  };

  //! \brief Obtains texture value for a property.
  double getValue(const Property& prop, const Vec3& X) const;

  std::map<std::string, Property> properties; //!< Map of available properties
};


//! \brief Class to use a property as a function.
template<class Base, class Value>
class PropertyFuncType : public Base {
public:
  //! \brief Constructor initializes the members.
  //! \param prop Name of property
  //! \param props Texture property container
  PropertyFuncType(const std::string& prop, const TextureProperties& props) :
    m_prop(prop), m_props(props) {}

  //! \brief Empty destructor.
  virtual ~PropertyFuncType() {}

  //! \brief Evaluate function in a point.
  //! \param X Position to evaluate in
  Value evaluate(const Vec3& X) const override
  {
    Value val;
    m_props.getProperty(m_prop, X, val);
    return val;
  }

protected:
  std::string m_prop; //!< Name of property
  const TextureProperties& m_props; //!< Texture properties container
};


using PropertyFunc = PropertyFuncType<RealFunc,Real>; //!< Convenience type alias for scalars
using PropertyVecFunc = PropertyFuncType<VecFunc,Vec3>; //!< Convenience type alias for vector

#endif
