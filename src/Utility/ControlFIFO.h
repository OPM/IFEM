// $Id$
//==============================================================================
//!
//! \file ControlFIFO.h
//!
//! \date Oct 7 2013
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Application control over a FIFO.
//!
//==============================================================================

#ifndef CONTROL_FIFO_H_
#define CONTROL_FIFO_H_

#include <string>
#include <map>

namespace tinyxml2 { class XMLElement; }


/*!
  \brief Callback for FIFO option handling.
*/

class ControlCallback
{
public:
  //! \brief Callback on receiving a XML control block.
  virtual void OnControl(const tinyxml2::XMLElement* context) = 0;
  //! \brief Returns context name for callback.
  virtual std::string GetContext() const = 0;
};


/*!
  \brief This class enables simple application control over a FIFO.

  \details A fifo is opened, and users can write instructions to the fifo
           in XML format. These are then processed between time steps.
*/

class ControlFIFO
{
public:
  //! \brief Default constructor.
  ControlFIFO() : fifo(-1) {}

  //! \brief The destructor tears down the opened fifo and removes the file.
  ~ControlFIFO();

  //! \brief Registers a callback handler.
  //! \param[in] callback The callback handler to register
  void registerCallback(ControlCallback& callback);

  //! \brief Opens the fifo and prepares for receiving.
  //! \param[in] name The name of the filesystem entry for the fifo
  bool open(const char* name = "ifem-control");

  //! \brief Polls for new data in the fifo.
  void poll();

private:
  std::string fifo_name; //!< Name of filesystem entry of our fifo
  int         fifo;      //!< fifo handle

  std::map<std::string,ControlCallback*> callbacks; //!< Registered callbacks
};

#endif
