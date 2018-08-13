// $Id$
//==============================================================================
//!
//! \file ControlFIFO.C
//!
//! \date Oct 7 2013
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Application control over a FIFO.
//!
//==============================================================================

#include "ControlFIFO.h"
#include "tinyxml.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>


ControlFIFO::~ControlFIFO ()
{
  if (fifo == -1)
    return;

  close(fifo);
  unlink(fifo_name.c_str());
}


bool ControlFIFO::open (const char* name)
{
  fifo_name = name;
  unlink(name);
#if !defined(__MINGW64__) || !defined(__MINGW32__)
  if (mkfifo(name, S_IWUSR | S_IRUSR | S_IRGRP | S_IROTH) == 0) {
    fifo = ::open(name, O_RDONLY | O_NONBLOCK);
    return true;
  }
  std::cerr <<" *** Error creating control fifo '"<< fifo_name <<"'\n"
            <<"     Application will not be externally controllable."
            << std::endl;
#endif
  return false;
}


void ControlFIFO::poll ()
{
  if (fifo == -1)
    return;

  char temp[2048];
  temp[0] = '\0';
  int len = read(fifo, temp, 2048);
  if (len < 0)
    return;

  temp[len] = '\0';
  if (!strlen(temp))
    return;

  TiXmlDocument doc;
  doc.Parse(temp);
  if (!doc.RootElement()) {
    std::cerr <<" *** Invalid control data received:\n"<< temp << std::endl;
    return;
  }

  TiXmlElement* elem = doc.RootElement()->FirstChildElement();
  for (; elem; elem = elem->NextSiblingElement())
    if (callbacks.find(elem->Value()) != callbacks.end())
      callbacks[elem->Value()]->OnControl(elem);
}


void ControlFIFO::registerCallback (ControlCallback& callback)
{
  callbacks.insert(std::make_pair(callback.GetContext(),&callback));
}
