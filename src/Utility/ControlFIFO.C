#include "ControlFIFO.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>


ControlFIFO::ControlFIFO() : fifo(-1)
{
}


ControlFIFO::~ControlFIFO()
{
  if (fifo > -1) {
    close(fifo);
    unlink(fifo_name.c_str());
  }
}


bool ControlFIFO::open(const std::string& name)
{
#if !defined(__MINGW64__) || !defined(__MINGW32__)
  fifo_name = name;
  unlink(fifo_name.c_str());
  if (mkfifo(fifo_name.c_str(), S_IWUSR | S_IRUSR | S_IRGRP | S_IROTH)) {

    std::cerr << "Error creating control fifo '" << fifo_name << "'" << std::endl
              << "=> Application will not be externally controllable" << std::endl;
    return false;
  }

  fifo = ::open("ifem-control", O_RDONLY | O_NONBLOCK);

  return true;
#endif
  return false;
}


void ControlFIFO::poll()
{
  if (fifo == -1)
    return;

  char temp[2048];
  temp[0] = '\0';

  int len = read(fifo, temp, 2048);
  if (len > -1)
    temp[len] = 0;
  else
    return;

  TiXmlDocument doc;

  if(!strlen(temp))
    return;

  if (strlen(temp) && !doc.Parse(temp) || !doc.RootElement()) {
    std::cerr << "Invalid control data received: " << std::endl
              << temp << std::endl;
    return;
  }
  TiXmlElement* elem = doc.RootElement()->FirstChildElement();
  while(elem) {
    CallbackMap::iterator it;
    if ((it=callbacks.find(elem->Value())) != callbacks.end())
      callbacks[elem->Value()]->OnControl(elem);
    elem = elem->NextSiblingElement();
  }
}


void ControlFIFO::registerCallback(ControlCallback& callback)
{
  callbacks.insert(std::make_pair(callback.GetContext(), &callback));
}
