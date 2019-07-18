// $Id$
//==============================================================================
//!
//! \file LogStream.C
//!
//! \date Mar 26 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Log stream class.
//!
//==============================================================================

#include "LogStream.h"


utl::LogStream::LogStream(std::ostream& out, int ppid, int mypid) :
  m_out(&out), m_ppid(ppid), m_pid(mypid)
{
}


utl::LogStream& utl::LogStream::operator<<(LogStream::StandardEndLine manip)
{
  if (m_pid == m_ppid && m_out)
    manip(*m_out);

  for (auto extra : m_extra)
    manip(*extra);

  return *this;
}


void utl::LogStream::addExtraLog(std::ostream* extra, bool clear)
{
  std::shared_ptr<std::ostream> file(extra);
  this->addExtraLog(file,clear);
}


void utl::LogStream::addExtraLog(std::shared_ptr<std::ostream> extra, bool clear)
{
  if (clear)
    m_extra.clear();

  m_extra.push_back(extra);
}


int utl::LogStream::precision(int streamsize)
{
  int result = streamsize;
  if (m_out) {
    result = m_out->precision();
    m_out->precision(streamsize);
  }
  for (auto it : m_extra)
    it->precision(streamsize);

  return result;
}


void utl::LogStream::flush()
{
  if (m_out)
    m_out->flush();
  for (auto it : m_extra)
    it->flush();
}


std::ios_base::fmtflags utl::LogStream::flags(std::ios_base::fmtflags flags)
{
  std::ios_base::fmtflags result = flags;
  if (m_out) {
    result = m_out->flags();
    m_out->flags(flags);
  }
  for (auto it : m_extra)
    it->flags(flags);

  return result;
}
