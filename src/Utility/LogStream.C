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

namespace utl {

nullstream NullStream;


LogStream::LogStream(std::ostream& out, int ppid, int mypid) :
  m_out(&out), m_ppid(ppid), m_pid(mypid)
{
}


LogStream& LogStream::operator<<(LogStream::StandardEndLine manip)
{
  if (m_pid == m_ppid && m_out)
    manip(*m_out);

  for (auto extra : m_extra)
    manip(*extra);

  return *this;
}


LogStream& LogStream::operator=(const LogStream& log2)
{
  m_out = log2.m_out;
  m_extra = log2.m_extra;
  m_ppid = log2.m_ppid;
  m_pid = log2.m_pid;

  return *this;
}


void LogStream::addExtraLog(std::shared_ptr<std::ostream>& extra, bool clear)
{
  if (clear)
    m_extra.clear();

  m_extra.push_back(extra);
}


int LogStream::precision() const
{
  return m_out?m_out->precision():6;
}


int LogStream::precision(int streamsize)
{
  int result = streamsize;
  if (m_out)
    result = m_out->precision(streamsize);
  for (auto it : m_extra)
    it->precision(streamsize);

  return result;
}


void LogStream::flush()
{
  if (m_out)
    m_out->flush();
  for (auto it : m_extra)
    it->flush();
}


std::ios_base::fmtflags LogStream::flags(std::ios_base::fmtflags flags)
{
  std::ios_base::fmtflags result = flags;
  if (m_out)
    result = m_out->flags(flags);
  for (auto it : m_extra)
    it->flags(flags);

  return result;
}

}
