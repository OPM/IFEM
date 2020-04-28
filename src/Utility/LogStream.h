// $Id$
//==============================================================================
//!
//! \file LogStream.h
//!
//! \date Mar 26 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Log stream class.
//!
//==============================================================================

#ifndef LOG_STREAM_H_
#define LOG_STREAM_H_

#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>

namespace utl {

//! \brief Logging stream class.
class LogStream
{
public:
  //! \brief Default constructor.
  //! \param out The output stream to wrap
  explicit LogStream(std::ostream* out) : m_out(out) { m_ppid = m_pid = 0; }
  //! \brief Constructor initializing the output stream from a reference.
  //! \param out The output stream to wrap
  //! \param ppid The PID to print on
  //! \param mypid The PID of this process
  LogStream(std::ostream& out, int ppid = 0, int mypid = 0);

  //! \brief Nullifies the output stream.
  void setNull() { m_out = nullptr; }
  //! \brief Sets the output stream.
  void setStream(std::ostream& out) { m_out = &out; }
  //! \brief Sets PIDs for the stream.
  void setPIDs(int ppid, int mypid) { m_ppid = ppid; m_pid = mypid; }

  //! \brief Adds an extra logging stream.
  void addExtraLog(std::ostream* extra, bool clear = false);
  //! \brief Adds an extra logging stream.
  void addExtraLog(std::shared_ptr<std::ostream> extra, bool clear = false);
  //! \brief Drop an extra logging stream.
  void removeExtraLog(std::shared_ptr<std::ostream> extra);

  //! \brief Write data to stream
  template<typename T>
  LogStream& operator<<(const T& data) { return this->write(data); }

  //! \brief Write data to stream
  template<typename T>
  LogStream& write(const T& data)
  {
    if (m_ppid == m_pid && m_out)
      *m_out << data;
    for (auto extra : m_extra)
      *extra << data;

    return *this;
  }

  //! \brief Set precision of output
  int precision(int streamsize);
  //! \brief Get current precision
  int precision() const { return m_out ? m_out->precision() : 6; }

  //! \brief Check state of stream
  bool good() const { return m_out ? m_out->good() : true; }

  //! \brief Flush streams
  void flush();

  //! \brief Obtain stream flags
  std::ios_base::fmtflags flags() const { return m_out?m_out->flags():std::ios_base::fmtflags(); }
  //! \brief Set stream flags
  std::ios_base::fmtflags flags(std::ios_base::fmtflags fmtfl);

  //! \brief Assignment operator
  LogStream& operator=(const LogStream&);

  //! \brief This is the type of std::cout.
  typedef std::basic_ostream<char, std::char_traits<char> > CoutType;
  //! \brief Function pointer type.
  typedef CoutType& (*StandardEndLine)(CoutType&);

  //! \brief Handling of the std::endl manipulator with friends
  LogStream& operator << (StandardEndLine manip);

protected:
  std::ostream* m_out; //!< Main output stream
  std::vector<std::shared_ptr<std::ostream>> m_extra; //!< Extra output streams
  int m_ppid; //!< PID to print on
  int m_pid;  //!< This process' PID
};

}

#endif
