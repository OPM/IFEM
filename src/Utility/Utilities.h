// $Id$
//==============================================================================
//!
//! \file Utilities.h
//!
//! \date May 29 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Various utility methods.
//!
//==============================================================================

#ifndef _UTILITIES_H
#define _UTILITIES_H

#include "matrix.h"
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>

class TiXmlElement;
class TiXmlNode;

typedef std::map<int,int> IntMap; //!< Convenience type


namespace utl
{
  //! \brief Parses a character string into an integer or an integer range.
  //! \param values The integer value(s) is/are appended to this vector
  //! \param[in] argv Character string with integer data
  //!
  //! \details An integer range is recognised through the syntax \a i:j.
  void parseIntegers(std::vector<int>& values, const char* argv);

  //! \brief Parses a (possibly graded) sequence of knot values.
  //! \param xi The knot value(s) is/are appended to this vector
  //!
  //! \details The method uses the strtok function to parse a sequence of
  //! space-separated numbers and assumes that (at least one) call to strtok
  //! with the first argument different from nullptr has been made before
  //! invoking this method.
  bool parseKnots(std::vector<Real>& xi);

  //! \brief Reads one line, ignoring comment lines and leading blanks.
  //! \details The data read is kept in an internal static buffer.
  //! \param is File stream to read from
  //! \returns Pointer to the static buffer containg the data read.
  char* readLine(std::istream& is);
  //! \brief Ignores comment lines and blank lines from an input stream.
  //! \param is File stream to read from
  bool ignoreComments(std::istream& is);

  //! \brief Extracts a boolean attribute value from the specified XML-element.
  //! \param[in] xml Pointer to XML-element to extract from
  //! \param[in] att The attribute tag
  //! \param[out] val The attribute value
  //! \return \e true if the attribute \a att is found in \a xml,
  //! otherwise \e false
  bool getAttribute(const TiXmlElement* xml, const char* att, bool& val);
  //! \brief Extracts an integer attribute value from the specified XML-element.
  //! \param[in] xml Pointer to XML-element to extract from
  //! \param[in] att The attribute tag
  //! \param[out] val The attribute value
  //! \return \e true if the attribute \a att is found in \a xml,
  //! otherwise \e false
  bool getAttribute(const TiXmlElement* xml, const char* att, int& val);
  //! \brief Extracts a size_t attribute value from the specified XML-element.
  //! \param[in] xml Pointer to XML-element to extract from
  //! \param[in] att The attribute tag
  //! \param[out] val The attribute value
  //! \return \e true if the attribute \a att is found in \a xml,
  //! otherwise \e false
  bool getAttribute(const TiXmlElement* xml, const char* att, size_t& val);
  //! \brief Extracts a real attribute value from the specified XML-element.
  //! \param[in] xml Pointer to XML-element to extract from
  //! \param[in] att The attribute tag
  //! \param[out] val The attribute value
  //! \return \e true if the attribute \a att is found in \a xml,
  //! otherwise \e false
  bool getAttribute(const TiXmlElement* xml, const char* att, Real& val);
  //! \brief Extracts a string attribute value from the specified XML-element.
  //! \param[in] xml Pointer to XML-element to extract from
  //! \param[in] att The attribute tag
  //! \param[out] val The attribute value
  //! \param[in] toLower If \e true, convert return string to lower case
  //! \return \e true if the attribute \a att is found in \a xml,
  //! otherwise \e false
  bool getAttribute(const TiXmlElement* xml, const char* att, std::string& val,
                    bool toLower = false);
  //! \brief Returns the value (if any) of the specified XML-node.
  //! \param[in] xml Pointer to XML-node to extract the value from
  //! \param[in] tag The name of the XML-element to extract the value from
  const char* getValue(const TiXmlNode* xml, const char* tag);

  //! \brief Parses a sequence of knot values from the specified XML-node.
  //! \param[in] xml Pointer to XML-node to extract from
  //! \param xi The knot value(s) is/are appended to this vector
  bool parseKnots(const TiXmlNode* xml, std::vector<Real>& xi);

  //! \brief Transforms the integer value \a num into a unique range.
  //! \details This method is invoked on a series of (non-unique) values.
  //! When a value that has not been encountered before is detected,
  //! it is transformed into \a runner+1 and \a runner is incremented.
  //! \param num The integer value to transform, updated value on output
  //! \param runner The last new value assigned
  //! \param old2new Mapping from old values to new values
  //! return \e true if the value of \a num changed
  bool renumber(int& num, int& runner, IntMap& old2new);

  //! \brief Transforms the integer value \a num according to a given mapping.
  //! \details In this method the \a old2new mapping is not updated.
  //! \param num The integer value to transform, updated value on output
  //! \param[in] old2new Mapping from old values to new values
  //! \param[in] msg If \e true, give error message if given value not found
  //! \return \a true if \a num was found in the \a old2new mapping
  bool renumber(int& num, const IntMap& old2new, bool msg = false);

  //! \brief Compresses the rows of a 2D array based on given scatter indices.
  //! \param[in] index Scatter indices of the columns that should be retained
  //! \param[in] nr Number of rows in the 2D array
  //! \param[in] in The input array stored column-wise in a 1D array
  //! \param[out] out The output array stored column-wise in a 1D array
  //! \param[in] offset_in Optional start offset for the \a in vector
  int gather(const std::vector<int>& index, size_t nr,
             const std::vector<Real>& in, std::vector<Real>& out,
             size_t offset_in = 0);

  //! \brief Compresses the rows of a 2D array based on given scatter indices.
  //! \param[in] index Scatter indices of the columns that should be retained
  //! \param[in] nr Number of rows in the 2D array
  //! \param[in] in The input array stored column-wise in a 1D array
  //! \param[out] out The output array stored as a 2D matrix
  //! \param[in] offset_in Optional start offset for the \a in vector
  int gather(const std::vector<int>& index, size_t nr,
             const utl::vector<Real>& in, utl::matrix<Real>& out,
             size_t offset_in = 0);

  //! \brief Compresses a row of a 2D array based on given scatter indices.
  //! \param[in] index Scatter indices of the columns that should be retained
  //! \param[in] ir Index of the row to compress
  //! \param[in] nr Number of rows in the 2D array
  //! \param[in] in The input array stored column-wise in a 1D array
  //! \param[out] out The output array stored column-wise in a 1D array
  //! \param[in] offset_in Optional start offset for the \a in vector
  //! \param[in] shift_idx Optional constant shift in the scatter indices
  int gather(const std::vector<int>& index, size_t ir, size_t nr,
             const std::vector<Real>& in, std::vector<Real>& out,
             size_t offset_in = 0, int shift_idx = 0);

  //! \brief Searches for a real value in an ordered array of reals.
  //! \param[in] a The array of real values
  //! \param[in] v The value to search for
  //! \return 0-based index of the closest array value.
  //! \note It is assumed that the array \a a is sorted in encreasing order on
  //! input. If this is not the case the method will deliver incorrect result.
  size_t find_closest(const std::vector<Real>& a, Real v);

  //! \brief Prints a string stream to an output stream.
  //! \param out The output stream to print to
  //! \param[in] str The string stream to print
  //! \param[in] pid The PID of this process (only process 0 will print)
  //!
  //! \details This method makes sure to flush the output stream
  //! across all processes up front.
  void printSyncronized(std::ostream& out,
                        const std::stringstream& str, int pid);

  //! \brief Right-justifies the input string to the given total \a width.
  std::string adjustRight(size_t width, const std::string& s,
                          const std::string& suffix = " : ");

  //! \brief Splits an integer into its (unique) digits in ascending order.
  std::set<int> getDigits(int value);

  //! \brief Returns a const iterator to the entry with value \a iVal.
  IntMap::const_iterator findValue(const IntMap& iMap, int iVal);
  //! \brief Returns the key corresponding to the value \a iVal.
  //! \details If not in the map, the value \a iVal is returned.
  int findKey(const IntMap& iMap, int iVal);

  //! \brief Merges integer array \a a2 into array \a a1.
  //! \details Does not require the arrays to be sorted.
  //! The values of \a a2 not already in \a a1 are appended to \a a1.
  void merge(std::vector<int>& a1, const std::vector<int>& a2);
  //! \brief Merges real array \a a2 into array \a a1 based on array indices.
  //! \details Does not require the arrays to be sorted.
  //! The values of \a a2 not already in \a a1 are appended to \a a1.
  void merge(std::vector<Real>& a1, const std::vector<Real>& a2,
             const std::vector<int>& k1, const std::vector<int>& k2);
}

#endif
