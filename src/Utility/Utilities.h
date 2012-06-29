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
#include <iostream>
#include <vector>
#include <map>

class TiXmlElement;


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
  //! with the first argument different from NULL has been made before
  //! invoking this method.
  bool parseKnots(std::vector<real>& xi);

  //! \brief Reads one line, ignoring comment lines and leading blanks.
  //! \details The data read is kept in an internal static buffer.
  //! \param is File stream to read from
  //! \returns Pointer to the static buffer containg the data read.
  char* readLine(std::istream& is);
  //! \brief Ignores comment lines and blank lines from an input stream.
  //! \param is File stream to read from
  bool ignoreComments(std::istream& is);

  //! \brief Extracts an integer attribute value from the specified XML-element.
  //! \param[in] xml Pointer to XML-element to extract from
  //! \param[in] att The attribute tag
  //! \param[out] val The attribute value
  //! \return \e true if the attribute \a att is found in \a child,
  //! otherwise \e false
  bool getAttribute(const TiXmlElement* xml, const char* att, int& val);
  //! \brief Extracts an unsigned integer attribute value from the specified XML-element.
  //! \param[in] xml Pointer to XML-element to extract from
  //! \param[in] att The attribute tag
  //! \param[out] val The attribute value
  //! \return \e true if the attribute \a att is found in \a child,
  //! otherwise \e false
  bool getAttribute(const TiXmlElement* xml, const char* att, size_t& val);
  //! \brief Extracts a real attribute value from the specified XML-element.
  //! \param[in] xml Pointer to XML-element to extract from
  //! \param[in] att The attribute tag
  //! \param[out] val The attribute value
  //! \return \e true if the attribute \a att is found in \a child,
  //! otherwise \e false
  bool getAttribute(const TiXmlElement* xml, const char* att, real& val);
  //! \brief Extracts a string attribute value from the specified XML-element.
  //! \param[in] xml Pointer to XML-element to extract from
  //! \param[in] att The attribute tag
  //! \param[out] val The attribute value
  //! \param[in] toLower If \e true, convert return string to lower case
  //! \return \e true if the attribute \a att is found in \a child,
  //! otherwise \e false
  bool getAttribute(const TiXmlElement* xml, const char* att, std::string& val,
                    bool toLower = false);
  //! \brief Returns the value (if any) of the specified XML-element.
  //! \param[in] xml Pointer to XML-element to extract the value from
  //! \param[in] tag The name of the XML-element to extract the value from
  const char* getValue(const TiXmlElement* xml, const char* tag);

  //! \brief Parses a sequence of knot values from the specified XML-element.
  //! \param[in] xml Pointer to XML-element to extract from
  //! \param xi The knot value(s) is/are appended to this vector
  bool parseKnots(const TiXmlElement* xml, std::vector<real>& xi);

  //! \brief Transforms the integer value \a num into a unique range.
  //! \details This method is invoked on a series of (non-unique) values.
  //! When a value that has not been encountered before is detected,
  //! it is transformed into \a runner+1 and \a runner is incremented.
  //! \param num The integer value to transform, updated value on output
  //! \param runner The last new value assigned
  //! \param old2new Mapping from old values to new values
  //! return \e true if the value of \a num changed
  bool renumber(int& num, int& runner, std::map<int,int>& old2new);

  //! \brief Transforms the integer value \a num according to a given mapping.
  //! \details In this method the \a old2new mapping is not updated.
  //! \param num The integer value to transform, updated value on output
  //! \param[in] old2new Mapping from old values to new values
  //! \param[in] msg If \e true, give error message if given value not found
  //! \return \a true if \a num was found in the \a old2new mapping
  bool renumber(int& num, const std::map<int,int>& old2new, bool msg = false);

  //! \brief Compresses the rows of a 2D array based on given scatter indices.
  //! \param[in] index Scatter indices of the columns that should be retained
  //! \param[in] nr Number of rows in the 2D array
  //! \param[in] in The input array stored column-wise in a 1D array
  //! \param[out] out The output array stored column-wise in a 1D array
  //! \param[in] offset_in Optional start offset for the \a in vector
  int gather(const std::vector<int>& index, size_t nr,
	     const std::vector<real>& in, std::vector<real>& out,
	     size_t offset_in = 0);

  //! \brief Compresses the rows of a 2D array based on given scatter indices.
  //! \param[in] index Scatter indices of the columns that should be retained
  //! \param[in] nr Number of rows in the 2D array
  //! \param[in] in The input array stored column-wise in a 1D array
  //! \param[out] out The output array stored as a 2D matrix
  //! \param[in] offset_in Optional start offset for the \a in vector
  int gather(const std::vector<int>& index, size_t nr,
	     const utl::vector<real>& in, utl::matrix<real>& out,
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
             const std::vector<real>& in, std::vector<real>& out,
             size_t offset_in = 0, int shift_idx = 0);

  //! \brief Returns the number of monomials in Pascal's triangle.
  //! \param[in] p Polynomial order (>= 0)
  //! \param[in] nsd Number of spatial dimensions (2 or 3)
  size_t Pascal(int p, unsigned short int nsd);
  //! \brief Evaluates the monomials of Pascal's triangle in 2D for order \a p.
  void Pascal(int p, real x, real y, std::vector<real>& phi);
  //! \brief Evaluates the monomials of Pascal's triangle in 3D for order \a p.
  void Pascal(int p, real x, real y, real z, std::vector<real>& phi);

  //! \brief Searches for a real value in an ordered array of reals.
  //! \param[in] a The array of real values
  //! \param[in] v The value to search for
  //! \return 0-based index of the closest array value.
  //! \note It is assumed that the array \a a is sorted in encreasing order on
  //! input. If this is not the case the method will deliver incorrect result.
  size_t find_closest(const std::vector<real>& a, real v);
}

#endif
