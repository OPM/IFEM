// $Id$
//==============================================================================
//!
//! \file Utilities.C
//!
//! \date May 29 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Various utility methods.
//!
//==============================================================================

#include "Utilities.h"
#include "tinyxml.h"
#include <cstdlib>
#include <algorithm>

#ifdef HAVE_MPI
#include <mpi.h>
#endif


void utl::parseIntegers (std::vector<int>& values, const char* argv)
{
  if (!argv) return;

  char* endp = nullptr;
  int endVal = 0;
  values.push_back(strtol(argv,&endp,10));
  if (endp && *endp == ':')
    endVal = strtol(endp+1,&endp,10);
  while (values.back() < endVal)
    values.push_back(values.back()+1);
}


bool utl::parseKnots (std::vector<Real>& xi)
{
  char* cline = strtok(nullptr," ");
  if (!cline)
    return false;
  else if (isalpha(cline[0]))
  {
    // Geometric grading
    bool  biased = toupper(cline[0]) == 'B';
    int    ru    = atoi(strtok(nullptr," "));
    double alpha = atof(strtok(nullptr," "));
    double xi1   = (cline = strtok(nullptr," ")) ? atof(cline) : 0.0;
    double xi2   = (cline = strtok(nullptr," ")) ? atof(cline) : 1.0;
    if (xi1 < 0.0 || xi2 <= xi1 || xi2 > 1.0 || ru < 1)
      return false;
    if (biased && ru > 1 && alpha != 1.0)
      alpha = pow(alpha,1.0/Real(ru));

    double x = xi1;
    double D = xi2 - xi1;
    D *= (alpha == 1.0 ? 1.0/Real(ru+1) : (1.0-alpha)/(1.0-pow(alpha,ru+1)));
    if (xi1 > 0.0) xi.push_back(xi1);
    for (int i = 0; i < ru; i++)
    {
      x += D;
      D *= alpha;
      xi.push_back(x);
    }
    if (xi2 < 1.0) xi.push_back(xi2);
  }
  else
  {
    // Explicit specification of knots
    xi.push_back(atof(cline));
    while ((cline = strtok(nullptr," ")))
      xi.push_back(atof(cline));
  }

  return !xi.empty();
}


char* utl::readLine (std::istream& is)
{
  static char buf[1024];
  for (size_t n = 0; is.good(); n = strlen(buf))
  {
    for (size_t i = 0; i < n; i++)
      if (!isspace(buf[i]))
        if (buf[i] == '#')
          break; // comment line - skip
        else
          return buf+i; // first non-blank character

    is.getline(buf,sizeof(buf));
  }

  return 0;
}


bool utl::ignoreComments (std::istream& is)
{
  int c = ' ';
  while (isspace(c))
    c = is.get();

  while (c == '#')
  {
    is.ignore(256,'\n');
    c = is.get();
  }
  is.putback(c);

  return is.good();
}


bool utl::getAttribute (const TiXmlElement* xml, const char* att, bool& val)
{
  if (!xml || !xml->Attribute(att))
    return false;

  const char* value = xml->Attribute(att);
  if (!strcasecmp(value,"true") || !strcasecmp(value,"on"))
    val = true;
  else if (!strcasecmp(value,"false") || !strcasecmp(value,"off"))
    val = false;
  else if (value[0] == '1' || !strcasecmp(value,"yes"))
    val = true;
  else if (value[0] == '0' || !strcasecmp(value,"no"))
    val = false;
  else
    return false;

  return true;
}


bool utl::getAttribute (const TiXmlElement* xml, const char* att, int& val)
{
  if (xml && xml->Attribute(att))
    val = atoi(xml->Attribute(att));
  else
    return false;

  return true;
}


bool utl::getAttribute (const TiXmlElement* xml, const char* att, size_t& val)
{
  if (xml && xml->Attribute(att))
    val = atoi(xml->Attribute(att));
  else
    return false;

  return true;
}


bool utl::getAttribute (const TiXmlElement* xml, const char* att, Real& val)
{
  if (xml && xml->Attribute(att))
    val = atof(xml->Attribute(att));
  else
    return false;

  return true;
}


bool utl::getAttribute (const TiXmlElement* xml, const char* att,
                        std::string& val, bool toLower)
{
  if (!xml || !xml->Attribute(att))
    return false;

  val = xml->Attribute(att);

  if (toLower)
    for (size_t i = 0; i < val.size(); i++)
      val[i] = tolower(val[i]);

  return true;
}


/*!
  This method accepts two alternative ways of specifying the value \e myValue :
  \verbatim <name>myValue</name> \endverbatim and
  \verbatim <name value="myValue"/> \endverbatim
*/

const char* utl::getValue (const TiXmlNode* xml, const char* tag)
{
  if (xml->Value() && !strcasecmp(xml->Value(),tag))
  {
    const TiXmlElement* elm = dynamic_cast<const TiXmlElement*>(xml);
    if (elm && elm->Attribute("value"))
      return elm->Attribute("value");
    else if (xml->FirstChild())
      return xml->FirstChild()->Value();
  }

  return nullptr;
}


bool utl::parseKnots (const TiXmlNode* xml, std::vector<Real>& xi)
{
  if (!xml->FirstChild())
    return false;

  std::string xiVal("xi ");
  xiVal += xml->FirstChild()->Value();
  strtok(const_cast<char*>(xiVal.c_str())," ");
  return utl::parseKnots(xi);
}


bool utl::renumber (int& num, int& runner, IntMap& old2new)
{
  IntMap::iterator it = old2new.find(num);
  if (it == old2new.end())
    it = old2new.insert(std::make_pair(num,++runner)).first;

  if (num == it->second) return false;

  num = it->second;
  return true;
}


bool utl::renumber (int& num, const IntMap& old2new, bool msg)
{
  IntMap::const_iterator it = old2new.find(num);
  if (it == old2new.end())
  {
    if (msg)
      std::cerr <<" *** utl::renumber: Old value "<< num
                <<" does not exist in old2new mapping"<< std::endl;
    return false;
  }

  num = it->second;
  return true;
}


int utl::gather (const std::vector<int>& index, size_t nr,
                 const std::vector<Real>& in, std::vector<Real>& out,
                 size_t offset_in)
{
  int outside = 0;
  out.resize(nr*index.size());
  const Real* data = &in.front() + offset_in;
  Real* outVec = &out.front();
  for (size_t i = 0; i < index.size(); i++, outVec += nr)
    if (index[i] >= 0 && offset_in+nr*index[i] < in.size())
      memcpy(outVec,data+nr*index[i],nr*sizeof(Real));
    else if (index[i] >= 0)
      outside++;

  return outside;
}


int utl::gather (const std::vector<int>& index, size_t nr,
                 const utl::vector<Real>& in, utl::matrix<Real>& out,
                 size_t offset_in)
{
  int outside = 0;
  out.resize(nr,index.size());
  const Real* data = &in.front() + offset_in;
  for (size_t i = 0; i < index.size(); i++)
    if (index[i] >= 0 && offset_in+nr*index[i] < in.size())
      out.fillColumn(1+i,data+nr*index[i]);
    else if (index[i] >= 0)
      outside++;

  return outside;
}


int utl::gather (const std::vector<int>& index, size_t ir, size_t nr,
                 const std::vector<Real>& in, std::vector<Real>& out,
                 size_t offset_in, int shift_idx)
{
  if (ir >= nr) return index.size();

  int outside = 0;
  out.resize(index.size());
  offset_in += ir;
  for (size_t i = 0; i < index.size(); i++)
    if (index[i] >= shift_idx)
    {
      size_t ip = offset_in + nr*(index[i]-shift_idx);
      if (ip < in.size())
        out[i] = in[ip];
      else
        outside++;
    }
    else if (index[i] >= 0)
      outside++;

  return outside;
}


size_t utl::find_closest (const std::vector<Real>& a, Real v)
{
  // The lower_bound function uses binary search to find the index of the value.
  // Thus, this works only when the vector a is sorted in increasing order.
  size_t i = std::lower_bound(a.begin(),a.end(),v) - a.begin();
  if (i > 0 && i < a.size())
    return a[i]-v < v-a[i-1] ? i : i-1;
  else if (i == 0)
    return 0;
  else
    return i-1;
}


void utl::printSyncronized (std::ostream& out, const std::stringstream& str,
                            int pid)
{
  out << std::flush;
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if (pid == 0)
    out << str.str();
}


std::string utl::adjustRight (size_t width, const std::string& s,
                              const std::string& suffix)
{
  if (s.size() >= width) return s + suffix;

  static std::string blank(32,' ');
  return blank.substr(0,width-s.size()) + s + suffix;
}


std::set<int> utl::getDigits (int num)
{
  std::set<int> digits;

  for (; num != 0; num /= 10)
    digits.insert(num%10);

  return digits;
}


namespace utl
{
  /*!
    \brief A helper class used to search for values in an integer map.
  */
  class cmpInt
  {
    int myValue; //!< The integer value to search for
  public:
    //! The constructor initializes the value to search for.
    cmpInt(int value) : myValue(value) {}
    //! \brief Returns \e true if \a value.second equals \a myValue.
    bool operator()(const std::pair<int,int>& value) const
    {
      return value.second == myValue;
    }
  };
}


IntMap::const_iterator utl::findValue (const IntMap& iMap, int iVal)
{
  return std::find_if(iMap.begin(),iMap.end(),cmpInt(iVal));
}


int utl::findKey (const IntMap& iMap, int iVal)
{
  IntMap::const_iterator it = utl::findValue(iMap,iVal);
  return it == iMap.end() ? iVal : it->first;
}


void utl::merge (std::vector<int>& a1, const std::vector<int>& a2)
{
  for (size_t i = 0; i < a2.size(); i++)
    if (std::find(a1.begin(),a1.end(),a2[i]) == a1.end())
      a1.push_back(a2[i]);
}


void utl::merge (std::vector<Real>& a1, const std::vector<Real>& a2,
                 const std::vector<int>& k1, const std::vector<int>& k2)
{
  for (size_t i = 0; i < k2.size(); i++)
    if (std::find(k1.begin(),k1.end(),k2[i]) == k1.end())
      a1.push_back(a2[i]);
}


void utl::interleave(const std::vector<Real>& v1, const std::vector<Real>& v2,
                     std::vector<Real>& out, size_t n1, size_t n2)

{
  out.resize(v1.size()+v2.size());
  auto it_v1 = v1.begin();
  auto it_v2 = v2.begin();
  for (auto it_out = out.begin(); it_out != out.end();
            it_out += n1+n2, it_v1 += n1, it_v2 += n2) {
    std::copy(it_v1, it_v1+n1, it_out);
    std::copy(it_v2, it_v2+n2, it_out+n1);
  }
}
