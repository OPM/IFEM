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

#ifdef USE_OPENMP
#include <omp.h>
#endif


void utl::parseIntegers (std::vector<int>& values, const char* argv)
{
  if (!argv) return;

  char* endp = 0;
  int endVal = 0;
  values.push_back(strtol(argv,&endp,10));
  if (endp && *endp == ':')
    endVal = strtol(endp+1,&endp,10);
  while (values.back() < endVal)
    values.push_back(values.back()+1);
}


bool utl::parseKnots (std::vector<real>& xi)
{
  char* cline = strtok(NULL," ");
  if (toupper(cline[0]) == 'G')
  {
    // Geometric grading
    int    ru    = atoi(strtok(NULL," "));
    double alpha = atof(strtok(NULL," "));
    double xi1   = (cline = strtok(NULL," ")) ? atof(cline) : 0.0;
    double xi2   = (cline = strtok(NULL," ")) ? atof(cline) : 1.0;
    if (xi1 < 0.0 || xi2 <= xi1 || xi2 > 1.0 || ru < 1)
      return false;

    double D1 = 0.0;
    double D2 = (xi2-xi1);
    D2 *= (alpha <= 1.0 ? 1.0/real(ru+1) : (1.0-alpha)/(1.0-pow(alpha,ru+1)));
    if (xi1 > 0.0) xi.push_back(xi1);
    for (int i = 0; i < ru; i++)
    {
      xi.push_back(xi1+D1+D2);
      D1 = D2;
      if (alpha > 1.0) D2 = alpha*D1;
    }
    if (xi2 < 1.0) xi.push_back(xi2);
  }
  else
  {
    // Explicit specification of knots
    xi.push_back(atof(cline));
    while ((cline = strtok(NULL," ")))
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


bool utl::getAttribute (const TiXmlElement* xml, const char* att, int& val)
{
  if (xml->Attribute(att))
    val = atoi(xml->Attribute(att));
  else
    return false;

  return true;
}


bool utl::getAttribute (const TiXmlElement* xml, const char* att, real& val)
{
  if (xml->Attribute(att))
    val = atof(xml->Attribute(att));
  else
    return false;

  return true;
}


bool utl::getAttribute (const TiXmlElement* xml, const char* att,
			std::string& val, bool toLower)
{
  if (!xml->Attribute(att))
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

const char* utl::getValue (const TiXmlElement* xml, const char* tag)
{
  if (xml->Value() && !strcasecmp(xml->Value(),tag))
  {
    if (xml->Attribute("value"))
      return xml->Attribute("value");
    else if (xml->FirstChild())
      return xml->FirstChild()->Value();
  }

  return NULL;
}


bool utl::parseKnots (const TiXmlElement* xml, std::vector<real>& xi)
{
  std::string xiVal("xi ");
  xiVal += xml->FirstChildElement()->Value();
  strtok(const_cast<char*>(xiVal.c_str())," ");
  return utl::parseKnots(xi);
}


bool utl::renumber (int& num, int& runner, std::map<int,int>& old2new)
{
  std::map<int,int>::iterator it = old2new.find(num);
  if (it == old2new.end())
    it = old2new.insert(std::make_pair(num,++runner)).first;

  if (num == it->second) return false;

  num = it->second;
  return true;
}


bool utl::renumber (int& num, const std::map<int,int>& old2new, bool msg)
{
  std::map<int,int>::const_iterator it = old2new.find(num);
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
		 const std::vector<real>& in, std::vector<real>& out,
		 size_t offset_in)
{
  int outside = 0;
  out.resize(nr*index.size());
  const real* data = &in.front() + offset_in;
  real* outVec = &out.front();
  for (size_t i = 0; i < index.size(); i++, outVec += nr)
    if (index[i] >= 0 && offset_in+nr*index[i] < in.size())
      memcpy(outVec,data+nr*index[i],nr*sizeof(real));
    else
      outside++;

  return outside;
}


int utl::gather (const std::vector<int>& index, size_t nr,
		 const utl::vector<real>& in, utl::matrix<real>& out,
		 size_t offset_in)
{
  int outside = 0;
  out.resize(nr,index.size());
  const real* data = &in.front() + offset_in;
  for (size_t i = 0; i < index.size(); i++)
    if (index[i] >= 0 && offset_in+nr*index[i] < in.size())
      out.fillColumn(1+i,data+nr*index[i]);
    else
      outside++;

  return outside;
}


int utl::gather (const std::vector<int>& index, size_t ir, size_t nr,
                 const std::vector<real>& in, std::vector<real>& out,
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
    else
      outside++;

  return outside;
}


size_t utl::Pascal (int p, unsigned short int nsd)
{
  size_t nM = 1;
  for (int q = 1; q <= p; q++)
    for (int i = q; i >= 0; i--)
      if (nsd == 2)
	nM++;
      else for (int j = q; j >= 0; j--)
        if (i+j <= q) nM++;

  return nM;
}


void utl::Pascal (int p, real x, real y, std::vector<real>& phi)
{
  phi.clear();
  phi.reserve(Pascal(p,2));
  phi.push_back(real(1));
  for (int q = 1; q <= p; q++)
    for (int i = q; i >= 0; i--)
    {
      int k, j = q-i;
      real a = real(1);
      for (k = 0; k < i; k++) a *= x;
      for (k = 0; k < j; k++) a *= y;
      phi.push_back(a);
    }
}


void utl::Pascal (int p, real x, real y, real z, std::vector<real>& phi)
{
  phi.clear();
  phi.reserve(Pascal(p,3));
  phi.push_back(real(1));
  for (int q = 1; q <= p; q++)
    for (int i = q; i >= 0; i--)
      for (int j = q; j >= 0; j--)
	if (i+j <= q)
	{
	  int l, k = q-i-j;
	  real a = real(1);
	  for (l = 0; l < i; l++) a *= x;
	  for (l = 0; l < j; l++) a *= y;
	  for (l = 0; l < k; l++) a *= z;
	  phi.push_back(a);
	}
}


size_t utl::find_closest (const std::vector<real>& a, real v)
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


void utl::calcThreadGroups(int nel1, int nel2, ThreadGroups& result)
{
  int threads=1;
  int groups=1;
  int stripsize=1;
  int remainder=0;
  int dir=0, mul=1;
#ifdef USE_OPENMP
  threads = omp_get_max_threads();
  if (threads > 1)
    groups = 2;

  int els;
  int s1 = nel1/(groups*threads);
  int s2 = nel2/(groups*threads);
  int r1 = nel1-(s1*groups*threads);
  int r2 = nel2-(s2*groups*threads);
  if (r1*nel2 < r2*nel1) {
    stripsize = s1;
    dir = 0;
    els = nel1;
    mul = 1;
  } else {
    stripsize = s2;
    els = nel2;
    dir = 1;
    mul = nel1;
  }

  if (stripsize < 2 && groups > 1) {
    std::cerr << __FUNCTION__ << ": Warning: too many threads available." << std::endl
              << "Reducing to a suitable amount" << std::endl;
    while (((stripsize = els/(groups*threads)) < 2) && threads > 1)
      threads--;
    if (threads == 1)
      groups=1;
    stripsize = els/(groups*threads);
  }
  remainder = els-(stripsize*groups*threads);

#if SP_DEBUG > 1
  std::cout << "we have " << threads << " threads available" << std::endl;
  std::cout << "nel1 " << nel1 << std::endl;
  std::cout << "nel2 " << nel2 << std::endl;
  std::cout << "stripsize " << stripsize << std::endl;
  std::cout << "# of strips " << els/stripsize << std::endl;
  std::cout << "remainder " << remainder << std::endl;
#endif
#endif
  result.resize(groups);

  if (groups == 1) {
    result[0].resize(1);
    for (int i=0;i<nel1*nel2;++i)
      result[0][0].push_back(i);
  } else {
    std::vector< std::vector<int> > stripsizes;
    stripsizes.resize(2);
    stripsizes[0].resize(threads,stripsize);
    stripsizes[1].resize(threads,stripsize);
    int r=0;
    for (int i=0;i<remainder && r < remainder;++i) {
      stripsizes[1][threads-1-i]++;
      r++;
      if (r < remainder) {
        stripsizes[0][threads-1-i]++;
        r++;
      }
    }
    std::vector< std::vector<int> > startelms;
    startelms.resize(2);
    int offs=0;
    for (int i=0;i<threads;++i) {
      startelms[0].push_back(offs*mul);
      offs += stripsizes[0][i];
      startelms[1].push_back(offs*mul);
      offs += stripsizes[1][i];
    }
    for (size_t g=0;g<result.size();++g) { // loop over groups
      result[g].resize(threads);
      for (int t=0;t<threads;++t) { // loop over threads
        int maxx = dir==0?stripsizes[g][t]:nel1;
        int maxy = dir==1?stripsizes[g][t]:nel2;
        for (int i2=0; i2 < maxy; ++i2) { // loop in y direction
          for (int i1=0;i1<maxx; ++i1) {
            int iEl = startelms[g][t]+i1+i2*nel1;
            result[g][t].push_back(iEl);
          }
        }
      }
    }
  }

#if defined(USE_OPENMP) && SP_DEBUG > 1
  for (size_t i=0;i<result.size();++i) {
    std::cout << "group " << i << std::endl;
    for (size_t j=0;j<result[i].size();++j) {
      std::cout << "\t thread " << j << ": ";
      for (size_t k=0;k<result[i][j].size();++k)
        std::cout << result[i][j][k] << " ";
      std::cout << std::endl;
    }
  }
#endif
}


void utl::calcThreadGroups(int nel1, int nel2, int nel3, ThreadGroups& result)
{
  int threads=1;
  int groups=1;
  int stripsize=1;
  int remainder=0;
  int dir=0, mul=1;
#ifdef USE_OPENMP
  threads = omp_get_max_threads();
  if (threads > 1)
    groups = 2;

  int els;
  int s1 = nel1/(groups*threads);
  int s2 = nel2/(groups*threads);
  int s3 = nel3/(groups*threads);
  int r1 = nel1-(s1*groups*threads);
  int r2 = nel2-(s2*groups*threads);
  int r3 = nel3-(s3*groups*threads);
  if (r1*nel2*nel3 < r2*nel1*nel3 && r1*nel2*nel3 < r3*nel1*nel2 ) {
    // strips along x axis
    stripsize = s1;
    dir = 0;
    els = nel1;
    mul = 1;
  } else if (r2*nel1*nel3 < r1*nel2*nel3 && r2*nel1*nel3 < r3*nel1*nel2 ) {
    // strips along y axis
    stripsize = s2;
    els = nel2;
    dir = 1;
    mul = nel1;
  } else {
    // strips along z axis
    stripsize = s3;
    els = nel3;
    dir = 2;
    mul = nel1*nel2;
  }

  if (stripsize < 2 && groups > 1) {
    std::cerr << __FUNCTION__ << ": Warning: too many threads available." << std::endl
              << "Reducing to a suitable amount" << std::endl;
    while (((stripsize = els/(groups*threads)) < 2) && threads > 1)
      threads--;
    if (threads == 1)
      groups=1;
    stripsize = els/(groups*threads);
  }
  remainder = els-(stripsize*groups*threads);

#if SP_DEBUG > 1
  std::cout << "we have " << threads << " threads available" << std::endl;
  std::cout << "nel1 " << nel1 << std::endl;
  std::cout << "nel2 " << nel2 << std::endl;
  std::cout << "nel3 " << nel3 << std::endl;
  std::cout << "stripsize " << stripsize << std::endl;
  std::cout << "# of strips " << els/stripsize << std::endl;
  std::cout << "remainder " << remainder << std::endl;
#endif
#endif

  result.resize(groups);
  if (groups == 1) {
    result[0].resize(1);
    for (int i=0;i<nel1*nel2*nel3;++i)
      result[0][0].push_back(i);
  } else {
    std::vector< std::vector<int> > stripsizes;
    stripsizes.resize(2);
    stripsizes[0].resize(threads,stripsize);
    stripsizes[1].resize(threads,stripsize);
    int r=0;
    for (int i=0;i<remainder && r < remainder;++i) {
      stripsizes[1][threads-1-i]++;
      r++;
      if (r < remainder) {
        stripsizes[0][threads-1-i]++;
        r++;
      }
    }
    std::vector< std::vector<int> > startelms;
    startelms.resize(2);
    int offs=0;
    for (int i=0;i<threads;++i) {
      startelms[0].push_back(offs*mul);
      offs += stripsizes[0][i];
      startelms[1].push_back(offs*mul);
      offs += stripsizes[1][i];
    }
    for (size_t g=0;g<result.size();++g) { // loop over groups
      result[g].resize(threads);
      for (int t=0;t<threads;++t) { // loop over threads
        int maxx = dir==0?stripsizes[g][t]:nel1;
        int maxy = dir==1?stripsizes[g][t]:nel2;
        int maxz = dir==2?stripsizes[g][t]:nel3;
        for (int i3=0; i3 < maxz; ++i3) {
          for (int i2=0; i2 < maxy; ++i2) { // loop in y direction
            for (int i1=0; i1< maxx; ++i1) {
              int iEl = startelms[g][t]+i1+i2*nel1+i3*nel1*nel2;
              result[g][t].push_back(iEl);
            }
          }
        }
      }
    }
  }

#if defined(USE_OPENMP) && SP_DEBUG > 1
  for (size_t i=0;i<result.size();++i) {
    std::cout << "group " << i << std::endl;
    for (size_t j=0;j<result[i].size();++j) {
      std::cout << "\t thread " << j << " (" << result[i][j].size() << "): ";
      for (size_t k=0;k<result[i][j].size();++k)
        std::cout << result[i][j][k] << " ";
      std::cout << std::endl;
    }
  }
#endif
}


void utl::mapThreadGroups(ThreadGroups& result, const std::vector<int>& map)
{
  for (size_t l=0;l<result.size();++l)
    for (size_t k=0;k<result[l].size();++k)
      for (size_t j=0;j<result[l][k].size();++j)
        result[l][k][j] = map[result[l][k][j]];
}
