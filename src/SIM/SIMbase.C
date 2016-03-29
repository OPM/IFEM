// $Id$
//==============================================================================
//!
//! \file SIMbase.C
//!
//! \date Dec 08 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for NURBS-based FEM simulators.
//!
//==============================================================================

#include "SIMbase.h"
#include "SIMoptions.h"
#include "ASMs2DC1.h"
#include "ASMunstruct.h"
#ifdef HAS_PETSC
#include "SAMpatchPara.h"
#else
#include "SAMpatch.h"
#endif
#include "IntegrandBase.h"
#include "AlgEqSystem.h"
#include "LinSolParams.h"
#include "EigSolver.h"
#include "GlbNorm.h"
#include "ElmNorm.h"
#include "AnaSol.h"
#include "Vec3Oper.h"
#include "Functions.h"
#include "Profiler.h"
#include "Utilities.h"
#include "HDF5Writer.h"
#include "IFEM.h"
#include "tinyxml.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iterator>


bool SIMbase::preserveNOrder  = false;
bool SIMbase::ignoreDirichlet = false;


SIMbase::SIMbase (IntegrandBase* itg) : g2l(&myGlb2Loc)
{
  isRefined = false;
  nsd = 3;
  myProblem = itg;
  mySol = nullptr;
  myEqSys = nullptr;
  mySam = nullptr;
  mySolParams = nullptr;
  nGlPatches = 0;
  nIntGP = nBouGP = 0;

  MPCLess::compareSlaveDofOnly = true; // to avoid multiple slave definitions
}


SIMbase::~SIMbase ()
{
#ifdef SP_DEBUG
  IFEM::cout <<"\nEntering SIMbase destructor"<< std::endl;
#endif

  for (IntegrandMap::iterator it = myInts.begin(); it != myInts.end(); ++it)
    if (it->second != myProblem) delete it->second;

  if (myProblem)   delete myProblem;
  if (mySol)       delete mySol;
  if (myEqSys)     delete myEqSys;
  if (mySam)       delete mySam;
  if (mySolParams) delete mySolParams;

  for (PatchVec::iterator i1 = myModel.begin(); i1 != myModel.end(); i1++)
    delete *i1;

  myModel.clear();
  this->SIMbase::clearProperties();

#ifdef SP_DEBUG
  IFEM::cout <<"Leaving SIMbase destructor"<< std::endl;
#endif
}


void SIMbase::clearProperties ()
{
  for (PatchVec::iterator i1 = myModel.begin(); i1 != myModel.end(); i1++)
    (*i1)->clear(true); // retain the geometry only
  for (SclFuncMap::iterator i2 = myScalars.begin(); i2 != myScalars.end(); i2++)
    delete i2->second;
  for (VecFuncMap::iterator i3 = myVectors.begin(); i3 != myVectors.end(); i3++)
    delete i3->second;
  for (TracFuncMap::iterator i4 = myTracs.begin(); i4 != myTracs.end(); i4++)
    delete i4->second;

  myPatches.clear();
  myGlb2Loc.clear();
  myScalars.clear();
  myVectors.clear();
  myTracs.clear();
  myProps.clear();
  myInts.clear();
}


bool SIMbase::parseGeometryTag (const TiXmlElement* elem)
{
  IFEM::cout <<"  Parsing <"<< elem->Value() <<">"<< std::endl;

  if (!strcasecmp(elem->Value(),"patchfile") && elem->FirstChild()) {
    if (myModel.empty()) {
      const char* file = elem->FirstChild()->Value();
      IFEM::cout <<"\tReading data file "<< file << std::endl;
      std::ifstream isp(file);
      this->readPatches(isp,myModel,"\t");

      if (myModel.empty())
      {
        std::cerr <<" *** SIMbase::parse: No patches read"<< std::endl;
        return false;
      }
    }
  }

  else if (!strcasecmp(elem->Value(),"nodefile") && elem->FirstChild()) {
    if (!this->createFEMmodel()) return false;

    const char* file = elem->FirstChild()->Value();
    if (strstr(file,".gno")) {
      std::ifstream isn(file);
      if (!isn.good()) {
        std::cerr <<" *** SIMbase::read: Failure opening input file "
                  << file << std::endl;
        return false;
      }

      IFEM::cout <<"\tReading data file "<< file << std::endl;
      this->readNodes(isn);
    }
    else if (strstr(file,".hdf5")) {
      IFEM::cout <<"\tReading global node numbers from "<< file << std::endl;
      HDF5Writer hdf5(file,ProcessAdm(),true,true);
      const char* field = elem->Attribute("field");
      for (int i = 1; i <= nGlPatches; i++)
      {
        std::vector<int> nodes;
        ASMbase* pch = this->getPatch(this->getLocalPatchIndex(i));
        if (pch && hdf5.readVector(0, field ? field : "node numbers", i, nodes))
          pch->setNodeNumbers(nodes);
      }
      hdf5.closeFile(0, true);
    }
  }

  else if (!strcasecmp(elem->Value(),"partitioning")) {
    int proc = 0;
    if (!utl::getAttribute(elem,"procs",proc))
      return false;
    else if (proc != adm.getNoProcs()) // silently ignore
      return true;
    IFEM::cout <<"\tNumber of partitions: "<< proc << std::endl;

    const TiXmlElement* part = elem->FirstChildElement("part");
    for (; part; part = part->NextSiblingElement("part")) {
      int first = -2, last = -2;
      utl::getAttribute(part,"proc",proc);
      utl::getAttribute(part,"lower",first);
      utl::getAttribute(part,"upper",last);
      if (last > nGlPatches) nGlPatches = last;
      if (proc == adm.getProcId()) {
        myPatches.reserve(last-first+1);
        for (int j = first; j <= last && j > -1; j++)
          myPatches.push_back(j);
      }
    }

    // If equal number of blocks per processor
    if (myPatches.empty() && utl::getAttribute(elem,"nperproc",proc)) {
      for (int j = 1; j <= proc; j++)
        myPatches.push_back(adm.getProcId()*proc+j);
      nGlPatches = adm.getNoProcs()*proc;
    }
  }

  else if (!strcasecmp(elem->Value(),"topologysets")) {
    std::string name, type;
    const TiXmlElement* set = elem->FirstChildElement("set");
    for (; set; set = set->NextSiblingElement("set"))
      if (utl::getAttribute(set,"name",name)) {
        int idim = 3;
        utl::getAttribute(set,"type",type,true);
        if (type == "volume")
          idim = 3;
        else if (type == "face" || type == "surface")
          idim = 2;
        else if (type == "edge" || type == "curve")
          idim = 1;
        else if (type == "vertex" || type == "point")
          idim = 0;
        else
          utl::getAttribute(set,"dimension",idim);
        if (idim > 0 && utl::getAttribute(set,"closure",type,true))
          if (type == "open") idim = -idim; // i.e. excluding its boundary

        TopEntity& top = myEntitys[name];
        const TiXmlElement* item = set->FirstChildElement("item");
        for (; item; item = item->NextSiblingElement("item")) {
          int patch = 1;
          utl::getAttribute(item,"patch",patch);
          if ((patch = this->getLocalPatchIndex(patch)) > 0) {
            if (abs(idim) == (int)this->getNoParamDim())
              top.insert(TopItem(patch,0,idim));
            else if (item->FirstChild())
            {
              std::string value(item->FirstChild()->Value());
              char* cval = strtok(const_cast<char*>(value.c_str())," ");
              for (; cval; cval = strtok(nullptr," "))
                top.insert(TopItem(patch,atoi(cval),idim));
            }
          }
        }
      }

    if (!myEntitys.empty()) {
      IFEM::cout <<"\tTopology sets: ";
      TopologySet::const_iterator it;
      for (it = myEntitys.begin(); it != myEntitys.end(); it++) {
        if (it != myEntitys.begin()) IFEM::cout <<"\t               ";
        IFEM::cout << it->first;
        for (const auto& it2 : it->second)
          IFEM::cout << it2;
        IFEM::cout << std::endl;
      }
    }
  }

  return true;
}


//! \brief Integer value flagging global axes in boundary conditions.
#define GLOBAL_AXES     -1
//! \brief Integer value flagging local axes in boundary conditions.
#define LOCAL_AXES      -2
//! \brief Integer value flagging local projected axes in boundary conditions.
#define LOCAL_PROJECTED -3

bool SIMbase::parseBCTag (const TiXmlElement* elem)
{
  IFEM::cout <<"  Parsing <"<< elem->Value() <<">"<< std::endl;

  if (!strcasecmp(elem->Value(),"propertyfile") && elem->FirstChild()) {
    const char* file = elem->FirstChild()->Value();
    std::ifstream isp(file);
    if (!isp)
    {
      std::cerr <<" *** SIMbase::parseBCTag: Failure opening input file "
                << file << std::endl;
      return false;
    }
    IFEM::cout <<"\tReading data file "<< file << std::endl;
    while (isp.good())
    {
      Property p;
      int ldim, lindx = 0;
      isp >> p.pindx >> p.patch >> ldim;
      if (ldim < (int)this->getNoParamDim()) isp >> lindx;

      // We always require the item indices to be 1-based
      p.ldim = ldim;
      p.lindx = 1+lindx;
      p.patch = this->getLocalPatchIndex(1+p.patch);
      if (p.patch > 0 && isp.good())
        myProps.push_back(p);
    }
  }

  else if (!strcasecmp(elem->Value(),"propertycodes")) {
    const TiXmlElement* code = elem->FirstChildElement("code");
    for (; code; code = code->NextSiblingElement()) {
      int icode = 0;
      utl::getAttribute(code,"value",icode);
      const TiXmlElement* patch = code->FirstChildElement("patch");
      for (; patch; patch = patch->NextSiblingElement("patch")) {
        Property p;
        int ival = 0;
        p.pindx = icode;
        if (utl::getAttribute(patch,"face",ival))
          p.ldim = 2;
        else if (utl::getAttribute(patch,"edge",ival))
          p.ldim = 1;
        else if (utl::getAttribute(patch,"vertex",ival))
          p.ldim = 0;
        if (ival > 0)
          p.lindx = ival;
        if (utl::getAttribute(patch,"index",ival) && ival > 0)
          if (p.lindx > 0 && (p.patch = this->getLocalPatchIndex(ival)) > 0)
            myProps.push_back(p);
      }
    }
  }

  else if (!strcasecmp(elem->Value(),"neumann")) {
    std::string set, type;
    utl::getAttribute(elem,"set",set);
    utl::getAttribute(elem,"type",type,true);
    int code = this->getUniquePropertyCode(set);
    if (code == 0) utl::getAttribute(elem,"code",code);
    IFEM::cout <<"\tNeumann code "<< code;
    if (type == "anasol") {
      IFEM::cout <<" (analytic)";
      this->setPropertyType(code,Property::NEUMANN_ANASOL);
    }
    else if (type == "generic") {
      IFEM::cout <<" (generic)";
      if (elem->FirstChild()) {
        int ndir = 0;
        utl::getAttribute(elem,"direction",ndir);
        if (ndir > 0) IFEM::cout <<" direction "<< ndir;
        this->setPropertyType(code,Property::ROBIN);
        this->setNeumann(elem->FirstChild()->Value(),"expression",ndir,code);
      }
      else
        this->setPropertyType(code,Property::NEUMANN_GENERIC);
    }
    else {
      int ndir = 0;
      utl::getAttribute(elem,"direction",ndir);
      IFEM::cout <<" direction "<< ndir;
      if (!type.empty()) IFEM::cout <<" ("<< type <<")";
      if (elem->FirstChild())
        this->setNeumann(elem->FirstChild()->Value(),type,ndir,code);
      else
        this->setNeumann("0.0",type,ndir,code);
    }
    IFEM::cout << std::endl;
  }

  else if (!strcasecmp(elem->Value(),"dirichlet") && !ignoreDirichlet) {
    int comp = 0;
    int basis = 1;
    std::string set, type, axes;
    // long and short form supported, short prioritized
    utl::getAttribute(elem,"component",comp);
    utl::getAttribute(elem,"comp",comp);
    utl::getAttribute(elem,"set",set);
    utl::getAttribute(elem,"type",type,true);
    utl::getAttribute(elem,"axes",axes,true);
    utl::getAttribute(elem,"basis",basis);
    int code = this->getUniquePropertyCode(set,comp);
    if (code == 0) utl::getAttribute(elem,"code",code);
    if (axes == "local projected")
      comp = LOCAL_PROJECTED;
    else if (axes == "local")
      comp = LOCAL_AXES;
    else
      comp = GLOBAL_AXES;
    if (comp <= LOCAL_AXES)
    {
      // Check for definition of first tangent direction (relevant for 3D only)
      utl::getAttribute(elem,"t1",axes,true);
      if (axes[0] >= 'x' && axes[0] <= 'z')
        comp -= 10*(axes[0]-'w');
    }
    const TiXmlNode* dval = elem->FirstChild();
    if (type == "anasol") {
      this->setPropertyType(code,Property::DIRICHLET_ANASOL,-1,basis);
      IFEM::cout <<"\tDirichlet code "<< code <<": (analytic)"<< std::endl;
    }
    // this is a horrible hack
    else if (!dval || (type != "expression" &&
                       type != "field" && atof(dval->Value()) == 0.0)) {
      this->setPropertyType(code,Property::DIRICHLET,comp,basis);
      IFEM::cout <<"\tDirichlet code "<< code <<": (fixed)"<< std::endl;
    }
    else if (dval) {
      this->setPropertyType(code,Property::DIRICHLET_INHOM,comp,basis);
      IFEM::cout <<"\tDirichlet code "<< code;
      if (!type.empty()) IFEM::cout <<" ("<< type <<")";
      myScalars[abs(code)] = utl::parseRealFunc(dval->Value(),type);
      IFEM::cout << std::endl;
    }
  }

  return true;
}


static bool noDumpDataYet = true; //!< To read only once in adaptive loops


bool SIMbase::parseOutputTag (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"dump_lhs_matrix") &&
      strcasecmp(elem->Value(),"dump_rhs_vector") &&
      strcasecmp(elem->Value(),"dump_sol_vector"))
    return opt.parseOutputTag(elem);

  if (elem->FirstChild() && noDumpDataYet)
  {
    DumpData dmp;
    std::string format;
    utl::getAttribute(elem,"format",format);
    utl::getAttribute(elem,"step",dmp.step);
    dmp.format = format[0];
    dmp.fname = elem->FirstChild()->Value();
    if (toupper(elem->Value()[5]) == 'R')
      rhsDump.push_back(dmp);
    else if (toupper(elem->Value()[5]) == 'L')
      lhsDump.push_back(dmp);
    else
      solDump.push_back(dmp);
  }

  return true;
}


bool SIMbase::parse (const TiXmlElement* elem)
{
  bool result = true;
  if (!strcasecmp(elem->Value(),"discretization"))
    result = opt.parseDiscretizationTag(elem);
  else if (!strcasecmp(elem->Value(),"linearsolver"))
  {
    std::string solver;
    if (utl::getAttribute(elem,"class",solver,true))
      opt.setLinearSolver(solver);
  }
  else if (!strcasecmp(elem->Value(),"eigensolver"))
    utl::getAttribute(elem,"mode",opt.eig);
  else if (!strcasecmp(elem->Value(),"postprocessing"))
    noDumpDataYet = lhsDump.empty() && rhsDump.empty();
  else if (!strcasecmp(elem->Value(),"initialcondition"))
  {
    std::string field, file;
    if (utl::getAttribute(elem,"field",field))
    {
      ICInfo info(field);
      std::string type("file");
      utl::getAttribute(elem,"type",type);
      if (type == "file") {
        utl::getAttribute(elem,"file",file);
        utl::getAttribute(elem,"file_field",info.file_field);
        utl::getAttribute(elem,"file_level",info.file_level);
        utl::getAttribute(elem,"geo_level",info.geo_level);
      }
      else { // function
        // long and short from supported, short prioritized
        utl::getAttribute(elem,"component",info.component);
        utl::getAttribute(elem,"comp",info.component);
        info.file_field = type;
        info.function = utl::getValue(elem,"initialcondition");
        file = "nofile";
      }
      utl::getAttribute(elem,"level",info.sim_level);
      int basis = 1;
      utl::getAttribute(elem,"basis",basis);
      info.basis = basis;

      IFEM::cout <<"\tInitial condition";
      if (info.component > -1)
        IFEM::cout <<" function: "<< info.function;
      else
        IFEM::cout <<" file: "<< file;

      IFEM::cout <<"\n\tField name: \""<< info.sim_field;
      if (info.component == -1)
        IFEM::cout <<"\" (on file \""<< info.file_field <<"\")";
      else
        IFEM::cout <<" (component "<< info.component <<" basis "<< basis <<")";

      IFEM::cout <<"\n\tTime level: "<< info.sim_level;
      if (info.component == -1)
        IFEM::cout <<" (on file "<< info.file_level
                   <<" with basis "<< info.geo_level <<")";

      IFEM::cout << std::endl;

      myICs[file].push_back(info);
    }
    else
      result = false;
  }

  // Create the default geometry of no patchfile is specified
  if (myModel.empty() && !strcasecmp(elem->Value(),"geometry"))
    if (this->getNoParamDim() > 0 && !elem->FirstChildElement("patchfile"))
    {
      IFEM::cout <<"  Using default linear geometry basis on unit domain [0,1]";
      if (this->getNoParamDim() > 1) IFEM::cout <<"^"<< this->getNoParamDim();
      IFEM::cout << std::endl;
      myModel.resize(1,this->createDefaultGeometry(elem));
    }

  if (!strcasecmp(elem->Value(),"linearsolver")) {
    if (!mySolParams)
      mySolParams = new LinSolParams(nsd);
    result &= mySolParams->read(elem);
  }

  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if (!strcasecmp(elem->Value(),"geometry"))
      result &= this->parseGeometryTag(child);
    else if (!strcasecmp(elem->Value(),"boundaryconditions"))
      result &= this->parseBCTag(child);
    else if (!strcasecmp(elem->Value(),"linearsolver"))
      result &= this->parseLinSolTag(child);
    else if (!strcasecmp(elem->Value(),"eigensolver"))
      result &= opt.parseEigSolTag(child);
    else if (!strcasecmp(elem->Value(),"postprocessing"))
      result &= this->parseOutputTag(child);
    else if (!strcasecmp(elem->Value(),"discretization"))
      result &= opt.parseDiscretizationTag(child);

  return result;
}


int SIMbase::parseMaterialSet (const TiXmlElement* elem, int mindex)
{
  std::string set;
  utl::getAttribute(elem,"set",set);
  int code = this->getUniquePropertyCode(set,0);

  if (code == 0)
    utl::getAttribute(elem,"code",code);

  if (code > 0)
    this->setPropertyType(code,Property::MATERIAL,mindex);

  return code;
}


const char** SIMbase::getPrioritizedTags () const
{
  // Tags to be parsed first, and in the order specified
  static const char* special[] = { "discretization", "geometry", 0 };
  return special;
}


bool SIMbase::parse (char* keyWord, std::istream& is)
{
  char* cline = 0;
  if (!strncasecmp(keyWord,"PATCHES",7))
  {
    int npatch = atoi(keyWord+7);
    if (myModel.empty())
    {
      IFEM::cout <<"\nNumber of patches: "<< npatch << std::endl;
      for (int i = 0; i < npatch && (cline = utl::readLine(is)); i++)
      {
	std::ifstream is(cline);
	if (is.good())
	{
	  IFEM::cout <<"\nReading patch file "<< cline << std::endl;
	  ASMbase* pch = this->readPatch(is,i);
	  if (pch)
	    myModel.push_back(pch);
	}
	else
	  std::cerr <<" *** SIMbase: Failure opening patch file"
		    << cline << std::endl;
      }

      if ((int)myModel.size() < npatch)
      {
	std::cerr <<" *** SIMbase::parse: Expected "<< npatch
		  <<" patches but could read only "<< myModel.size()
		  << std::endl;
	return false;
      }
    }
    else // just read through the npatch next lines without doing anything
      for (int i = 0; i < npatch && utl::readLine(is); i++);
  }

  else if (!strncasecmp(keyWord,"PATCHFILE",9))
  {
    if (myModel.empty())
    {
      size_t i = 9; while (i < strlen(keyWord) && isspace(keyWord[i])) i++;
      IFEM::cout <<"\nReading data file "<< keyWord+i << std::endl;
      std::ifstream isp(keyWord+i);
      this->readPatches(isp,myModel);

      if (myModel.empty())
      {
	std::cerr <<" *** SIMbase::parse: No patches read"<< std::endl;
	return false;
      }
    }
  }

  else if (!strncasecmp(keyWord,"NODEFILE",8))
  {
    if (!this->createFEMmodel()) return false;

    bool oneBasedIdx = keyWord[8] == '1';
    size_t i = (oneBasedIdx || keyWord[8] == '0') ? 9 : 8;
    while (i < strlen(keyWord) && isspace(keyWord[i])) i++;

    std::stringstream filenames(keyWord+i);
    for (size_t basis = 0; filenames.good(); basis++)
    {
      std::string filename;
      filenames >> filename;
      if (basis == 0 && filenames.good()) basis = 1;

      std::ifstream isn(filename.c_str());
      if (!isn.good())
      {
	std::cerr <<" *** SIMbase::read: Failure opening input file "
		  << filename << std::endl;
	return false;
      }

      if (basis < 2) IFEM::cout <<"\n";
      IFEM::cout <<"Reading data file "<< filename;
      if (basis > 0) IFEM::cout <<" (basis "<< basis <<")";
      IFEM::cout << std::endl;

      while (isn.good())
      {
	int patch = 0;
	isn >> patch;
	if (!oneBasedIdx) ++patch;
	int pid = this->getLocalPatchIndex(patch);
	if (pid < 0) return false;

	if (!this->readNodes(isn,pid-1,basis,oneBasedIdx))
	{
	  std::cerr <<" *** SIMbase::parse: Failed to assign node numbers"
		    <<" for patch "<< patch << std::endl;
	  return false;
	}
      }
    }
  }

  else if (!strncasecmp(keyWord,"PROPERTYFILE",12))
  {
    bool oneBasedIdx = keyWord[12] == '1';
    size_t i = (oneBasedIdx || keyWord[12] == '0') ? 13 : 12;
    while (i < strlen(keyWord) && isspace(keyWord[i])) i++;

    std::ifstream isp(keyWord+i);
    if (!isp)
    {
      std::cerr <<" *** SIMbase::read: Failure opening input file "
		<< std::string(keyWord+i) << std::endl;
      return false;
    }

    IFEM::cout <<"\nReading data file "<< keyWord+i << std::endl;
    while (isp.good())
    {
      Property p;
      int ldim, lindx = 0;
      isp >> p.pindx >> p.patch >> ldim;
      if (ldim < (int)this->getNoParamDim()) isp >> lindx;

      if (!oneBasedIdx)
      {
	// We always require the item indices to be 1-based
	++p.patch;
	++lindx;
      }

      p.ldim = ldim;
      p.lindx = lindx;
      p.patch = this->getLocalPatchIndex(p.patch);
      if (p.patch > 0 && isp.good())
	myProps.push_back(p);
    }
  }

  else if (!strncasecmp(keyWord,"DIRICHLET",9))
  {
    if (ignoreDirichlet) return true; // Ignore all boundary conditions

    int ndir = atoi(keyWord+9);
    IFEM::cout <<"\nNumber of Dirichlet properties: "<< ndir << std::endl;
    for (int i = 0; i < ndir && (cline = utl::readLine(is)); i++)
    {
      int code = atoi(strtok(cline," "));
      double d = (cline = strtok(nullptr," ")) ? atof(cline) : 0.0;
      IFEM::cout <<"\tDirichlet code "<< code <<": ";
      if (d == 0.0)
      {
	this->setPropertyType(code,Property::DIRICHLET);
	IFEM::cout <<"(fixed)";
      }
      else if (code > 0)
      {
	this->setPropertyType(-code,Property::DIRICHLET_INHOM);

	cline = strtok(nullptr," ");
	myScalars[code] = const_cast<RealFunc*>(utl::parseRealFunc(cline,d));
      }
      else
	IFEM::cout <<"(ignored)";
      IFEM::cout << std::endl;
    }
  }

  else if (!strncasecmp(keyWord,"LINEARSOLVER",12))
    this->readLinSolParams(is,atoi(keyWord+12));

  else if (!strncasecmp(keyWord,"PARTITIONING",12))
  {
    int nproc = atoi(keyWord+12);
    IFEM::cout <<"\nNumber of partitions: "<< nproc << std::endl;

    nGlPatches = 0;
    for (int i = 0; i < nproc && (cline = utl::readLine(is)); i++)
    {
      int proc  = atoi(strtok(cline," "));
      int first = atoi(strtok(nullptr," "));
      int last  = atoi(strtok(nullptr," "));

      if (last > nGlPatches) nGlPatches = last;

      if (proc == myPid && last >= first)
      {
	myPatches.reserve(last-first+1);
	for (int j = first; j <= last; j++)
	  myPatches.push_back(j);
      }
    }

    if (myPatches.empty())
    {
      std::cerr <<" *** SIMbase::parse: No partitioning input for processor: "
		<< myPid << std::endl;
      return false;
    }
  }

  else if (!strncasecmp(keyWord,"VISUALIZATION",13))
  {
    cline = utl::readLine(is);
    if ((cline = strtok(cline," "))) opt.format = atoi(cline);
    if ((cline = strtok(nullptr," "))) opt.saveInc = atoi(cline);
  }

#ifdef SP_DEBUG
  // Since the same input file might be parsed by several substep solvers,
  // warnings on ignored keywords are issued when compiled in debug mode only.
  else if (isalpha(keyWord[0]))
    std::cerr <<" *** SIMbase::parse: Unknown keyword: "<< keyWord << std::endl;
#endif

  return true;
}


void SIMbase::readLinSolParams (std::istream& is, int npar)
{
  if (!mySolParams)
    mySolParams = new LinSolParams(nsd);

  mySolParams->read(is,npar);
}


bool SIMbase::parseLinSolTag (const TiXmlElement* elem)
{
  if (!strcasecmp(elem->Value(),"class"))
    if (elem->FirstChild())
      opt.setLinearSolver(elem->FirstChild()->Value());

  return true;
}


bool SIMbase::createFEMmodel (char resetNumb)
{
  if (resetNumb)
  {
    ASMstruct::resetNumbering();
    ASMunstruct::resetNumbering();
  }

  for (size_t i = 0; i < myModel.size(); i++)
  {
    myModel[i]->setGauss(opt.nGauss[0]); // in the case of immersed boundaries,
    // the number of Gauss quadrature points must be known at this point

    if (myModel[i]->isShared() && myModel[i]->hasXNodes())
    {
      // This patch shares its FE data with another patch, but has been assigned
      // additional nodes due to constraints in local coordinate systems,
      // Lagrange multipliers, etc. Need therefore to unshare the FE data.
      ASMbase* newPatch = myModel[i]->cloneUnShared();
      if (newPatch)
      {
        IFEM::cout <<"\tNote: Unsharing FE data of P"<< 1+i << std::endl;
        newPatch->setGlobalNodeNums(myModel[i]->getMyNodeNums());
        delete myModel[i];
        myModel[i] = newPatch;
      }
    }
    else if (!myModel[i]->generateFEMTopology())
      return false;
    else if (myModel[i]->isShared() && resetNumb == 'Y')
      myModel[i]->setGlobalNodeNums(IntVec());
  }

  if (nGlPatches == 0 && (!adm.isParallel() || adm.getNoProcs() == 1))
    nGlPatches = myModel.size();

#ifdef HAS_PETSC
  // When PETSc is used, we need to retain all DOFs in the equation system.
  // The fixed DOFs (if any) will receive a homogeneous constraint instead.
  ASMbase::fixHomogeneousDirichlet = opt.solver != SystemMatrix::PETSC;
#else
  ASMbase::fixHomogeneousDirichlet = true;
#endif

  return true;
}


int SIMbase::getLocalPatchIndex (int patchNo) const
{
  if (patchNo < 1 || (patchNo > nGlPatches && nGlPatches > 0))
  {
    std::cerr <<" *** SIMbase::getLocalPatchIndex: Patch number "<< patchNo
	      <<" out of range [1,"<< nGlPatches <<"]"<< std::endl;
    return -1;
  }
  else if (myPatches.empty() || nProc == 1)
    return patchNo;

  for (size_t i = 0; i < myPatches.size(); i++)
    if (myPatches[i] == patchNo)
      return 1+i;

  return 0;
}


ASMbase* SIMbase::getPatch (size_t idx) const
{
  return idx > 0 && idx <= myModel.size() ? myModel[idx-1] : nullptr;
}


bool SIMbase::preprocess (const IntVec& ignored, bool fixDup)
{
  if (myModel.empty())
    return true; // Empty simulator, nothing to preprocess

  if (mySam && !isRefined)
  {
    std::cerr <<" *** SIMbase::preprocess: Logic error, invoked more than once"
              <<" for "<< (myHeading.empty() ? this->getName() : myHeading)
              << std::endl;
    return false;
  }

  if (nProc > 1 && myPatches.empty() && adm.isParallel())
  {
    std::cerr <<" *** SIMbase::preprocess: No partitioning information for "
              << nProc <<" processors found."<< std::endl;
    return false;
  }

  static int substep = 10;
  this->printHeading(substep);

  // Perform some sub-class specific pre-preprocessing, if any
  this->preprocessA();

  // Create the classical FE data structures
  if (!this->createFEMmodel('Y')) return false;

  PatchVec::const_iterator mit;
  IntVec::const_iterator it;
  ASMbase* pch;
  size_t patch;

  // Erase all patches that should be ignored in the analysis
  for (it = ignored.begin(); it != ignored.end(); ++it)
    if ((pch = this->getPatch(*it)))
      pch->clear();

  // If material properties are specified for at least one patch, assign the
  // property code 999999 to all patches with no material property code yet
  PatchVec pchWthMat;
  for (PropertyVec::const_iterator p = myProps.begin(); p != myProps.end(); ++p)
    if (p->pcode == Property::MATERIAL && (pch = this->getPatch(p->patch)))
      if (!pch->empty()) pchWthMat.push_back(pch);

  if (!pchWthMat.empty())
    for (mit = myModel.begin(), patch = 1; mit != myModel.end(); ++mit, patch++)
      if (std::find(pchWthMat.begin(),pchWthMat.end(),*mit) == pchWthMat.end())
	myProps.push_back(Property(Property::MATERIAL,999999,patch,
				   (*mit)->getNoParamDim()));

  if (fixDup)
  {
    // Check for duplicated nodes (missing topology)
    int nDupl = 0;
    std::map<Vec3,int> globalNodes;
    for (mit = myModel.begin(), patch = 1; mit != myModel.end(); ++mit, patch++)
      if (!(*mit)->empty())
      {
	IFEM::cout <<"   * Checking Patch "<< patch << std::endl;
	for (size_t node = 1; node <= (*mit)->getNoNodes(); node++)
	{
	  Vec3 X((*mit)->getCoord(node));
	  std::map<Vec3,int>::const_iterator xit = globalNodes.find(X);
	  if (xit == globalNodes.end())
	    globalNodes.insert(std::make_pair(X,(*mit)->getNodeID(node)));
	  else if ((*mit)->mergeNodes(node,xit->second))
	    nDupl++;
	}
      }
    if (nDupl > 0)
      IFEM::cout <<"   * "<< nDupl <<" duplicated nodes merged."<< std::endl;
  }

#if SP_DEBUG > 2
  typedef std::pair<int,int> Ipair;
  typedef std::vector<Ipair> Ipairs;
  std::map<int,Ipairs> nodeInfo;
  for (mit = myModel.begin(), patch = 1; mit != myModel.end(); ++mit, patch++)
    if (!(*mit)->empty())
      for (size_t n = 1; n <= (*mit)->getNoNodes(); n++)
        nodeInfo[(*mit)->getNodeID(n)].push_back(std::make_pair(patch,n));

  std::map<int,Ipairs>::const_iterator nit;
  for (nit = nodeInfo.begin(); nit != nodeInfo.end(); ++nit)
    if (nit->second.size() > 1)
    {
      std::cout <<"\nConnectivity for node "<< nit->first <<":";
      for (size_t n = 0; n < nit->second.size(); n++)
        std::cout <<" P"<< nit->second[n].first <<","<< nit->second[n].second;
    }
  std::cout << std::endl;
#endif

  // Renumber the nodes to account for overlapping nodes and erased patches.
  // In parallel simulations, the resulting global-to-local node number mapping
  // will map the global node numbers to local node numbers on the current
  // processor. In serial simulations, the global-to-local mapping will be unity
  // unless the original global node number sequence had "holes" due to
  // duplicated nodes and/or erased patches.
  int ngnod = 0;
  int renum = 0;
  if (preserveNOrder)
  {
    renum = ASMbase::renumberNodes(myModel,myGlb2Loc);
    ngnod = g2l->size();
  }
  else for (mit = myModel.begin(); mit != myModel.end(); ++mit)
    renum += (*mit)->renumberNodes(myGlb2Loc,ngnod);

  if (renum > 0)
    IFEM::cout <<"\nRenumbered "<< renum <<" nodes."<< std::endl;

  for (mit = myModel.begin(); mit != myModel.end(); ++mit)
    (*mit)->renumberNodes(*g2l);
  ASMs2DC1::renumberNodes(*g2l);

  // Perform specialized preprocessing before the assembly initialization.
  // This typically involves the system-level Lagrange multipliers, etc.
  if (!this->preprocessBeforeAsmInit(ngnod))
    return false;

  IFEM::cout <<"\nResolving Dirichlet boundary conditions"<< std::endl;
  ASMstruct::resetNumbering(ngnod); // to account for possibly added nodes

  // Process the Dirichlet boundary conditions in the order of increasing
  // dimension, such that vertex definitions override definitions on edges,
  // and edge definitions override definitions on faces
  size_t nprop = 0;
  int code, dofs, ierr = 0, iwar = 0;
  PropertyVec::const_iterator q;
  for (unsigned char dim = 0; nprop < myProps.size(); dim++)
    for (q = myProps.begin(); q != myProps.end(); ++q)
      if (abs(q->ldim) == dim || (dim == 2 && abs(q->ldim) > 3))
      {
        nprop++;
        code = q->pindx;
        dofs = abs(code%1000000);
        switch (q->pcode) {
        case Property::DIRICHLET:
          code = 0;
        case Property::DIRICHLET_INHOM:
          break;

        case Property::UNDEFINED:
          ++iwar;
#ifdef SP_DEBUG
          std::cout <<"  ** SIMbase::preprocess: Undefined property set, code="
                    << q->pindx <<" Patch="<< q->patch <<" Item="
                    << (int)q->lindx <<" "<< (int)q->ldim <<"D"<< std::endl;
#endif
        default:
          dofs = 0;
          break;
        }

        if (dofs > 0)
          if (this->addConstraint(q->patch,q->lindx,q->ldim,dofs,code,ngnod,q->basis))
            IFEM::cout << std::endl;
          else
            ++ierr;
      }

  if (iwar > 0)
    std::cerr <<"\n  ** SIMbase::preprocess: Warning: "<< iwar
              <<" undefined property sets were detected.\n";
  if (ierr > 0)
    std::cerr <<"\n *** SIMbase::preprocess: Error: "<< ierr
              <<" invalid Dirichlet properties were detected.\n"
              <<"     Please check your model, execution aborts..."
              << std::endl;
  else if (iwar > 0)
    std::cerr <<"     Please verify your model, execution continues..."
              << std::endl;

  // Compute the set of all MPCs over the whole model. This will also merge
  // multiple constraint equations defined on interfaces between patches.
  MPCSet allMPCs;
  ASMbase::mergeAndGetAllMPCs(myModel,allMPCs);

  // Set initial values for the inhomogeneous dirichlet conditions, if any,
  // and compute coupling coefficients for the C1-continuity constraints
  if (!this->initDirichlet()) return false;

  // Resolve possibly chaining of the MPC equations
  if (!allMPCs.empty())
    ASMbase::resolveMPCchains(allMPCs,this->hasTimeDependentDirichlet());

  // Generate element groups for multi-threading
  bool silence = msgLevel < 1 || (msgLevel < 2 && myModel.size() > 1);
  for (mit = myModel.begin(); mit != myModel.end() && myProblem; ++mit)
    if (!(*mit)->empty())
      (*mit)->generateThreadGroups(*myProblem,silence);

  for (q = myProps.begin(); q != myProps.end(); ++q)
    if (q->pcode == Property::NEUMANN ||
        q->pcode == Property::NEUMANN_GENERIC ||
        q->pcode == Property::ROBIN)
      this->generateThreadGroups(*q,silence);

  // Preprocess the result points
  this->preprocessResultPoints();

  // Initialize data structures for the algebraic system
  if (mySam) delete mySam;
#ifdef HAS_PETSC
  if (opt.solver == SystemMatrix::PETSC)
    mySam = new SAMpatchPara(*g2l,adm);
  else
    mySam = new SAMpatch();
#else
  mySam = new SAMpatch();
#endif
  if (!static_cast<SAMpatch*>(mySam)->init(myModel,ngnod))
    return false;

  if (!myProblem)
  {
    std::cerr <<"\n *** SIMbase::preprocess(): No problem integrand for the "
              << this->getName() <<"-simulator."<< std::endl;
    return false;
  }

  if (opt.solver == SystemMatrix::PETSC)
    for (mit = myModel.begin(); mit != myModel.end() && myProblem; ++mit)
      if ((*mit)->end_BC() != (*mit)->begin_BC())
      {
        std::cerr <<"\n *** SIMbase::preprocess: Patch "<< (*mit)->idx+1
                  <<"  has boundary condition codes (eliminated DOFs).\n"
                  <<"       This is not supported when PETSc is used,"
                  <<" the input file must be corrected."<< std::endl;
        return false;
      }

  // Now perform the sub-class specific final preprocessing, if any
  return this->preprocessB() && ierr == 0;
}


void SIMbase::generateThreadGroups (const Property& p, bool silence)
{
  ASMbase* pch = this->getPatch(p.patch);
  if (pch && abs(p.ldim)+1 == pch->getNoParamDim())
    pch->generateThreadGroups(p.lindx,silence);
}


/*!
  \brief A helper class used by SIMbase::getUniquePropertyCode.
  \details The class is just an unary function that checks whether a Property
  object has a specified integer code.
*/

class hasCode
{
  int myCode; //!< The property code to compare with

public:
  //! \brief Constructor initializing the property code to search for.
  hasCode(int code) : myCode(abs(code)) {}
  //! \brief Returns \e true if the Property \a p has the code \a myCode
  bool operator()(const Property& p) { return abs(p.pindx) == myCode; }
};


int SIMbase::getUniquePropertyCode (const std::string& setName, int code)
{
  if (setName.empty()) return 0;

  int cinc = code < 0 ? -1000000 : 1000000;
  if (code == 0) code = cinc;
  PropertyVec::const_iterator pit = myProps.begin();
  for (int trial = 0; pit != myProps.end(); trial++)
  {
    if (trial > 0) code += cinc;
    pit = std::find_if(myProps.begin(),myProps.end(),hasCode(code));
  }

  return this->createPropertySet(setName,code) ? code : 0;
}


bool SIMbase::createPropertySet (const std::string& setName, int pc)
{
  if (setName.empty()) return true;

  TopologySet::const_iterator tit = myEntitys.find(setName);
  if (tit == myEntitys.end())
  {
    std::cerr <<" *** SIMbase::createPropertySet: Undefined topology set \""
              << setName <<"\""<< std::endl;
    return false;
  }

  if (std::find_if(myProps.begin(),myProps.end(),hasCode(pc)) != myProps.end())
  {
    std::cerr <<" *** SIMbase::createPropertySet: Duplicated property code "
              << pc << std::endl;
    return false;
  }

  // Create the actual property objects that are used during simulation
  TopEntity::const_iterator top;
  for (top = tit->second.begin(); top != tit->second.end(); ++top)
    myProps.push_back(Property(Property::UNDEFINED,pc,
                               top->patch,top->idim,top->item));

  return true;
}


/*!
  Negative values on \a pindex are used to flag the use of local axes for
  Dirichlet boundary conditions as follows:
  \a pindex=-2 : Compute local axes directions directly at the Greville points.
  \a pindex=-3 : Compute local axes directions by projecting the local tangent
  direction onto the spline basis, to obtain directions at control points.
*/

size_t SIMbase::setPropertyType (int code, Property::Type ptype,
                                 int pindex, char basis)
{
  size_t nDefined = 0;
  for (PropertyVec::iterator p = myProps.begin(); p != myProps.end(); ++p)
    if (abs(p->pindx) == abs(code) && p->pcode == Property::UNDEFINED)
      if (p->patch > 0 && p->patch <= myModel.size())
      {
        ++nDefined;
        p->pcode = ptype;
        p->basis = basis;

        if (ptype == Property::MATERIAL && pindex >= 0)
          p->pindx = pindex; // Index to material property container
        else if (ptype == Property::NEUMANN_ANASOL ||
                 ptype == Property::DIRICHLET_ANASOL)
        {
          // Let all analytical boundary condition properties receive the same
          // property code, because there can only be one analytical solution
          for (PropertyVec::iterator q = myProps.begin(); q != p; ++q)
            if (ptype == q->pcode)
            {
              p->pindx = abs(q->pindx);
              break;
            }
        }
        else if (ptype >= Property::DIRICHLET && pindex <= LOCAL_AXES)
        {
          p->lindx *= -1; // flag the use of local axis directions
          if (abs(p->ldim) == 2 && pindex < 10)
          {
            // Flag first local tangent direction ['x','z']
            p->ldim = (p->ldim/2)*('w'-pindex/10);
            pindex %= 10;
          }
          if (pindex == LOCAL_PROJECTED)
            p->lindx -= 10; // enable projection of the local axes definitions
          preserveNOrder = true; // because extra nodes might be added
        }

        if (p->ldim != 0 && p->pindx > 0 && code < 0) // flag direct evaluation
          if (ptype >= Property::DIRICHLET_INHOM)
            p->pindx = -p->pindx;
      }

  return nDefined;
}


size_t SIMbase::setVecProperty (int code, Property::Type ptype,
                                VecFunc* field, int pflag)
{
  if (field) myVectors[abs(code)] = field;
  return this->setPropertyType(code,ptype,pflag);
}


bool SIMbase::setTracProperty (int code, Property::Type ptype,
                               TractionFunc* field)
{
  if (field) myTracs[code] = field;
  return this->setPropertyType(code,ptype);
}


bool SIMbase::setNeumann (const std::string& prop, const std::string& type,
			  int direction, int code)
{
  if (direction == 0 && this->getNoFields() == 1)
  {
    RealFunc* f = utl::parseRealFunc(prop,type);
    if (!f)
      return false;
    else
      myScalars[code] = f;
  }
  else if (direction < 0 || this->getNoFields() == 1)
  {
    VecFunc* f = utl::parseVecFunc(prop,type);
    if (!f)
      return false;
    else
      myVectors[code] = f;
  }
  else
  {
    TractionFunc* f = utl::parseTracFunc(prop,type,direction);
    if (!f)
      return false;
    else
      myTracs[code] = f;
  }

  return this->setPropertyType(code,Property::NEUMANN);
}


VecFunc* SIMbase::getVecFunc (size_t patch, Property::Type ptype) const
{
  for (PropertyVec::const_iterator p = myProps.begin(); p != myProps.end(); ++p)
    if (p->patch == patch)
      if (p->pcode == ptype || ptype == Property::UNDEFINED)
      {
        VecFuncMap::const_iterator it = myVectors.find(p->pindx);
        if (it != myVectors.end()) return it->second;
      }

  return nullptr;
}


RealFunc* SIMbase::getSclFunc (int code) const
{
  SclFuncMap::const_iterator it = myScalars.find(code);
  return it == myScalars.end() ? nullptr : it->second;
}


bool SIMbase::initSystem (int mType, size_t nMats, size_t nVec, bool withRF)
{
  if (!mySam) return true; // Silently ignore when no algebraic system

#if SP_DEBUG > 2
  mySam->print(std::cout);
  std::string heading("\n\nNodal coordinates for Patch 1");
  size_t n, i, j = heading.size()-1;
  typedef std::pair<int,int> Ipair;
  typedef std::vector<Ipair> Ipairs;
  std::map<int,Ipairs> nodeInfo;
  std::map<int,Ipairs>::const_iterator it;
  for (i = 0; i < myModel.size() && i < 9; i++, heading[j]++)
    if (!myModel[i]->empty())
    {
      myModel[i]->printNodes(std::cout,heading.c_str());
      for (n = 1; n <= myModel[i]->getNoNodes(); n++)
        nodeInfo[myModel[i]->getNodeID(n)].push_back(std::make_pair(i+1,n));
    }

  for (it = nodeInfo.begin(); it != nodeInfo.end(); ++it)
    if (it->second.size() > 1)
    {
      std::cout <<"\nConnectivity for node "<< it->first <<":";
      for (n = 0; n < it->second.size(); n++)
        std::cout <<" P"<< it->second[n].first <<","<< it->second[n].second;
    }
  std::cout << std::endl;
#endif

  if (myEqSys) delete myEqSys;
  myEqSys = new AlgEqSystem(*mySam,adm);

  // Workaround SuperLU bug for tiny systems
  if (mType == SystemMatrix::SPARSE && this->getNoElms(true) < 3)
  {
    std::cerr <<" ** System too small for SuperLU, falling back to Dense."
              << std::endl;
    mType = SystemMatrix::DENSE;
  }

  return myEqSys->init(static_cast<SystemMatrix::Type>(mType),
                       mySolParams, nMats, nVec, withRF,
                       myProblem->getLinearSystemType(), opt.num_threads_SLU);
}


bool SIMbase::setAssociatedRHS (size_t iMat, size_t iVec)
{
  if (!myEqSys) return false;

  return myEqSys->setAssociatedVector(iMat,iVec);
}


bool SIMbase::setMode (int mode, bool resetSol)
{
  if (myInts.empty())
    myInts.insert(std::make_pair(0,myProblem));

  for (IntegrandMap::iterator it = myInts.begin(); it != myInts.end(); ++it)
    if (it->second)
    {
      it->second->setMode((SIM::SolutionMode)mode);
      if (resetSol) it->second->resetSolution();
    }
    else
    {
      std::cerr <<" *** SIMbase::setMode: No integrand yet";
      if (it->first > 0) std::cerr <<", code="<< it->first;
      std::cerr << std::endl;
      return false;
    }

  return true;
}


void SIMbase::setIntegrationPrm (unsigned short int i, double prm)
{
  if (!myInts.empty())
    for (IntegrandMap::iterator it = myInts.begin(); it != myInts.end(); ++it)
      it->second->setIntegrationPrm(i,prm);
  else if (myProblem)
    myProblem->setIntegrationPrm(i,prm);
  else
    std::cerr <<"  ** SIMbase::setIntegrationPrm: myProblem not set yet."
              << std::endl;
}


void SIMbase::setQuadratureRule (size_t ng, bool redimBuffers, bool printQP)
{
  size_t nInterfaceGP = nIntGP = nBouGP = 0;
  for (size_t i = 0; i < myModel.size(); i++)
    if (!myModel[i]->empty())
    {
      myModel[i]->setGauss(ng);
      // Count the interior integration points
      myModel[i]->getNoIntPoints(nIntGP,nInterfaceGP);
    }

  if (!myProblem) return;

  for (PropertyVec::const_iterator p = myProps.begin(); p != myProps.end(); ++p)
    if (p->pcode == Property::NEUMANN && myProblem->hasBoundaryTerms())
    {
      // Account for possibly more than one Neumann property on a boundary
      bool notCounted = true;
      for (PropertyVec::const_iterator q = myProps.begin(); q != p; ++q)
	if (q->patch == p->patch && q->lindx == p->lindx && q->pcode==p->pcode)
	  notCounted = false;

      if (notCounted) // Count the boundary integration points
        myModel[p->patch-1]->getNoBouPoints(nBouGP,abs(p->ldim),p->lindx);
    }

  // Let the integrands know how many integration points in total we have
  if (redimBuffers)
    for (IntegrandMap::iterator it = myInts.begin(); it != myInts.end(); ++it)
      it->second->initIntegration(nIntGP,nBouGP);

  if (printQP)
  {
    IFEM::cout <<"Number of quadrature points "<< nIntGP;
    if (nInterfaceGP > 0) IFEM::cout <<" "<< nInterfaceGP;
    if (nBouGP > 0) IFEM::cout <<" "<< nBouGP;
    IFEM::cout << std::endl;
  }
}


void SIMbase::printProblem () const
{
  if (myProblem)
  {
    IFEM::cout <<"\nProblem definition:"<< std::endl;
    myProblem->printLog();
  }

#if SP_DEBUG > 1
  std::cout <<"\nProperty mapping:";
  for (PropertyVec::const_iterator p = myProps.begin(); p != myProps.end(); ++p)
    std::cout <<"\n"<< p->pcode <<" "<< p->pindx <<" "<< p->patch
              <<" "<< (int)p->lindx <<" "<< (int)p->ldim;
  std::cout << std::endl;
#endif
}


size_t SIMbase::getNoFields (int basis) const
{
  return myModel.empty() ? 0 : myModel.front()->getNoFields(basis);
}


size_t SIMbase::getNoDOFs () const
{
  return mySam ? mySam->getNoDOFs() : 0;
}


size_t SIMbase::getNoNodes (bool unique, int basis) const
{
  size_t nnod = 0;
  if (unique && mySam)
    nnod = mySam->getNoNodes(basis < 1 ? 'A' : (basis < 2 ? 'D' : 'N'+basis));
  else
    for (size_t i = 0; i < myModel.size(); i++)
      nnod += myModel[i]->getNoNodes(basis);

  return nnod;
}


size_t SIMbase::getNoElms (bool includeXelms) const
{
  if (mySam && includeXelms)
    return mySam->getNoElms();

  size_t noElms = 0;
  for (size_t i = 0; i < myModel.size(); i++)
    noElms += myModel[i]->getNoElms(false,includeXelms);

  return noElms;
}


size_t SIMbase::getNoSolutions () const
{
  return myProblem ? myProblem->getNoSolutions() : 0;
}


size_t SIMbase::getNoRHS () const
{
  return myEqSys ? myEqSys->getNoRHS() : 1;
}


char SIMbase::getNoBasis () const
{
  size_t result = myModel.empty() ? 0 : myModel.front()->getNoBasis();
#ifdef SP_DEBUG
  for (size_t i = 1; i < myModel.size(); i++)
    assert(myModel[i]->getNoBasis() == result);
#endif

  return result;
}


bool SIMbase::hasTimeDependentDirichlet () const
{
  for (size_t i = 0; i < myModel.size(); i++)
    if (myModel[i]->hasTimeDependentDirichlet(myScalars,myVectors))
      return true;

  return false;
}


bool SIMbase::initDirichlet (double time)
{
  if (time == 0.0)
    for (size_t i = 0; i < myModel.size(); i++)
      if (!myModel[i]->initConstraints())
	return false;

  Vector dummy;
  return this->updateDirichlet(time,&dummy);
}


bool SIMbase::updateDirichlet (double time, const Vector* prevSol)
{
  if (prevSol)
    for (size_t i = 0; i < myModel.size(); i++)
#ifdef HAVE_MPI
      if (!myModel[i]->updateDirichlet(myScalars,myVectors,time,g2l))
#else
      if (!myModel[i]->updateDirichlet(myScalars,myVectors,time))
#endif
	return false;

  SAMpatch* pSam = dynamic_cast<SAMpatch*>(mySam);
  return pSam ? pSam->updateConstraintEqs(myModel,prevSol) : true;
}


bool SIMbase::updateGrid (const Vector& displ)
{
  if (displ.empty()) return true; // No displacements (yet), totally fine

  bool ok = true;
  Vector locdisp;
  for (size_t i = 0; i < myModel.size() && ok; i++)
  {
    myModel[i]->extractNodeVec(displ,locdisp,myModel[i]->getNoSpaceDim(),-1);
    ok = myModel[i]->updateCoords(locdisp);
  }

  return ok;
}


bool SIMbase::updateGrid (const std::string& field)
{
  const Vector* displ = this->getDependentField(field);
  if (displ) return this->updateGrid(*displ);

  std::cerr <<" *** SIMbase::updateGrid: No such field \""<< field
	    <<"\" registered for \""<< this->getName() <<"\"."<< std::endl;
  return false;
}


void SIMbase::getBoundaryNodes (int pcode, IntVec& glbNodes, Vec3Vec* XYZ) const
{
  glbNodes.clear();
  if (XYZ) XYZ->clear();

  ASMbase* pch;
  size_t node;
  for (PropertyVec::const_iterator p = myProps.begin(); p != myProps.end(); ++p)
    if (abs(p->pindx) == pcode && (pch = this->getPatch(p->patch))) {
      if (abs(p->ldim)+1 == pch->getNoParamDim()) {
        // The boundary is of one dimension lower than the patch
        IntVec nodes;
        pch->getBoundaryNodes(abs(p->lindx),nodes);
        for (const auto& it : nodes)
          if (std::find(glbNodes.begin(),glbNodes.end(),it) == glbNodes.end()) {
            glbNodes.push_back(it);
            if (XYZ) {
              if ((node = pch->getNodeIndex(it,true)))
                XYZ->push_back(pch->getCoord(node));
              else
                XYZ->push_back(Vec3());
            }
          }
      }
      else if (pch->getNoParamDim() == abs(p->ldim)) {
        // The boundary and the patch are of same dimension
        for (node = 1; node <= pch->getNoNodes(); node++)
        {
          glbNodes.push_back(pch->getNodeID(node));
          if (XYZ) XYZ->push_back(pch->getCoord(node));
        }
      }
    }
}


int SIMbase::findClosestNode (const Vec3& X) const
{
  if (myModel.empty()) return -1;

  ASMbase* closestPch = nullptr;
  std::pair<size_t,double> closest(0,1.0e99);
  for (size_t i = 0; i < myModel.size(); i++)
  {
    std::pair<size_t,double> node = myModel[i]->findClosestNode(X);
    if (node.first > 0 && node.second < closest.second)
    {
      closest = node;
      closestPch = myModel[i];
    }
  }
  if (!closestPch) return -2;

#ifdef SP_DEBUG
  std::cout <<"SIMbase::findClosestNode("<< X <<") -> Node "<< closest.first
            <<" in Patch "<< closestPch->idx+1 <<" distance="<< closest.second
            << std::endl;
#endif

  return closestPch->getNodeID(closest.first);
}


bool SIMbase::assembleSystem (const TimeDomain& time, const Vectors& prevSol,
			      bool newLHSmatrix, bool poorConvg)
{
  PROFILE1("Element assembly");

  bool ok = true;
  bool isAssembling = (myProblem->getMode() != SIM::INIT &&
                       myProblem->getMode() != SIM::RECOVERY);
  if (isAssembling)
    myEqSys->initialize(newLHSmatrix);

  // Loop over the integrands
  IntegrandMap::const_iterator it;
  for (it = myInts.begin(); it != myInts.end() && ok; ++it)
  {
    if (msgLevel > 1)
      IFEM::cout <<"\n\nProcessing integrand associated with code "<< it->first
                << std::endl;

    GlobalIntegral& sysQ = it->second->getGlobalInt(myEqSys);
    if (&sysQ != myEqSys && isAssembling)
      sysQ.initialize(newLHSmatrix);

    if (!prevSol.empty())
      it->second->initIntegration(time,prevSol.front(),poorConvg);

    // Loop over the different material regions, integrating interior
    // coefficient matrix terms for the patch associated with each material
    size_t lp = 0;
    ASMbase* pch = nullptr;
    PropertyVec::const_iterator p, p2;
    if (it->second->hasInteriorTerms())
    {
      for (p = myProps.begin(); p != myProps.end() && ok; ++p)
        if (p->pcode == Property::MATERIAL &&
            (it->first == 0 || it->first == p->pindx))
          if (!(pch = this->getPatch(p->patch)))
          {
            std::cerr <<" *** SIMbase::assembleSystem: Patch index "<< p->patch
                      <<" out of range [1,"<< myModel.size() <<"]"<< std::endl;
            ok = false;
          }
          else if (this->initMaterial(p->pindx))
          {
            lp = p->patch;
            if (msgLevel > 1)
              IFEM::cout <<"\nAssembling interior matrix terms for P"<< lp
                         << std::endl;
            ok &= this->initBodyLoad(lp);
            ok &= this->extractPatchSolution(it->second,prevSol,lp-1);
            ok &= pch->integrate(*it->second,sysQ,time);
          }
          else
            ok = false;

      if (lp == 0 && it->first == 0)
        // All patches refer to the same material, and we assume it has been
        // initialized during input processing (thus no initMaterial call here)
        for (size_t k = 0; k < myModel.size() && ok; k++)
        {
          lp = k+1;
          if (msgLevel > 1)
            IFEM::cout <<"\nAssembling interior matrix terms for P"<< lp
                       << std::endl;
          ok &= this->initBodyLoad(lp);
          ok &= this->extractPatchSolution(it->second,prevSol,k);
          ok &= myModel[k]->integrate(*it->second,sysQ,time);
        }
    }

    // Assemble contributions from the Neumann boundary conditions
    // and other boundary integrals (Robin properties, contact, etc.)
    if (it->second->hasBoundaryTerms() && myEqSys->getVector())
      for (p = myProps.begin(); p != myProps.end() && ok; ++p)
        if ((p->pcode == Property::NEUMANN && it->first == 0) ||
            ((p->pcode == Property::NEUMANN_GENERIC ||
              p->pcode == Property::ROBIN) && it->first == p->pindx))
        {
          if (!(pch = this->getPatch(p->patch)))
          {
            std::cerr <<" *** SIMbase::assembleSystem: Patch index "<< p->patch
                      <<" out of range [1,"<< myModel.size() <<"]"<< std::endl;
            ok = false;
            break;
          }

          for (p2 = myProps.begin(); p2 != myProps.end() && ok; ++p2)
            if (p2->pcode == Property::MATERIAL && p->patch == p2->patch)
              if (!(ok = this->initMaterial(p2->pindx)))
                std::cerr <<" *** SIMbase::assembleSystem: Failed to initialize"
                          <<" material for patch "<< p2->patch << std::endl;

          if (abs(p->ldim)+1 == pch->getNoParamDim())
          {
            if (p->pcode == Property::NEUMANN_GENERIC ||
                this->initNeumann(p->pindx))
            {
              if (msgLevel > 1)
                IFEM::cout <<"\nAssembling Neumann matrix terms for boundary "
                           << (int)p->lindx <<" on P"<< p->patch << std::endl;
              if (p->patch != lp)
                ok &= this->extractPatchSolution(it->second,prevSol,p->patch-1);
              ok &= pch->integrate(*it->second,p->lindx,sysQ,time);
              lp = p->patch;
            }
            else
              ok = false;
          }
          else if (abs(p->ldim) == 1 && pch->getNoParamDim() == 3)
          {
            if (p->pcode == Property::NEUMANN_GENERIC ||
                this->initNeumann(p->pindx))
            {
              if (msgLevel > 1)
                IFEM::cout <<"\nAssembling Neumann matrix terms for edge "
                           << (int)p->lindx <<" on P"<< p->patch << std::endl;
              if (p->patch != lp)
                ok &= this->extractPatchSolution(it->second,prevSol,p->patch-1);
              ok &= pch->integrateEdge(*it->second,p->lindx,sysQ,time);
              lp = p->patch;
            }
            else
              ok = false;
          }
        }

    if (ok) ok = this->assembleDiscreteTerms(it->second,time);
    if (ok && &sysQ != myEqSys && isAssembling)
      ok = sysQ.finalize(newLHSmatrix);
  }
  if (ok && isAssembling)
    ok = myEqSys->finalize(newLHSmatrix);

  if (!ok)
    std::cerr <<" *** SIMbase::assembleSystem: Failure.\n"<< std::endl;

  return ok;
}


bool SIMbase::extractLoadVec (Vector& loadVec) const
{
  // Expand load vector from equation ordering to DOF-ordering
  SystemVector* b = myEqSys->getVector();
  if (!b || !mySam->expandSolution(*b,loadVec,0.0))
    return false;

#if SP_DEBUG > 1
  std::cout <<"\nLoad vector:"<< loadVec;
#endif
  return true;
}


bool SIMbase::applyDirichlet (Vector& glbVec) const
{
  return mySam->applyDirichlet(glbVec);
}


bool SIMbase::solveSystem (Vector& solution, int printSol,
                           const char* compName, bool newLHS, size_t idxRHS)
{
  SystemMatrix* A = myEqSys->getMatrix();
  SystemVector* b = myEqSys->getVector(idxRHS);
  if (!A) std::cerr <<" *** SIMbase::solveSystem: No LHS matrix"<< std::endl;
  if (!b) std::cerr <<" *** SIMbase::solveSystem: No RHS vector"<< std::endl;
  if (!A || !b) return false;

  // Dump system matrix to file, if requested
  std::vector<DumpData>::iterator it;
  for (it = lhsDump.begin(); it != lhsDump.end(); ++it)
    if (it->doDump()) {
      IFEM::cout <<"\nDumping system matrix to file "<< it->fname << std::endl;
      std::ofstream os(it->fname.c_str());
      os << std::setprecision(17);
      SystemMatrix* M = myEqSys->getMatrix(0);
      char matName[] = {'A'};
      for (int i = 0; M; M = myEqSys->getMatrix(++i), ++matName[0])
        M->dump(os,it->format,matName); // label matrices as A,B,C,...
    }

  // Dump right-hand-side vector to file, if requested
  for (it = rhsDump.begin(); it != rhsDump.end(); ++it)
    if (it->doDump()) {
      IFEM::cout <<"\nDumping RHS vector to file "<< it->fname << std::endl;
      std::ofstream os(it->fname.c_str());
      os << std::setprecision(17);
      SystemVector* c = myEqSys->getVector(0);
      char vecName[] = {'b'};
      for (int i = 0; c; c = myEqSys->getVector(++i), ++vecName[0])
        c->dump(os,it->format,vecName); // label vectors as b,c,d,...
    }

  // Solve the linear system of equations
  bool status = true;
  double rCond = 0.0;
  if (msgLevel > 1)
  {
    IFEM::cout <<"\nSolving the equation system ..."<< std::endl;
    PROFILE1("Equation solving");
    status = A->solve(*b,newLHS,&rCond);
  }
  else
  {
    PROFILE1("Equation solving");
    status = A->solve(*b,newLHS);
  }
  if (rCond > 0.0)
    IFEM::cout <<"\tCondition number: "<< 1.0/rCond << std::endl;

  // Dump solution vector to file, if requested
  for (it = solDump.begin(); it != solDump.end() && status; ++it)
    if (it->doDump()) {
      IFEM::cout <<"\nDumping solution vector to file "<< it->fname << std::endl;
      std::ofstream os(it->fname.c_str());
      os << std::setprecision(17);
      b->dump(os,it->format,"b");
    }

  // Expand solution vector from equation ordering to DOF-ordering
  if (status)
    status = mySam->expandSolution(*b,solution);

  if (printSol > 0 && status)
    this->printSolutionSummary(solution,printSol,compName);

  return status;
}


bool SIMbase::solveMatrixSystem (Vectors& solution, int printSol,
                                 const char* compName)
{
  solution.resize(myEqSys->getNoRHS());
  for (size_t i = 0; i < solution.size(); i++)
    if (!this->solveSystem(solution[i],printSol,compName,i==0,i))
      return false;
    else
      printSol = 0; // Print summary only for the first solution

  return true;
}


void SIMbase::printSolutionSummary (const Vector& solution, int printSol,
                                    const char* compName,
                                    std::streamsize outPrec)
{
  // Compute and print solution norms
  const size_t nf = this->getNoFields(1);
  size_t iMax[nf];
  double dMax[nf];
  double dNorm = this->solutionNorms(solution,dMax,iMax,nf);

  int oldPrec = adm.cout.precision();
  if (outPrec > 0)
    adm.cout << std::setprecision(outPrec);

  if (compName)
    adm.cout <<"\n >>> Solution summary <<<\n\nL2-norm            : ";
  else
    adm.cout <<"  Primary solution summary: L2-norm         : ";
  adm.cout << utl::trunc(dNorm);

  if (nf == 1 && utl::trunc(dMax[0]) != 0.0)
  {
    if (compName)
      adm.cout <<"\nMax "<< compName <<"   : ";
    else
      adm.cout <<"\n                            Max value       : ";
    adm.cout << dMax[0] <<" node "<< iMax[0];
  }
  else if (nf > 1)
  {
    char D = 'X';
    for (size_t d = 0; d < nf; d++, D=='Z' ? D='x' : D++)
      if (utl::trunc(dMax[d]) != 0.0)
      {
        if (compName)
          adm.cout <<"\nMax "<< D <<'-'<< compName <<" : ";
        else
          adm.cout <<"\n                            Max "<< D <<"-component : ";
        adm.cout << dMax[d] <<" node "<< iMax[d];
      }
  }
  adm.cout << std::endl;
  adm.cout << std::setprecision(oldPrec);

  // Print entire solution vector if it is small enough
  if (mySam->getNoEquations() < printSol)
  {
    adm.cout <<"\nSolution vector:";
    for (int inod = 1; inod <= mySam->getNoNodes(); inod++)
    {
      adm.cout <<"\nNode "<< inod <<":";
      std::pair<int,int> dofs = mySam->getNodeDOFs(inod);
      for (int d = dofs.first-1; d < dofs.second; d++)
        adm.cout <<" "<< utl::trunc(solution[d]);
    }
    adm.cout << std::endl;
  }
#if SP_DEBUG > 2
  else
    std::cout <<"\nSolution vector:"<< *myEqSys->getVector();
#endif
}


void SIMbase::printNorms (const Vectors& norms, size_t w) const
{
  if (norms.empty()) return;

  NormBase* norm = this->getNormIntegrand();
  const Vector& n = norms.front();

  IFEM::cout <<"Energy norm"
             << utl::adjustRight(w-11,norm->getName(1,1)) << n(1)
             <<"\nExternal energy"
             << utl::adjustRight(w-15,norm->getName(1,2)) << n(2);

  if (mySol)
    IFEM::cout <<"\nExact norm"
               << utl::adjustRight(w-10,norm->getName(1,3)) << n(3)
               <<"\nExact error"
               << utl::adjustRight(w-11,norm->getName(1,4)) << n(4)
               <<"\nExact relative error (%) : "
               << 100.0*n(4)/n(3);

  IFEM::cout << std::endl;
  delete norm;
}


void SIMbase::getWorstDofs (const Vector& u, const Vector& r,
                            size_t nWorst, double eps,
                            std::map<std::pair<int,int>,RealArray>& worst) const
{
  size_t i;
  RealArray data(3);
  std::multimap<double,size_t> energy;

  // Compute the energy at each DOF and insert into a map sorted on the energy
  for (i = 0; i < u.size() && i < r.size(); i++)
    energy.insert(std::make_pair(fabs(u[i]*r[i]),i+1));

  // Pick the nWorst highest energies from the back of the map
  std::multimap<double,size_t>::reverse_iterator rit = energy.rbegin();
  for (i = 0; i < nWorst && rit != energy.rend(); i++, ++rit)
    if (rit->first > eps)
    {
      data[0] = rit->first;
      data[1] = u(rit->second);
      data[2] = r(rit->second);
      worst[mySam->getNodeAndLocalDof(rit->second)] = data;
    }
}


char SIMbase::getNodeType (int inod) const
{
  return mySam ? mySam->getNodeType(inod) : ' ';
}


Vec4 SIMbase::getNodeCoord (int inod) const
{
  Vec4 Xnod;
  size_t node = 0;
  for (PatchVec::const_iterator it = myModel.begin(); it != myModel.end(); ++it)
    if ((node = (*it)->getNodeIndex(inod,true)))
    {
      Xnod = (*it)->getCoord(node);
      if (myModel.size() > 1)
        Xnod.idx = (*it)->idx; // Store patch index, if multi-patch model
      break;
    }

  return Xnod;
}


bool SIMbase::isFixed (int inod, int dof) const
{
#ifdef HAS_PETSC
  SAMpatchPara* pSam = dynamic_cast<SAMpatchPara*>(mySam);
  if (pSam) return pSam->isDirichlet(inod,dof);
#endif

  size_t node = 0;
  for (PatchVec::const_iterator it = myModel.begin(); it != myModel.end(); ++it)
    if ((node = (*it)->getNodeIndex(inod,true)))
      return (*it)->isFixed(node,dof,true);

  return true;
}


int SIMbase::getGlobalNode (int node) const
{
  std::map<int,int>::const_iterator it = utl::findValue(*g2l,node);
  return it != g2l->end() ? it->first : -1;
}


int SIMbase::getLocalNode (int node) const
{
  std::map<int,int>::const_iterator it = g2l->find(node);
  return it != g2l->end() ? it->second : -1;
}


SystemVector* SIMbase::getRHSvector (size_t idx, bool copy) const
{
  SystemVector* rhs = myEqSys->getVector(idx);
  return rhs && copy ? rhs->copy() : rhs;
}


void SIMbase::addToRHSvector (size_t idx, const SystemVector& vec, double scale)
{
  SystemVector* rhs = myEqSys->getVector(idx);
  if (!rhs || scale == 0.0) return;

  rhs->add(vec,scale);
}


void SIMbase::iterationNorms (const Vector& u, const Vector& r,
			      double& eNorm, double& rNorm, double& dNorm) const
{
  eNorm = mySam->dot(r,u,'A');
  rNorm = mySam->norm2(r,'D');
  dNorm = mySam->norm2(u,'D');
}


double SIMbase::solutionNorms (const Vector& x, double* inf,
			       size_t* ind, size_t nf, char type) const
{
  if (inf && ind && nf == 0) nf = nsd;

  for (size_t d = 0; d < nf; d++)
  {
    ind[d] = d+1;
    inf[d] = mySam->normInf(x,ind[d],type);
  }

  return mySam->normL2(x,type);
}


NormBase* SIMbase::getNormIntegrand () const
{
  return myProblem->getNormIntegrand(mySol);
}


ForceBase* SIMbase::getBoundaryForceIntegrand (const Vec3* X0) const
{
  return myProblem->getForceIntegrand(X0,mySol);
}


ForceBase* SIMbase::getNodalForceIntegrand () const
{
  return myProblem->getForceIntegrand();
}


bool SIMbase::solutionNorms (const TimeDomain& time,
			     const Vectors& psol, const Vectors& ssol,
			     Vectors& gNorm, Matrix* eNorm)
{
  PROFILE1("Norm integration");

  if (!mySam) return true; // Silently ignore when uninitialized system

  NormBase* norm = myProblem->getNormIntegrand(mySol);
  if (!norm)
  {
#ifdef SP_DEBUG
    std::cerr <<"  ** SIMbase::solutionNorms: No integrand."<< std::endl;
#endif
    return true; // Silently ignore when no norm integrand is provided
  }

  myProblem->initIntegration(time,psol.front());
  norm->initProjection(ssol.size());
  norm->initIntegration(nIntGP,nBouGP);

  // Number of recovered solution components
  size_t nCmp = ssol.empty() ? 0 : ssol.front().size() / mySam->getNoNodes();

#ifdef USE_OPENMP
  // When assembling in parallel, we must always do the norm summation
  // at the end in a serial loop, to avoid that the threads try to update
  // the same memory address simultaneously.
  Matrix dummy;
  if (!eNorm) eNorm = &dummy;
#endif

  // Initialize norm integral classes
  gNorm.resize(norm->getNoFields(0));
  size_t nNorms = 0;
  for (size_t j = 0; j < gNorm.size(); j++) {
    size_t nNrm = norm->getNoFields(1+j);
    gNorm[j].resize(nNrm,true);
    nNorms += nNrm;
  }

  GlbNorm globalNorm(gNorm,norm->getFinalOperation());
  LintegralVec elementNorms;
  if (eNorm)
  {
    eNorm->resize(nNorms,mySam->getNoElms(),true);
    elementNorms.reserve(eNorm->cols());
    for (size_t i = 0; i < eNorm->cols(); i++)
      elementNorms.push_back(new ElmNorm(eNorm->ptr(i),nNorms));
    norm->setLocalIntegrals(&elementNorms);
  }

  // Loop over the different material regions, integrating solution norm terms
  // for the patch domain associated with each material
  bool ok = true;
  size_t k, lp = 0;
  ASMbase* pch = nullptr;
  PropertyVec::const_iterator p;
  for (p = myProps.begin(); p != myProps.end() && ok; ++p)
    if (p->pcode == Property::MATERIAL)
      if (!(pch = this->getPatch(p->patch)))
	ok = false;
      else if (this->initMaterial(p->pindx))
      {
	lp = p->patch;
	ok = this->extractPatchSolution(psol,lp-1);
	for (k = 0; k < ssol.size(); k++)
	  if (!ssol[k].empty())
	    pch->extractNodeVec(ssol[k],norm->getProjection(k),nCmp);
	ok &= pch->integrate(*norm,globalNorm,time);
      }
      else
	ok = false;

  if (lp == 0)
    // All patches are referring to the same material, and we assume it has
    // been initialized during input processing (thus no initMaterial call here)
    for (size_t i = 0; i < myModel.size() && ok; i++)
    {
      ok = this->extractPatchSolution(psol,i);
      for (k = 0; k < ssol.size(); k++)
	if (!ssol[k].empty())
	  myModel[i]->extractNodeVec(ssol[k],norm->getProjection(k),nCmp);
      ok &= myModel[i]->integrate(*norm,globalNorm,time);
      lp = i+1;
    }

  if (norm->hasBoundaryTerms())
    for (p = myProps.begin(); p != myProps.end() && ok; ++p)
      if (p->pcode == Property::NEUMANN)
	if (!(pch = this->getPatch(p->patch)))
	  ok = false;

	else if (abs(p->ldim)+1 == pch->getNoParamDim())
	  if (this->initNeumann(p->pindx))
	  {
	    if (p->patch != lp)
	      ok = this->extractPatchSolution(psol,p->patch-1);
	    ok &= pch->integrate(*norm,p->lindx,globalNorm,time);
	    lp = p->patch;
	  }
	  else
	    ok = false;

	else if (abs(p->ldim)+2 == pch->getNoParamDim())
	  if (this->initNeumann(p->pindx))
	  {
	    if (p->patch != lp)
	      ok = this->extractPatchSolution(psol,p->patch-1);
	    ok &= pch->integrateEdge(*norm,p->lindx,globalNorm,time);
	    lp = p->patch;
	  }
	  else
	    ok = false;

  if (!ok) std::cerr <<" *** SIMbase::solutionNorms: Failure.\n"<< std::endl;

  // Clean up the dynamically allocated norm objects. This will also perform
  // the actual global norm assembly, in case the element norms are stored,
  // and always when multi-threading is used.
  for (k = 0; k < elementNorms.size(); k++)
  {
    globalNorm.assemble(elementNorms[k]);
    delete elementNorms[k];
  }

  // Add problem-dependent external norm contributions
  norm->addBoundaryTerms(gNorm,this->externalEnergy(psol));

  delete norm;

  for (k = 0; k < gNorm.size(); k++)
    adm.allReduceAsSum(gNorm[k]);

  return ok;
}


double SIMbase::externalEnergy (const Vectors& psol) const
{
  if (!myEqSys || !mySam)
    return 0.0;

  const Vector* reactionForces = myEqSys->getReactions();
  if (!reactionForces || psol.empty()) return 0.0;

  // Add norm contributions due to inhomogeneous Dirichlet boundary conditions.
  // That is, the path integral of the total solution vector times the
  // reaction forces at the prescribed DOFs.
  if (psol.size() == 1)
    return mySam->normReact(psol.front(),*reactionForces);

  static double extEnergy = 0.0;
  static Vector prevForces(reactionForces->size());
  extEnergy += mySam->normReact(psol[0]-psol[1],*reactionForces+prevForces);
  prevForces = *reactionForces;
  return extEnergy;
}


/*!
  The content of the output array \a RF is as follows:
  \f[
  RF[0] = \sum_{i=n}^{\rm nnod} \sum_{i=1}^{\rm nsd} f_n^i u_n^i
  \quad\quad\mbox{(energy)}\f]
  \f[
  RF[i] = \sum_{n=1}^{\rm nnod} f_n^i \quad\forall\quad i=1,\ldots,{\rm nsd}
  \f]
*/

bool SIMbase::getCurrentReactions (RealArray& RF, const Vector& psol) const
{
  const Vector* reactionForces = myEqSys->getReactions();
  if (!reactionForces || !mySam) return false;

  RF.resize(1+nsd);
  RF.front() = 2.0*mySam->normReact(psol,*reactionForces);
  for (size_t dir = 1; dir < RF.size(); dir++)
    RF[dir] = mySam->getReaction(dir,*reactionForces);

  return true;
}


bool SIMbase::getCurrentReactions (RealArray& RF, int pcode) const
{
  const Vector* reactionForces = myEqSys->getReactions();
  if (!reactionForces || !mySam) return false;

  IntVec glbNodes;
  this->getBoundaryNodes(pcode,glbNodes);

  RF.resize(nsd);
  for (unsigned char dir = 1; dir <= nsd; dir++)
    RF[dir] = mySam->getReaction(dir,*reactionForces,&glbNodes);

  return true;
}


bool SIMbase::systemModes (std::vector<Mode>& solution,
			   int nev, int ncv, int iop, double shift,
			   size_t iA, size_t iB)
{
  if (nev < 1 || ncv <= nev) return false;

  PROFILE1("Eigenvalue analysis");

  Vector eigVal;
  Matrix eigVec;
  if (nev > mySam->getNoEquations()) nev = mySam->getNoEquations();
  if (ncv > mySam->getNoEquations()) ncv = mySam->getNoEquations();

  // Solve the eigenvalue problem
  IFEM::cout <<"\nSolving the eigenvalue problem ..."<< std::endl;
  SystemMatrix* A = myEqSys->getMatrix(iA);
  SystemMatrix* B = myEqSys->getMatrix(iB);
#ifdef HAS_SLEPC
  // To interface SLEPC another interface is used
  bool ok = eig::solve(A,B,eigVal,eigVec,nev);
#else
  bool ok = eig::solve(A,B,eigVal,eigVec,nev,ncv,iop,shift);
#endif

  // Expand eigenvectors to DOF-ordering and print out eigenvalues
  bool freq = iop == 3 || iop == 4 || iop == 6;
  IFEM::cout <<"\n >>> Computed Eigenvalues <<<\n     Mode\t"
             << (freq ? "Frequency [Hz]" : "Eigenvalue");
  solution.resize(nev);
  for (int i = 1; i <= nev && ok; i++)
  {
    solution[i-1].eigNo = i;
    if (!mySam->expandVector(eigVec.getColumn(i),solution[i-1].eigVec))
      ok = false;
    else if (!freq)
      solution[i-1].eigVal = eigVal(i);
    else if (eigVal(i) < 0.0)
      solution[i-1].eigVal = -sqrt(-eigVal(i))*0.5/M_PI;
    else
      solution[i-1].eigVal =  sqrt( eigVal(i))*0.5/M_PI;

    IFEM::cout <<"\n     "<< i <<"\t\t"<< utl::trunc(solution[i-1].eigVal);
  }
  IFEM::cout << std::endl;
  return ok;
}


bool SIMbase::project (Matrix& ssol, const Vector& psol,
		       SIMoptions::ProjectionMethod pMethod,
		       const TimeDomain& time) const
{
  PROFILE1("Solution projection");

  if (msgLevel > 1)
    IFEM::cout <<"\nProjecting secondary solution ..."<< std::endl;

  ssol.clear();

  size_t i, j, n;
  size_t ngNodes = mySam->getNoNodes();

  Matrix values;
  Vector count(myModel.size() > 1 ? ngNodes : 0);

  if (pMethod == SIMoptions::DGL2 || pMethod == SIMoptions::CGL2)
  {
#ifdef USE_OPENMP
  if (!myModel.empty() && dynamic_cast<ASMunstruct*>(myModel.front()))
    omp_set_num_threads(1);
#endif

    // Reinitialize the integration point buffers within the integrands (if any)
    const_cast<SIMbase*>(this)->setQuadratureRule(opt.nGauss[0]);
    myProblem->initIntegration(time,psol);
  }

  for (i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    // Extract the primary solution control point values for this patch
    myProblem->initResultPoints(time.t);
    if (!this->extractPatchSolution(myProblem,Vectors(1,psol),i))
      return false;

    // Initialize material properties for this patch in case of multiple regions
    const_cast<SIMbase*>(this)->setPatchMaterial(i+1);

    // Project the secondary solution and retrieve control point values
    bool ok = false;
    switch (pMethod) {
    case SIMoptions::GLOBAL:
      if (msgLevel > 1 && i == 0)
        IFEM::cout <<"\tGreville point projection"<< std::endl;
      ok = myModel[i]->evalSolution(values,*myProblem);
      break;

    case SIMoptions::DGL2:
      if (msgLevel > 1 && i == 0)
        IFEM::cout <<"\tDiscrete global L2-projection"<< std::endl;
      ok = myModel[i]->globalL2projection(values,*myProblem);
      break;

    case SIMoptions::CGL2:
      if (msgLevel > 1 && i == 0)
        IFEM::cout <<"\tContinuous global L2-projection"<< std::endl;
      ok = myModel[i]->L2projection(values,*myProblem,time);
      break;

    case SIMoptions::SCR:
      if (msgLevel > 1 && i == 0)
        IFEM::cout <<"\tSuperconvergent recovery"<< std::endl;
      ok = myModel[i]->evalSolution(values,*myProblem,nullptr,'S');
      break;

    case SIMoptions::VDSA:
      if (msgLevel > 1 && i == 0)
        IFEM::cout <<"\tVariation diminishing projection"<< std::endl;
      ok = myModel[i]->evalSolution(values,*myProblem,nullptr,'A');
      break;

    case SIMoptions::QUASI:
      if (msgLevel > 1 && i == 0)
        IFEM::cout <<"\tQuasi interpolation"<< std::endl;
      ok = myModel[i]->evalSolution(values,*myProblem,nullptr,'L');
      break;

    case SIMoptions::LEASTSQ:
      if (msgLevel > 1 && i == 0)
	IFEM::cout <<"\tLeast squares projection"<< std::endl;
      ok = myModel[i]->evalSolution(values,*myProblem,nullptr,'W');
      break;

    default:
      std::cerr <<" *** SIMbase::project: Projection method "<< pMethod
                <<" not implemented."<< std::endl;
    }

    if (!ok)
    {
      std::cerr <<" *** SIMbase::project: Failure when projecting patch "
                << myModel[i]->idx+1 <<"."<< std::endl;
      return false;
    }

    size_t nComps = values.rows();
    size_t nNodes = values.cols();
    if (ssol.empty())
      ssol.resize(nComps,ngNodes);

    // Nodal averaging for nodes referred to by two or more patches
    // (these are typically the interface nodes)
    for (n = 1; n <= nNodes; n++)
      if (count.empty())
	ssol.fillColumn(myModel[i]->getNodeID(n),values.getColumn(n));
      else
      {
	int inod = myModel[i]->getNodeID(n);
	for (j = 1; j <= nComps; j++)
	  ssol(j,inod) += values(j,n);
	count(inod) ++;
      }
  }

  // Divide through by count(n) to get the nodal average at the interface nodes
  for (n = 1; n <= count.size(); n++)
    if (count(n) > 1.0)
      for (j = 1; j <= ssol.rows(); j++)
	ssol(j,n) /= count(n);

  return true;
}


bool SIMbase::project (Vector& values, const RealFunc* f,
                       int basis, int iField, int nFields, double time) const
{
  bool ok = true;
  for (size_t j = 0; j < myModel.size() && ok; j++)
  {
    if (myModel[j]->empty()) continue; // skip empty patches

    Vector loc_scalar;
    ok = myModel[j]->evaluate(f,loc_scalar,basis,time);

    if (nFields <= 1)
      ok &= myModel[j]->injectNodeVec(loc_scalar,values,1,basis);
    else
    {
      // Interleave 
      size_t i, k = iField;
      Vector loc_vector(loc_scalar.size()*nFields);
      myModel[j]->extractNodeVec(values,loc_vector,0,basis);
      for (i = 0; i < loc_scalar.size(); i++, k += nFields)
        loc_vector[k] = loc_scalar[i];
      ok &= myModel[j]->injectNodeVec(loc_vector,values,0,basis);
    }
  }

  return ok;
}


bool SIMbase::extractPatchSolution (IntegrandBase* problem,
                                    const Vectors& sol, size_t pindx) const
{
  ASMbase* pch = this->getPatch(pindx+1);
  if (!pch) return false;

  problem->initNodeMap(pch->getGlobalNodeNums());
  for (size_t i = 0; i < sol.size() && i < problem->getNoSolutions(); i++)
    if (!sol[i].empty())
      pch->extractNodeVec(sol[i],problem->getSolution(i),mySam->getMADOF());

  return this->extractPatchDependencies(problem,myModel,pindx);
}


size_t SIMbase::extractPatchSolution (const Vector& sol, Vector& vec,
                                      int pindx, unsigned char nndof) const
{
  ASMbase* pch = pindx >= 0 ? this->getPatch(pindx+1) : nullptr;
  if (!pch || sol.empty()) return 0;

  pch->extractNodeVec(sol,vec,nndof);

  return (nndof > 0 ? nndof : pch->getNoFields(1))*pch->getNoNodes(1);
}


bool SIMbase::injectPatchSolution (Vector& sol, const Vector& vec,
                                   int pindx, unsigned char nndof) const
{
  ASMbase* pch = pindx >= 0 ? this->getPatch(pindx+1) : nullptr;

  return pch ? pch->injectNodeVec(vec,sol,nndof) : false;
}


bool SIMbase::evalSecondarySolution (Matrix& field, int pindx) const
{
  ASMbase* pch = pindx >= 0 ? this->getPatch(pindx+1) : nullptr;
  if (!pch) return false;

  const_cast<SIMbase*>(this)->setPatchMaterial(pindx+1);
  return pch->evalSolution(field,*myProblem);
}


bool SIMbase::extractPatchElmRes (const Matrix& glbRes, Matrix& elRes,
				  int pindx) const
{
  ASMbase* pch = pindx >= 0 ? this->getPatch(pindx+1) : nullptr;
  if (!pch || glbRes.empty()) return false;

  pch->extractElmRes(glbRes,elRes);
  return true;
}


bool SIMbase::setPatchMaterial (size_t patch)
{
  for (PropertyVec::const_iterator p = myProps.begin(); p != myProps.end(); ++p)
    if (p->pcode == Property::MATERIAL && p->patch == patch)
      return this->initMaterial(p->pindx);

  return false;
}


bool SIMbase::refine (const LR::RefineData& prm, const char* fName)
{
  Vectors svec;
  return this->refine(prm,svec,fName);
}


bool SIMbase::refine (const LR::RefineData& prm,
                      Vector& sol, const char* fName)
{
  if (sol.empty())
    return this->refine(prm,fName);

  Vectors svec(1,sol);
  bool result = this->refine(prm,svec,fName);
  sol = svec.front();
  return result;
}


/*!
  \note Solution transfer is available for single-patch models only.
  Therefore, we assume there is a one-to-one correspondance between the
  patch-level and global-level node numbering, and we don't need to use
  ASMbase::[extract|inject]NodeVec(), which in any case would not work
  as we would have needed to regenerate the global node numbering data first.

  TODO: If/when going adaptive multi-patch, the solution transfer has to be
  done after preprocess().
*/

bool SIMbase::refine (const LR::RefineData& prm,
                      Vectors& sol, const char* fName)
{
  isRefined = false;
  ASMunstruct* pch = nullptr;
  for (size_t i = 0; i < myModel.size(); i++)
    if ((pch = dynamic_cast<ASMunstruct*>(myModel[i])))
    {
      if (isRefined && !sol.empty())
      {
        std::cerr <<" *** SIMbase::refine: Solution transfer is not implemented"
                  <<" for multi-patch models.\n";
        return false;
      }

      if (!pch->refine(prm,sol,fName))
        return false;

      isRefined = true;
    }

  return isRefined;
}
