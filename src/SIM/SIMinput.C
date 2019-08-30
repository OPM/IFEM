// $Id$
//==============================================================================
//!
//! \file SIMinput.C
//!
//! \date Sep 24 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Sub-class with functionality for model input and setup.
//!
//==============================================================================

#include "SIMinput.h"
#include "SIMoptions.h"
#include "ModelGenerator.h"
#include "ASMstruct.h"
#include "ASMunstruct.h"
#ifdef HAS_LRSPLINE
#include "ASMLRSpline.h"
#endif
#include "IntegrandBase.h"
#include "GlbL2projector.h"
#include "LinSolParams.h"
#include "DualField.h"
#include "Functions.h"
#include "FunctionSum.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include "HDF5Reader.h"
#include "IFEM.h"
#include "tinyxml.h"
#include <fstream>
#include <sstream>
#include <numeric>


std::istream* SIMinput::getPatchStream (const char* tag, const char* patch)
{
  if (!strcasecmp(tag+5,"file"))
  {
    IFEM::cout <<"\tReading data file "<< patch << std::endl;
    return new std::ifstream(patch);
  }
  else if (!strcasecmp(tag+5,"es") || strlen(tag) == 5)
  {
    IFEM::cout <<"\tReading inlined patch geometry definition"<< std::endl;
    // Replace all '\' and '|' characters in the string by newline '\n'
    for (size_t i = 0; i < strlen(patch); i++)
      if (patch[i] == '\\' || patch[i] == '|')
	const_cast<char&>(patch[i]) = '\n';
    return new std::stringstream(patch);
  }
  else
    return nullptr;
}


bool SIMinput::readPatches (std::istream& isp, const char* whiteSpace)
{
  unsigned char maxSpaceDim = 0;
  for (int pchInd = 0; isp.good(); pchInd++)
  {
    ASMbase* pch = this->readPatch(isp,pchInd,CharVec(),whiteSpace);
    if (pch)
    {
      myModel.push_back(pch);
      if (pch->getNoSpaceDim() > maxSpaceDim)
        maxSpaceDim = pch->getNoSpaceDim();
    }
    else if (this->getLocalPatchIndex(pchInd+1) > 0)
      return false;
  }

  // Reset number of space dimensions if all patches have less than nsd
  if (maxSpaceDim > 0 && maxSpaceDim < nsd)
  {
    IFEM::cout <<"  Resetting number of space dimensions to "<< (int)maxSpaceDim
               <<" to match patch file dimensionality."<< std::endl;
    nsd = maxSpaceDim;
  }

  return true;
}


bool SIMinput::parseGeometryTag (const TiXmlElement* elem)
{
  IFEM::cout <<"  Parsing <"<< elem->Value() <<">"<< std::endl;

  if (!strncasecmp(elem->Value(),"patch",5) && elem->FirstChild())
  {
    if (!myModel.empty())
      return true; // We already have a model, skip geometry definition

    const char* patch = elem->FirstChild()->Value();
    std::istream* isp = getPatchStream(elem->Value(),patch);
    if (isp)
    {
      this->readPatches(*isp,"\t");
      delete isp;
    }
    else
      return true;

    if (myModel.empty())
    {
      std::cerr <<" *** SIMinput::parse: No patches read."<< std::endl;
      return false;
    }
    if (myPatches.empty())
      nGlPatches = myModel.size();
  }

  else if (!strcasecmp(elem->Value(),"nodefile") && elem->FirstChild())
  {
    if (!this->createFEMmodel())
      return false;

    const char* file = elem->FirstChild()->Value();
    if (strstr(file,".gno"))
    {
      std::ifstream isn(file);
      if (!isn.good())
      {
        std::cerr <<" *** SIMinput::parse: Failure opening input file \""
                  << file <<"\"."<< std::endl;
        return false;
      }

      IFEM::cout <<"\tReading data file "<< file << std::endl;
      this->readNodes(isn);
    }
    else if (strstr(file,".hdf5"))
    {
      IFEM::cout <<"\tReading global node numbers from "<< file << std::endl;
      HDF5Reader hdf5(file,ProcessAdm());
      const char* field = elem->Attribute("field");
      for (int i = 1; i <= nGlPatches; i++)
      {
        IntVec nodes;
        ASMbase* pch = this->getPatch(i,true);
        std::stringstream str;
        str << "/0/" << (field ? field : "node numbers") << "/" << i;
        if (pch && hdf5.readVector(str.str(),nodes))
          pch->setNodeNumbers(nodes);
      }
    }
  }

  else if (!strcasecmp(elem->Value(),"partitioning"))
  {
    int proc = 0;
    if (!utl::getAttribute(elem,"procs",proc))
      return false;
    else if (proc != adm.getNoProcs()) // silently ignore
      return true;
    IFEM::cout <<"\tNumber of partitions: "<< proc << std::endl;

    const TiXmlElement* part = elem->FirstChildElement("part");
    if (part) nGlPatches = 0;
    for (; part; part = part->NextSiblingElement("part"))
    {
      int first = -2, last = -2;
      utl::getAttribute(part,"proc",proc);
      utl::getAttribute(part,"lower",first);
      utl::getAttribute(part,"upper",last);

      if (last > nGlPatches)
        nGlPatches = last;

      if (proc == adm.getProcId())
      {
        myPatches.reserve(last-first+1);
        for (int j = first; j <= last && j > -1; j++)
          myPatches.push_back(j);
      }
      for (int i = first; i <= last; i++)
        adm.dd.setPatchOwner(i,proc);
    }

    // If equal number of blocks per processor
    if (myPatches.empty() && utl::getAttribute(elem,"nperproc",proc))
    {
      for (int j = 1; j <= proc; j++)
        myPatches.push_back(adm.getProcId()*proc+j);
      nGlPatches = adm.getNoProcs()*proc;
      for (int i = 1; i <= nGlPatches; i++)
        adm.dd.setPatchOwner(i,(i-1)/proc);
    }

    std::string file;
    if (myPatches.empty() && utl::getAttribute(elem,"file",file))
    {
      std::ifstream ifs(file, std::ios_base::in|std::ios_base::binary);
       std::vector<int> elms;
      if (ifs.good())
      {
        IFEM::cout << "\tReading partitioning from file " << file;
        std::vector<int> elmOfs(adm.getNoProcs());

        size_t ofs = 0;
        for (int i = 0; i < adm.getNoProcs(); ++i) {
          ifs.read(reinterpret_cast<char*>(&elmOfs[i]), sizeof(int));
          if (i < adm.getProcId())
            ofs += elmOfs[i];
        }
        ifs.seekg(ofs*sizeof(int), std::ios_base::cur);
        elms.resize(elmOfs[adm.getProcId()]);
        for (size_t i = 0; i < elms.size(); ++i)
	  ifs.read(reinterpret_cast<char*>(&elms[i]), sizeof(int));

        IFEM::cout << ", size = " << elms.size() << std::endl;
      }
      bool save = false;
      utl::getAttribute(elem,"save",save);
      adm.dd.setElms(elms, save ? file : "");
    }
  }

  else if (!strcasecmp(elem->Value(),"topologysets"))
  {
    std::string name, type;
    const TiXmlElement* set = elem->FirstChildElement("set");
    for (; set; set = set->NextSiblingElement("set"))
      if (utl::getAttribute(set,"name",name))
      {
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
        else if (type == "nodes")
          idim = 4;
        else
          utl::getAttribute(set,"dimension",idim);
        if (idim > 0 && idim < 4 && utl::getAttribute(set,"closure",type,true))
          if (type == "open") idim = -idim; // i.e., excluding its boundary

        TopEntity& top = myEntitys[name];
        const TiXmlElement* item = set->FirstChildElement("item");
        for (; item; item = item->NextSiblingElement("item"))
        {
          int patch = 1;
          utl::getAttribute(item,"patch",patch);
          if ((patch = this->getLocalPatchIndex(patch)) > 0)
          {
            ASMbase* pch = nullptr;
            if (abs(idim) == (int)this->getNoParamDim())
              top.insert(TopItem(patch,0,idim));
            else if (idim == 4 && (pch = this->getPatch(patch)))
            {
              if (item->FirstChild())
                utl::parseIntegers(pch->getNodeSet(name),
                                   item->FirstChild()->Value());
              top.insert(TopItem(patch,pch->getNodeSetIdx(name),idim));
            }
            else if (item->FirstChild())
            {
              IntVec items;
              utl::parseIntegers(items,item->FirstChild()->Value());
              for (int item : items)
                top.insert(TopItem(patch,item,idim));
            }
          }
        }
      }

    if (!myEntitys.empty())
    {
      IFEM::cout <<"\tTopology sets: ";
      TopologySet::const_iterator it;
      for (it = myEntitys.begin(); it != myEntitys.end(); ++it)
      {
        if (it != myEntitys.begin())
          IFEM::cout <<"\t               ";
        IFEM::cout << it->first;
        for (const TopItem& it2 : it->second)
          IFEM::cout << it2;
        IFEM::cout << std::endl;
      }
    }
  }

  return true;
}


/*!
  This method can not be part of the SIMinput::parseGeometryTag() method,
  since the periodicities must be processed after the patches have been
  order-elevated and/or refined. The method is instead invoked from the
  dimension-specific sub-classes which handle the refinements.
*/

bool SIMinput::parsePeriodic (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"periodic") || !this->createFEMmodel())
    return false;

  int patch = 0, pedir = 1;
  utl::getAttribute(elem,"patch",patch);
  utl::getAttribute(elem,"dir",pedir);
  if (patch < 1 || patch > nGlPatches)
  {
    std::cerr <<" *** SIMinput::parse: Invalid patch index "
              << patch << std::endl;
    return false;
  }

  ASMbase* pch = this->getPatch(patch,true);
  if (pch)
  {
    IFEM::cout <<"\tPeriodic "<< char('H'+pedir) <<"-direction P"<< patch
               << std::endl;
    pch->closeBoundaries(pedir);
  }

  return true;
}


//! \brief Integer value flagging global axes in boundary conditions.
#define GLOBAL_AXES     -1
//! \brief Integer value flagging local axes in boundary conditions.
#define LOCAL_AXES      -2
//! \brief Integer value flagging local projected axes in boundary conditions.
#define LOCAL_PROJECTED -3

bool SIMinput::parseBCTag (const TiXmlElement* elem)
{
  IFEM::cout <<"  Parsing <"<< elem->Value() <<">"<< std::endl;

  if (!strcasecmp(elem->Value(),"propertyfile") && elem->FirstChild())
  {
    const char* file = elem->FirstChild()->Value();
    std::ifstream isp(file);
    if (!isp)
    {
      std::cerr <<" *** SIMinput::parseBCTag: Failure opening input file \""
                << file <<"\"."<< std::endl;
      return false;
    }
    IFEM::cout <<"\tReading data file "<< file << std::endl;
    while (isp.good())
    {
      Property p;
      int ldim, lindx = 0;
      isp >> p.pindx >> p.patch >> ldim;
      if (ldim < (int)this->getNoParamDim())
        isp >> lindx;

      // We always require the item indices to be 1-based
      p.ldim = ldim;
      p.lindx = 1+lindx;
      p.patch = this->getLocalPatchIndex(1+p.patch);
      if (p.patch > 0 && isp.good())
        myProps.push_back(p);
    }
  }

  else if (!strcasecmp(elem->Value(),"propertycodes"))
  {
    const TiXmlElement* code = elem->FirstChildElement("code");
    for (; code; code = code->NextSiblingElement())
    {
      int icode = 0;
      utl::getAttribute(code,"value",icode);
      const TiXmlElement* patch = code->FirstChildElement("patch");
      for (; patch; patch = patch->NextSiblingElement("patch"))
      {
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

  else if (!strcasecmp(elem->Value(),"neumann"))
  {
    const TiXmlNode* nval = elem->FirstChild();

    std::string set, type;
    utl::getAttribute(elem,"set",set);
    utl::getAttribute(elem,"type",type,true);

    int ndir = 0, order = 1;
    utl::getAttribute(elem,"direction",ndir);
    utl::getAttribute(elem,"order",order);

    size_t i = myProps.size();
    int code = this->getUniquePropertyCode(set);
    if (code == 0) utl::getAttribute(elem,"code",code);
    IFEM::cout <<"\tNeumann code "<< code;

    if (order > 1)
    {
      IFEM::cout <<" order "<< order;
      // Flag Neumann order by increasing the topology index by 10
      for (; i < myProps.size(); i++)
        myProps[i].lindx += 10*(order-1);
    }

    if (type == "anasol")
    {
      IFEM::cout <<" (analytic)";
      this->setPropertyType(code,Property::NEUMANN_ANASOL);
    }
    else if (type == "generic")
    {
      IFEM::cout <<" (generic)";
      this->setPropertyType(code,Property::NEUMANN_GENERIC);
    }
    else if (nval)
    {
      IFEM::cout <<" direction "<< ndir;
      if (!type.empty())
        IFEM::cout <<" ("<< type <<")";
      this->setNeumann(nval->Value(),type,ndir,code);
    }
    else
      this->setNeumann(std::string(),type,ndir,code);

    IFEM::cout << std::endl;
  }

  else if (!strcasecmp(elem->Value(),"dirichlet") && !ignoreDirichlet)
  {
    const TiXmlNode* dval = nullptr;
    int comp = 0, symm = 0, basis = 1;
    std::string set, type, axes;
    utl::getAttribute(elem,"set",set);
    utl::getAttribute(elem,"type",type,true);
    utl::getAttribute(elem,"basis",basis);
    // Handle some predefined property codes for symmtry-conditions (C1-patches)
    if (type == "symmxy" || type == "symmyx")
      comp = symm = 12000;
    else if (type == "symmyz" || type == "symmzy")
      comp = symm = 23000;
    else if (type == "symmzx" || type == "symmxz")
      comp = symm = 31000;
    else if (type == "clamped")
      comp = 123123;
    else
    {
      dval = elem->FirstChild();
      utl::getAttribute(elem,"axes",axes,true);
      utl::getAttribute(elem,"component",comp);
      utl::getAttribute(elem,"comp",comp);
    }
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

    IFEM::cout <<"\tDirichlet code "<< code;
    if (type == "anasol")
    {
      IFEM::cout <<": (analytic)";
      this->setPropertyType(code,Property::DIRICHLET_ANASOL,-1,basis);
    }
    else if (dval && (type == "expression" || type == "field" ||
                      atof(dval->Value()) != 0.0))
    {
      if (!type.empty())
        IFEM::cout <<" ("<< type <<")";
      this->setPropertyType(code,Property::DIRICHLET_INHOM,comp,basis);
      RealFunc* f = utl::parseRealFunc(dval->Value(),type);
      if (!f)
        return false;
      else
        myScalars[abs(code)] = f;
    }
    else
    {
      IFEM::cout <<": (fixed)";
      this->setPropertyType(code,Property::DIRICHLET,comp,basis);
    }
    if (symm)
    {
      code = this->getUniquePropertyCode(set,1+(1+symm/10000)%3);
      IFEM::cout <<"\n\tDirichlet code "<< code <<": (fixed)";
      this->setPropertyType(code,Property::DIRICHLET,comp,basis);
    }
    IFEM::cout << std::endl;
  }

  else if (!strcasecmp(elem->Value(),"robin"))
  {
    std::string set, type, prop;
    utl::getAttribute(elem,"set",set);
    utl::getAttribute(elem,"type",type,true);
    const TiXmlNode* nval = elem->FirstChild();
    if (nval) prop = nval->Value();

    int code = this->getUniquePropertyCode(set);
    if (code == 0) utl::getAttribute(elem,"code",code);

    IFEM::cout <<"\tRobin code "<< code;
    if (!type.empty())
      IFEM::cout <<" ("<< type <<")";
    this->setPropertyType(code,Property::ROBIN);
    this->setNeumann(prop,type,1-this->getNoFields(),code);
    IFEM::cout << std::endl;
  }

  return true;
}


bool SIMinput::parseICTag (const TiXmlElement* elem)
{
  std::string field;
  if (!utl::getAttribute(elem,"field",field))
    return false;

  ICInfo info(field);
  std::string type("file"), file("nofile");
  utl::getAttribute(elem,"type",type);
  if (type == "file")
  {
    if (!utl::getAttribute(elem,"file",file) && elem->FirstChild())
      file = elem->FirstChild()->Value();
    utl::getAttribute(elem,"file_field",info.file_field);
    utl::getAttribute(elem,"file_level",info.file_level);
    utl::getAttribute(elem,"file_basis",info.file_basis);
    utl::getAttribute(elem,"geo_level",info.geo_level);
  }
  else if (elem->FirstChild()) // function
  {
    info.function = elem->FirstChild()->Value();
    info.file_field = type;
    // Both long and short form supported, short prioritized
    utl::getAttribute(elem,"component",info.component);
    utl::getAttribute(elem,"comp",info.component);
  }

  int sim_level = -1;
  if (utl::getAttribute(elem,"level",sim_level) && sim_level > -1)
    info.sim_field += '1'+char(sim_level);

  utl::getAttribute(elem,"basis",info.basis);
  myICs[file].push_back(info);

  IFEM::cout <<"\tInitial condition";
  if (type == "file")
    IFEM::cout <<" file: "<< file;
  else
    IFEM::cout <<" "<< type <<" function: "<< info.function;

  IFEM::cout <<"\n\tField name: \""<< info.sim_field <<"\"";
  if (info.basis > 0)
    IFEM::cout <<" on basis "<< (int)info.basis;
  if (type == "file")
    IFEM::cout <<" (on file \""<< info.file_field
               <<"\" at time level "<< info.geo_level <<")";
  else
    IFEM::cout <<" (component "<< (int)info.component <<")";

  if (sim_level > -1 || type == "file")
    IFEM::cout <<"\n\tTime level: "<< sim_level;
  if (type == "file")
    IFEM::cout <<" (on file "<< info.file_level <<")";

  IFEM::cout << std::endl;
  return true;
}


bool SIMinput::parseLinSolTag (const TiXmlElement* elem)
{
  if (!strcasecmp(elem->Value(),"class"))
    if (elem->FirstChild())
      opt.setLinearSolver(elem->FirstChild()->Value());

  return true;
}


static bool noDumpDataYet = true; //!< To read only once in adaptive loops

bool SIMinput::parseOutputTag (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"dump_lhs_matrix") &&
      strcasecmp(elem->Value(),"dump_rhs_vector") &&
      strcasecmp(elem->Value(),"dump_sol_vector"))
    return opt.parseOutputTag(elem);

  if (elem->FirstChild() && noDumpDataYet)
  {
    IntVec steps;
    DumpData dmp;
    std::string format;
    utl::getAttribute(elem,"format",format);
    utl::getAttribute(elem,"step",steps);
    utl::getAttribute(elem,"eps",dmp.eps);
    dmp.step.insert(steps.begin(),steps.end());
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


FunctionBase* SIMinput::parseDualTag (const TiXmlElement* elem, int ftype)
{
  IFEM::cout <<"  Parsing <"<< elem->Value() <<">";

  int patch = 0, comp = 1;
  double depth = 1.0, width = 0.0;
  Vec3 X0, normal(1.0,0.0,0.0);
  utl::getAttribute(elem,"patch",patch);
  utl::getAttribute(elem,"comp",comp);
  utl::getAttribute(elem,"X0",X0);
  utl::getAttribute(elem,"normal",normal);
  utl::getAttribute(elem,"depth",depth);
  utl::getAttribute(elem,"width",width);

  // In adaptive analysis, reduce the depth by a factor of 0.5 until dmin
  double dmin = depth; int skip = 1;
  utl::getAttribute(elem,"dmin",dmin);
  utl::getAttribute(elem,"skip",skip);
  for (int ref = skip; ref <= isRefined && depth > dmin; ref += skip)
    depth *= 0.5;

  if (patch > 0)
    IFEM::cout <<"\n\tpatch  = "<< patch;
  IFEM::cout <<"\n\tX0     = "<< X0
             <<"\n\tnormal = "<< normal
             <<"\n\tdepth  = "<< depth <<", width = "<< width;
  if (ftype > 1)
    IFEM::cout <<"\n\tcomp   = "<< comp;
  IFEM::cout << std::endl;

  if (nsd == 2)
  {
    if (ftype == 1)
      extrFunc.push_back(new DualRealFunc(X0,normal,depth,width,patch));
    else
      extrFunc.push_back(new DualVecFunc(comp,X0,normal,depth,width,patch));
  }
  else if (nsd == 3)
  {
    Vec3 XZp(0.0,0.0,1.0);
    utl::getAttribute(elem,"XZp",XZp);
    if (ftype == 1)
      extrFunc.push_back(new DualRealFunc(X0,normal,XZp,depth,width,patch));
    else
      extrFunc.push_back(new DualVecFunc(comp,X0,normal,XZp,depth,width,patch));
  }

  double weight = 0.0;
  utl::getAttribute(elem,"weight",weight);
  if (weight != 0.0)
  {
    if (dualField)
      static_cast<FunctionSum*>(dualField)->add(extrFunc.back(),weight);
    else
      dualField = new FunctionSum(extrFunc.back(),weight);
  }

  return extrFunc.back();
}


bool SIMinput::parse (const TiXmlElement* elem)
{
  bool result = true;
  if (!strcasecmp(elem->Value(),"discretization"))
    result = opt.parseDiscretizationTag(elem);
  else if (!strcasecmp(elem->Value(),"restart"))
    result = opt.parseRestartTag(elem);
  else if (!strcasecmp(elem->Value(),"initialcondition"))
    result = this->parseICTag(elem);
  else if (!strcasecmp(elem->Value(),"linearsolver"))
  {
    std::string solver;
    if (utl::getAttribute(elem,"class",solver,true))
      opt.setLinearSolver(solver);
    if (utl::getAttribute(elem,"l2class",solver,true) && solver == "petsc")
      GlbL2::MatrixType = SystemMatrix::PETSC;
  }
  else if (!strcasecmp(elem->Value(),"eigensolver"))
    utl::getAttribute(elem,"mode",opt.eig);
  else if (!strcasecmp(elem->Value(),"postprocessing"))
    noDumpDataYet = lhsDump.empty() && rhsDump.empty();
  else
    result = this->SIMbase::parse(elem);

  // Create the default geometry if no patchfile is specified
  if (!strcasecmp(elem->Value(),"geometry") && this->getNoParamDim() > 0)
    if (!elem->FirstChildElement("patchfile") &&
        !elem->FirstChildElement("patches") &&
        !elem->FirstChildElement("patch"))
    {
      if (myModel.empty())
      {
        const TiXmlElement* part = elem->FirstChildElement("partitioning");
        for (; part; part = part->NextSiblingElement("partitioning"))
          result &= this->parseGeometryTag(part);
      }

      myGen = this->getModelGenerator(elem);
      if (!myGen)
        return false;

      if (myModel.empty())
        myGen->createGeometry(*this);
      if (myPatches.empty())
        nGlPatches = myModel.size();

      myGen->createTopologySets(*this);
    }

  // Check if a characteristic model size is specified
  if (!strcasecmp(elem->Value(),"geometry"))
    utl::getAttribute(elem,"modelsize",ASMbase::modelSize);

  if (!strcasecmp(elem->Value(),"linearsolver"))
  {
    if (!mySolParams)
    {
      if (myProblem)
        mySolParams = new LinSolParams(myProblem->getLinearSystemType());
      else
        mySolParams = new LinSolParams();
    }
    result &= mySolParams->read(elem);
    if (GlbL2::MatrixType == SystemMatrix::PETSC)
    {
      // For now use same solver parameters in the L2-projection
      myGl2Params = new LinSolParams(*mySolParams,LinAlg::SYMMETRIC);
      GlbL2::SolverParams = myGl2Params;
    }
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
    else if (!strcasecmp(elem->Value(),"console"))
      result &= this->opt.parseConsoleTag(child);
    else if (!strcasecmp(elem->Value(),"discretization"))
      result &= opt.parseDiscretizationTag(child);

  return result;
}


int SIMinput::parseMaterialSet (const TiXmlElement* elem, int mindex)
{
  std::string setName;
  utl::getAttribute(elem,"set",setName);
  int code = this->getUniquePropertyCode(setName);

  if (code == 0)
    utl::getAttribute(elem,"code",code);

  if (code > 0)
    if (this->setPropertyType(code,Property::MATERIAL,mindex))
      return mindex;

  return 0;
}


bool SIMinput::parseTopologySet (const TiXmlElement* elem, IntVec& patches)
{
  std::string setName;
  if (utl::getAttribute(elem,"set",setName))
  {
    auto tit = myEntitys.find(setName);
    if (tit == myEntitys.end())
    {
      std::cerr <<" *** SIMinput::parseTopologySet: Undefined topology set \""
                << setName <<"\"."<< std::endl;
      return false;
    }

    patches.clear();
    for (const TopItem& top : tit->second)
      if (top.idim == (short int)this->getNoParamDim())
        patches.push_back(top.patch);
    if (!patches.empty())
      return true;

    std::cerr <<" *** SIMinput::parseTopologySet: Invalid topology set \""
              << setName <<"\" (no patches in this set)."<< std::endl;
    return false;
  }

  int lowpatch = 1, uppatch = 1;
  if (utl::getAttribute(elem,"patch",lowpatch))
    uppatch = lowpatch;
  if (utl::getAttribute(elem,"lowerpatch",lowpatch))
    uppatch = myModel.size();
  utl::getAttribute(elem,"upperpatch",uppatch);

  if (lowpatch < 1 || uppatch > nGlPatches)
  {
    std::cerr <<" *** SIMinput::parseTopologySet: Invalid patch indices, lower="
              << lowpatch <<" upper="<< uppatch << std::endl;
    return false;
  }

  patches.resize(uppatch-lowpatch+1);
  std::iota(patches.begin(),patches.end(),lowpatch);
  return true;
}


const char** SIMinput::getPrioritizedTags () const
{
  // Tags to be parsed first, and in the order specified
  static const char* special[] = { "console", "discretization", "geometry", 0 };
  return special;
}


bool SIMinput::parse (char* keyWord, std::istream& is)
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
          std::cerr <<" *** SIMinput: Failure opening patch file \""
                    << cline <<"\"."<< std::endl;
      }

      if ((int)myModel.size() < npatch)
      {
        std::cerr <<" *** SIMinput::parse: Expected "<< npatch
                  <<" patches but could read only "<< myModel.size()
                  <<"."<< std::endl;
        return false;
      }
    }
    else // just read through the npatch next lines without doing anything
      for (int i = 0; i < npatch && utl::readLine(is); i++);
  }

  else if (!strncasecmp(keyWord,"PATCHFILE",9))
  {
    if (!myModel.empty())
      return true;

    size_t i = 9; while (i < strlen(keyWord) && isspace(keyWord[i])) i++;
    IFEM::cout <<"\nReading data file "<< keyWord+i << std::endl;
    std::ifstream isp(keyWord+i);
    this->readPatches(isp);

    if (myModel.empty())
    {
      std::cerr <<" *** SIMinput::parse: No patches read."<< std::endl;
      return false;
    }
  }

  else if (!strncasecmp(keyWord,"NODEFILE",8))
  {
    if (!this->createFEMmodel())
      return false;

    bool oneBasedIdx = keyWord[8] == '1';
    size_t i = (oneBasedIdx || keyWord[8] == '0') ? 9 : 8;
    while (i < strlen(keyWord) && isspace(keyWord[i])) i++;

    std::stringstream filenames(keyWord+i);
    for (size_t basis = 0; filenames.good(); basis++)
    {
      std::string filename;
      filenames >> filename;
      if (basis == 0 && filenames.good())
        basis = 1;

      std::ifstream isn(filename.c_str());
      if (!isn.good())
      {
        std::cerr <<" *** SIMinput::parse: Failure opening input file \""
                  << filename <<"\"."<< std::endl;
        return false;
      }

      if (basis < 2)
        IFEM::cout <<"\n";
      IFEM::cout <<"Reading data file "<< filename;
      if (basis > 0)
        IFEM::cout <<" (basis "<< basis <<")";
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
          std::cerr <<" *** SIMinput::parse: Failed to assign node numbers"
                    <<" for patch "<< patch <<"."<< std::endl;
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
      std::cerr <<" *** SIMinput::parse: Failure opening input file \""
                << std::string(keyWord+i) <<"\"."<< std::endl;
      return false;
    }

    IFEM::cout <<"\nReading data file "<< keyWord+i << std::endl;
    while (isp.good())
    {
      Property p;
      int ldim, lindx = 0;
      isp >> p.pindx >> p.patch >> ldim;
      if (ldim < (int)this->getNoParamDim())
        isp >> lindx;

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

  else if (!strncasecmp(keyWord,"PERIODIC",8))
  {
    if (!this->createFEMmodel())
      return false;

    int nper = atoi(keyWord+8);
    IFEM::cout <<"\nNumber of periodicities: "<< nper << std::endl;
    for (int i = 0; i < nper && (cline = utl::readLine(is)); i++)
    {
      int patch = atoi(strtok(cline," "));
      int pedir = (cline = strtok(nullptr," ")) ? atof(cline) : 1;
      if (patch < 1 || patch > (int)myModel.size())
      {
        std::cerr <<" *** SIMinput::parse: Invalid patch index "
                  << patch << std::endl;
        return false;
      }
      IFEM::cout <<"\tPeriodic "<< char('H'+pedir) <<"-direction P"<< patch
                 << std::endl;
      myModel[patch-1]->closeBoundaries(pedir);
    }
  }

  else if (!strncasecmp(keyWord,"LINEARSOLVER",12))
  {
    std::cerr <<" *** SIMinput::parse: Reading linear solver parameters"
              <<" from flat file format is depreciated."<< std::endl;
    return false;
  }

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

      if (last > nGlPatches)
        nGlPatches = last;

      if (proc == myPid && last >= first)
      {
        myPatches.reserve(last-first+1);
        for (int j = first; j <= last; j++)
          myPatches.push_back(j);
      }
    }

    if (myPatches.empty())
    {
      std::cerr <<" *** SIMinput::parse: No partitioning input for processor "
                << myPid <<"."<< std::endl;
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
    std::cerr <<"  ** SIMinput::parse: Unknown keyword \""
              << keyWord <<"\"."<< std::endl;
#endif

  return true;
}


bool SIMinput::createFEMmodel (char resetNumb)
{
  if (resetNumb) ASMbase::resetNumbering();

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

  return true;
}


int SIMinput::getUniquePropertyCode (const std::string& setName, int code)
{
  if (setName.empty()) return 0;

  int cinc = code < 0 ? -1000000 : 1000000;
  if (code == 0) code = cinc;
  PropertyVec::const_iterator pit = myProps.begin();
  for (int trial = 0; pit != myProps.end(); trial++)
  {
    if (trial > 0) code += cinc;
    pit = std::find_if(myProps.begin(),myProps.end(),
                       [code](const Property& p)
                       { return abs(p.pindx) == abs(code); });
  }

  return this->createPropertySet(setName,code) ? code : 0;
}


bool SIMinput::createPropertySet (const std::string& setName, int pc)
{
  if (setName.empty()) return true;

  auto tit = myEntitys.find(setName);
  if (tit == myEntitys.end())
  {
    std::cerr <<" *** SIMinput::createPropertySet: Undefined topology set \""
              << setName <<"\"."<< std::endl;
    return false;
  }

  auto pit = std::find_if(myProps.begin(),myProps.end(),
                         [pc](const Property& p){ return abs(p.pindx) == pc; });
  if (pit != myProps.end())
  {
    std::cerr <<" *** SIMinput::createPropertySet: Duplicated property code "
              << pc <<"."<< std::endl;
    return false;
  }

  // Create the actual property objects that are used during simulation
  for (const TopItem& top : tit->second)
    myProps.push_back(Property(Property::UNDEFINED,pc,
                               top.patch,top.idim,top.item));

  return true;
}


/*!
  Negative values on \a pindex are used to flag the use of local axes for
  Dirichlet boundary conditions as follows:
  \a pindex=-2 : Compute local axes directions directly at the Greville points.
  \a pindex=-3 : Compute local axes directions by projecting the local tangent
  direction onto the spline basis, to obtain directions at control points.
*/

size_t SIMinput::setPropertyType (int code, Property::Type ptype,
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


size_t SIMinput::setVecProperty (int code, Property::Type ptype,
                                 VecFunc* field, int pflag)
{
  if (field) myVectors[abs(code)] = field;
  return this->setPropertyType(code,ptype,pflag);
}


bool SIMinput::setTracProperty (int code, Property::Type ptype,
                                TractionFunc* field)
{
  if (field) myTracs[code] = field;
  return this->setPropertyType(code,ptype);
}


bool SIMinput::setNeumann (const std::string& prop, const std::string& type,
                           int direction, int code)
{
  if ((direction == 0 && this->getNoFields() == 1) || direction > nsd)
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


IntVec SIMinput::getFunctionsForElements (const IntVec& elements)
{
  IntSet functions;
#ifdef HAS_LRSPLINE
  for (ASMbase* pch : myModel)
  {
    ASMLRSpline* lrPch = dynamic_cast<ASMLRSpline*>(pch);
    if (lrPch) lrPch->getFunctionsForElements(functions,elements);
  }
#ifdef SP_DEBUG
  size_t j = 0, k = 0;
  std::cout <<"SIMinput::getFunctionsForElements: nel="<< elements.size();
  for (int e : elements)  std::cout << ((++j)%10 == 1 ? '\n' : ' ') << e;
  std::cout <<"\nSIMinput::getFunctionsForElements: nfn="<< functions.size();
  for (int f : functions) std::cout << ((++k)%10 == 1 ? '\n' : ' ') << f;
  std::cout << std::endl;
#endif
#endif
  return IntVec(functions.begin(),functions.end());
}


bool SIMinput::refine (const LR::RefineData& prm)
{
  Vectors svec;
  return this->refine(prm,svec);
}


bool SIMinput::refine (const LR::RefineData& prm, Vector& sol)
{
  if (sol.empty())
    return this->refine(prm);

  Vectors svec(1,sol);
  bool result = this->refine(prm,svec);
  sol = svec.front();
  return result;
}


bool SIMinput::refine (const LR::RefineData& prm, Vectors& sol)
{
  if (myModel.empty())
    return false;

  ASMunstruct* pch = nullptr;
  for (size_t i = 0; i < myModel.size(); i++)
    if (!(pch = dynamic_cast<ASMunstruct*>(myModel[i])))
    {
      std::cerr <<" *** SIMinput::refine: Model is not constructed from"
                <<" unstructured (ASMunstruct) patches."<< std::endl;
      return false;
    }
    else if (nGlPatches == 1)
    {
      // Single patch models only pass refinement call to the ASM level
      if (!pch->refine(prm,sol))
        return false;

      ++isRefined;
      return true;
    }

  if (!prm.errors.empty()) // refinement by true_beta
  {
    std::cerr <<" *** SIMinput::refine: True beta refinement not"
              <<" implemented for multi-patch models."<< std::endl;
    return false;
  }

  // Multi-patch models need to pass refinement indices over patch boundaries
  std::vector<LR::RefineData> prmloc(myModel.size(),LR::RefineData(prm));
  std::vector<IntSet> refineIndices(myModel.size());
  std::vector<IntSet> conformingIndices(myModel.size());
  for (size_t i = 0; i < myModel.size(); i++)
  {
    // Extract local indices from the vector of global indices
    int locId;
    for (int k : prm.elements)
      if ((locId = myModel[i]->getNodeIndex(k+1)) > 0)
        refineIndices[i].insert(locId-1);

    // fetch all boundary nodes covered (may need to pass this to other patches)
    pch = dynamic_cast<ASMunstruct*>(myModel[i]);
    IntVec bndry_nodes = pch->getBoundaryCovered(refineIndices[i]);

    // DESIGN NOTE: It is tempting here to use patch connectivity information.
    // However, this does not account (in the general case)
    // for cross-connections in L-shape geometries, i.e.,
    //
    // +-----+
    // | #1  |
    // |     |         patch #1 (edge 3) connected to patch #2 (edge 4)
    // +-----+-----+   patch #2 (edge 2) connected to patch #3 (edge 1)
    // | #2  | #3  |
    // |     |     |   we need to pass the corner index of patch #3 (vertex 3)
    // +-----+-----+   to patch #1 (vertex 2), but this connection is not
    //                 guaranteed to appear in the input file

    // for all boundary nodes, check if these appear on other patches
    for (int k : bndry_nodes)
    {
      int globId = myModel[i]->getNodeID(k+1);
      for (size_t j = 0; j < myModel.size(); j++)
        if (j != i && (locId = myModel[j]->getNodeIndex(globId)) > 0)
        {
          conformingIndices[j].insert(locId-1);
          conformingIndices[i].insert(k);
        }
    }
  }

  Vectors lsols;
  lsols.reserve(sol.size()*myModel.size());
  size_t ngNodes = this->getNoNodes(1);

  for (size_t i = 0; i < myModel.size(); i++)
  {
    pch = dynamic_cast<ASMunstruct*>(myModel[i]);
    pch->extendRefinementDomain(refineIndices[i],conformingIndices[i]);
    LR::RefineData prmloc(prm);
    prmloc.elements = IntVec(refineIndices[i].begin(),refineIndices[i].end());

    Vectors lsol(sol.size());
    for (size_t j = 0; j < sol.size(); j++)
      myModel[i]->extractNodeVec(sol[j], lsol[j], sol[j].size()/ngNodes);
    if (!pch->refine(prmloc,lsol))
      return false;
    for (const Vector& s : lsol)
      lsols.push_back(s);
  }
  sol = lsols;
  ++isRefined;
  return true;
}


bool SIMinput::setInitialCondition (SIMdependency* fieldHolder,
                                    const std::string& fileName,
                                    const InitialCondVec& info)
{
  HDF5Reader hdf5reader(fileName,adm);

  std::map<std::string,PatchVec> basisMap;

  // Loop over the initial conditions
  for (const ICInfo& it : info)
  {
    // Do we have this field?
    Vector* field = fieldHolder->getField(it.sim_field);
    if (!field) continue;

    // Load basis
    CharVec nf(1,this->getNoFields(it.basis));
    PatchVec& basisVec = basisMap[it.file_basis];
    std::stringstream str;
    str << it.geo_level << "/" << it.file_basis << "/basis";
    int nPatches = hdf5reader.getFieldSize(str.str());
    if (basisVec.empty())
      for (int i = 0; i < nPatches; i++)
        if (this->getLocalPatchIndex(i+1) > 0)
        {
          std::stringstream str, spg2;
          str << it.geo_level << "/" << it.file_basis << "/basis/" << i+1;
          std::string pg2;
          hdf5reader.readString(str.str(),pg2);
          spg2 << pg2;
          basisVec.push_back(this->readPatch(spg2,i,nf));
        }

    // Load result field, patch by patch
    for (int i = 0; i < nPatches; i++)
    {
      int p = this->getLocalPatchIndex(i+1);
      ASMbase* pch = this->getPatch(p);
      if (!pch) continue;

      Vector loc, newloc;
      std::stringstream str;
      str << it.file_level <<"/"
          << it.file_basis <<"/fields/"
          << it.file_field <<"/"<< i+1;
      hdf5reader.readVector(str.str(), loc);
      basisVec[p-1]->copyParameterDomain(pch);
      if (pch->evaluate(basisVec[p-1], loc, newloc, it.basis))
        pch->injectNodeVec(newloc, *field,
                           newloc.size()/pch->getNoNodes(it.basis), it.basis);
    }
  }

  // Clean up basis patches
  for (const auto& itb : basisMap)
    for (ASMbase* pch : itb.second)
      delete pch;

  return true;
}


bool SIMinput::setInitialConditions (SIMdependency* fieldHolder)
{
  if (!fieldHolder)
    fieldHolder = this;

  bool result = true;
  for (const auto& it : myICs)
    if (it.first != "nofile")
      result &= this->setInitialCondition(fieldHolder,it.first,it.second);
    else for (const ICInfo& ic : it.second)
    {
      // Do we have this field?
      Vector* field = fieldHolder->getField(ic.sim_field);
      if (!field) continue;

      IFEM::cout <<"\nLoading initial condition for \""<< ic.sim_field
                 <<"\" component "<< (int)ic.component
                 <<"\nfrom "<< ic.file_field <<" function";
      RealFunc* fn = utl::parseRealFunc(ic.function,ic.file_field);
      IFEM::cout << std::endl;
      result &= this->project(*field,fn,ic.basis,ic.component-1,
                              this->getNoFields(ic.basis));
      delete fn;
    }

  return result;
}


bool SIMinput::hasIC (const std::string& name) const
{
  for (const auto& it : myICs)
    for (const ICInfo& ic : it.second)
      if (ic.sim_field.compare(0,name.size(),name) == 0)
        return true;

  return false;
}


bool SIMinput::deSerialize (const std::map<std::string,std::string>&)
{
  std::cerr <<" *** SIMinput::deSerialize: Must be implemented in sub-class.\n"
            <<"     Restart not supported for "<< this->getName() << std::endl;
  return false;
}


IntMat SIMinput::getElmConnectivities () const
{
  IntMat neigh(this->getNoElms());
  for (const ASMbase* pch : this->getFEModel())
    pch->getElmConnectivities(neigh);

  for (auto& iface : myInterfaces) {
    if (iface.dim != static_cast<int>(this->getNoSpaceDim())-1)
      continue;

    const ASMbase* master = this->getPatch(iface.master);
    const ASMbase* slave = this->getPatch(iface.slave);
    const ASMstruct* pch = dynamic_cast<const ASMstruct*>(master);

    std::vector<int> sElms, mElms;
    slave->getBoundaryElms(iface.sidx,iface.orient,sElms);
    master->getBoundaryElms(iface.midx,0,mElms);

    DomainDecomposition::OrientIterator iter(slave, iface.orient, iface.sidx);
    auto it_n = iter.begin();
    auto it_m = mElms.begin();
    for (size_t i = 0; i < iter.size(); ++i, ++it_m, ++it_n) {
      if (pch) {
        neigh[sElms[*it_n]][iface.sidx-1] = *it_m;
        neigh[*it_m][iface.midx-1] = sElms[*it_n];
      } else {
        neigh[sElms[*it_n]].push_back(*it_m);
        neigh[*it_m].push_back(sElms[*it_n]);
      }
    }
  }

  return neigh;
}
