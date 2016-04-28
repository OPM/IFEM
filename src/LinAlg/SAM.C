// $Id$
//==============================================================================
//!
//! \file SAM.C
//!
//! \date Oct 1 2007
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of FE matrices into system matrices.
//!
//==============================================================================

#include "SAM.h"
#include "SystemMatrix.h"

#ifdef USE_F77SAM
#if defined(_WIN32)
#define elmeq_  ELMEQ
#define syseq_  SYSEQ
#define addev2_ ADDEV2
#define expand_ EXPAND
#elif defined(_AIX)
#define elmeq_  elmeq
#define syseq_  syseq
#define addev2_ addev2
#define expand_ expand
#endif

extern "C" {
//! \brief Finds the matrix of element equation numbers \a MEEN for an element.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void elmeq_(const int* madof, const int* mnpc, const int* mpmceq,
	    const int* meqn, const int& nenod, int* meen,
	    int& nedof, int& neslv, int &neprd);

//! \brief Determines the control matrix \a MEQN.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void syseq_(const int* msc, const int* mpmceq, const int* mmceq,
	    const int& lpu, int* mpar, int* meqn, int& ierr);

//! \brief Adds an element vector \a eS into the system vector \a sysRHS.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void addev2_(const Real* eS, const Real* ttcc, const int* mpar,
	     const int* madof, const int* meqn, const int* mpmnpc,
	     const int* mmnpc, const int* mpmceq, const int* mmceq,
	     const int& iel, const int& nedof, const int& lpu,
	     Real* sysRHS, int* work, int& ierr);

//! \brief Expands and rearranges a system vector from equation to DOF-order.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void expand_(const Real* solVec, const Real* ttcc,
	     const int* mpmceq, const int* mmceq, const int* meqn,
	     const Real& ff, const Real& fs, const int& ndof, const int& neq,
	     Real* sv);
}
#endif


SAM::SAM () : nnod(mpar[0]), nel(mpar[1]), ndof(mpar[2]),
	      nspdof(mpar[5]), nceq(mpar[6]), neq(mpar[10]),
	      nmmnpc(mpar[14]), nmmceq(mpar[15])
{
  // Initialize the parameters array to zero
  memset(mpar,0,sizeof(mpar));

  // Initialize the array pointers
  mpmnpc = 0;
  mmnpc  = 0;
  madof  = 0;
  msc    = 0;
  mpmceq = 0;
  mmceq  = 0;
  ttcc   = 0;
  minex  = 0;
  meqn   = 0;
}


SAM::~SAM ()
{
  // Deallocate all dynamically allocated arrays
  if (mpmnpc) delete[] mpmnpc;
  if (mmnpc)  delete[] mmnpc;
  if (madof)  delete[] madof;
  if (msc)    delete[] msc;
  if (mpmceq) delete[] mpmceq;
  if (mmceq)  delete[] mmceq;
  if (ttcc)   delete[] ttcc;
  if (minex)  delete[] minex;
  if (meqn)   delete[] meqn;
}


std::pair<int,int> SAM::getNodeAndLocalDof (int idof) const
{
  for (int n = 1; n <= nnod; n++)
    if (madof[n] > idof)
      return std::make_pair(minex ? minex[n-1] : n, idof-madof[n-1]+1);

  return std::make_pair(0,0);
}


//! \brief Global stream operator printing a \a (node,localDOF) pair.

std::ostream& operator<< (std::ostream& os, const std::pair<int,int>& p)
{
  return os <<'('<< p.first <<','<< p.second <<')';
}


void SAM::print (std::ostream& os) const
{
  os <<"\n\nSAM::mpar: "<< mpar[0];
  for (int i = 1; i < 30; i++)
    os << (i%10 ? " " : "\n           ") << mpar[i];
  os << std::endl;

  if (mmnpc && mpmnpc)
  {
    os <<"\n\nElement --> Nodes";
    for (int e = 0; e < nel; e++)
    {
      os <<'\n'<< e+1 <<":";
      for (int i = mpmnpc[e]; i < mpmnpc[e+1]; i++)
        if (mmnpc[i-1] > 0)
          os <<' '<< (minex ? minex[mmnpc[i-1]-1] : mmnpc[i-1]);
        else if (mmnpc[i-1] < 0)
          os <<' '<< (minex ? -minex[-mmnpc[i-1]-1] : mmnpc[i-1]);
    }
    os << std::endl;
  }

  if (ttcc && mmceq && mpmceq)
  {
    os <<"\n\nConstraint Equations";
    for (int i = 0; i < nceq; i++)
    {
      int ip = mpmceq[i]-1;
      os <<'\n'<< i+1 <<": "<< this->getNodeAndLocalDof(mmceq[ip]) <<" =";
      if (ttcc[ip] || ip+2 >= mpmceq[i+1]) os <<' '<< ttcc[ip];
      for (int jp = ip+1; jp < mpmceq[i+1]-1; jp++)
      {
	if (ttcc[ip] || jp > ip+1) os <<" +";
	os <<' '<< ttcc[jp] <<'*'<< this->getNodeAndLocalDof(mmceq[jp]);
      }
    }
    os << std::endl;
  }

  if (meqn && madof)
  {
    char code[7];
    os <<"\n\nNode --> Equations";
    for (int n = 0; n < nnod; n++)
    {
      int dof, i = 0;
      os <<'\n'<< (minex ? minex[n] : n+1) <<":";
      for (dof = madof[n]; dof < madof[n+1]; dof++, i++)
      {
	os <<' '<< meqn[dof-1];
	code[i] = msc[dof-1] > 0 ? '0' + msc[dof-1] : '0';
      }
      code[i] = '\0';
      os <<'\t'<< code;
      if ((size_t)n < nodeType.size())
	os <<'\t'<< nodeType[n];
    }
    os << std::endl;
  }
}


bool SAM::initSystemEquations ()
{
  if (!msc && ndof > 0) return false;
  if (!mpmceq && nceq > 0) return false;

  // Initialize the DOF-to-equation connectivity array
  int i, j, ierr = 0;
  meqn = new int[ndof];
  memset(meqn,0,ndof*sizeof(int));
#ifdef USE_F77SAM
  syseq_(msc,mpmceq,mmceq,6,mpar,meqn,ierr);
  for (i = j = 0; i < ndof; i++)
    if (msc[i] == 0) msc[i] = --j; // reaction force indices
#else
  int ndof1  = 0;
  int ndof2  = 0;
  int npdof  = 0;
  int nddof  = 0;
  for (i = 0; i < ndof; i++)
    if (msc[i] == 0)
      msc[i] = -(++nspdof); // reaction force indices
    else if (msc[i] == 1)
      ndof1++;
    else if (msc[i] == 2)
      ndof2++;
    else
      ierr--;

  if (ierr < 0)
    std::cerr <<"SAM::initSystemEquations: The MSC array has "<< -ierr
	      <<" invalid values, it should only contain 0, 1 or 2"<< std::endl;

  for (i = 1; i <= nceq; i++)
  {
    int ip = mpmceq[i-1];
    int jp = mpmceq[i]-1;
    int idof = mmceq[ip-1];
    if (idof < 1 || idof > ndof)
    {
      ierr--;
      std::cerr <<"SAM::initSystemEquations: idof = "<< idof
		<<" is out of range [1,"<< ndof <<"]"<< std::endl;
      continue;
    }
    else if (msc[idof-1] > 0)
    {
      ierr--;
      std::cerr <<"SAM::initSystemEquations: Invalid status code "<< msc[idof-1]
		<<" for slave dof "<< idof << std::endl;
    }
    else if (jp == ip)
      npdof++; // prescribed DOF
    else if (jp > ip)
    {
      nddof++; // slave DOF
      while (ip++ < jp)
      {
	int jdof = mmceq[ip-1];
	if (jdof < 1 || jdof > ndof)
	  ierr--;
	else if (msc[jdof-1] < 1)
	  ierr--;
      }
    }
    else
    {
      ierr--;
      break;
    }
    meqn[idof-1] = -i;
  }

  mpar[3] = ndof1;
  mpar[4] = ndof2;
  mpar[7] = nspdof - nceq;
  mpar[8] = npdof;
  mpar[9] = nddof;
  neq = ndof1 + ndof2;

  i = 1;
  j = ndof1+1;
  for (int idof = 0; idof < ndof; idof++)
    if (msc[idof] == 1)
      meqn[idof] = i++;
    else if (msc[idof] == 2)
      meqn[idof] = j++;
#endif

  if (ierr == 0) return true;

  std::cerr <<"SAM::initSystemEquations: Failure "<< ierr << std::endl;
  return false;
}


int SAM::getNoNodes (char dofType) const
{
  if (dofType == 'A')
    return nnod;

  int n = 0;
  for (int i = 0; i < nnod; i++)
    if (nodeType[i] == dofType) n++;

  return n;
}


int SAM::getMaxDofCouplings () const
{
  std::vector<IntSet> dofc;
  if (!this->getDofCouplings(dofc))
    return 0;

  size_t maxdofc = 0;
  for (int i = 0; i < neq; i++)
    if (dofc[i].size() > maxdofc)
      maxdofc = dofc[i].size();

  return maxdofc;
}


bool SAM::getNoDofCouplings (IntVec& nnz) const
{
  nnz.clear();

  // Find the set of DOF couplings for each free DOF
  std::vector<IntSet> dofc;
  if (!this->getDofCouplings(dofc))
    return false;

  // Find total number of DOF couplings or non-zeroes in the system matrix
  nnz.reserve(neq);
  for (int i = 0; i < neq; i++)
    nnz.push_back(dofc[i].size());

  return true;
}


bool SAM::getDofCouplings (IntVec& irow, IntVec& jcol) const
{
  irow.clear();
  jcol.clear();

  // Find the set of DOF couplings for each free DOF
  std::vector<IntSet> dofc;
  if (!this->getDofCouplings(dofc))
    return false;

  irow.resize(neq+1);
  irow.front() = 0;

  // Find total number of dof couplings or non-zeroes in the system matrix
  int i;
  for (i = 0; i < neq; i++)
    irow[i+1] = irow[i] + dofc[i].size();

  jcol.reserve(irow.back());
  for (i = 0; i < neq; i++)
    jcol.insert(jcol.end(),dofc[i].begin(),dofc[i].end());

  return true;
}


bool SAM::getDofCouplings (std::vector<IntSet>& dofc) const
{
  dofc.resize(neq);
  for (int iel = 1; iel <= nel; iel++)
  {
    IntVec meen;
    if (!this->getElmEqns(meen,iel))
      return false;

    size_t i, j;
    int ieq, ip, ipmceq1, ipmceq2, jeq, jp, jpmceq1, jpmceq2;
    for (j = 0; j < meen.size(); j++)
      if ((jeq = meen[j]) > 0)
      {
        dofc[jeq-1].insert(jeq);
        for (i = 0; i < j; i++)
          if ((ieq = meen[i]) > 0)
          {
            dofc[ieq-1].insert(jeq);
            dofc[jeq-1].insert(ieq);
          }
      }
      else if (jeq < 0)
      {
        jpmceq1 = mpmceq[-jeq-1];
        jpmceq2 = mpmceq[-jeq]-1;
        for (jp = jpmceq1; jp < jpmceq2; jp++)
          if (mmceq[jp] > 0)
          {
            jeq = meqn[mmceq[jp]-1];
            for (i = 0; i < meen.size(); i++)
              if ((ieq = meen[i]) > 0)
              {
                dofc[ieq-1].insert(jeq);
                dofc[jeq-1].insert(ieq);
              }
              else if (ieq < 0)
              {
                ipmceq1 = mpmceq[-ieq-1];
                ipmceq2 = mpmceq[-ieq]-1;
                for (ip = ipmceq1; ip < ipmceq2; ip++)
                  if (mmceq[ip] > 0)
                  {
                    ieq = meqn[mmceq[ip]-1];
                    dofc[ieq-1].insert(jeq);
                  }
              }
          }
      }
  }

  return true;
}


bool SAM::initForAssembly (SystemMatrix& sysK, SystemVector& sysRHS,
			   Vector* reactionForces, bool dontLockSP) const
{
  sysK.initAssembly(*this,dontLockSP);
  return this->initForAssembly(sysRHS,reactionForces);
}


bool SAM::initForAssembly (SystemMatrix& sysM) const
{
  sysM.initAssembly(*this);
  return neq > 0 ? true : false;
}


bool SAM::initForAssembly (SystemVector& sysRHS, Vector* reactionForces) const
{
  sysRHS.redim(neq);
  sysRHS.init();
  if (reactionForces)
    reactionForces->resize(nspdof,true);

  return neq > 0 ? true : false;
}


bool SAM::assembleSystem (SystemMatrix& sysK, SystemVector& sysRHS,
			  const Matrix& eK, int iel,
			  Vector* reactionForces) const
{
  if (reactionForces)
  {
    Vector eS;
    IntVec meen;
    if (!this->getElmEqns(meen,iel,eK.rows()))
      return false;

    // Add (appropriately weighted) terms corresponding to constrained
    // (dependent and prescribed) DOFs in eK into reactionForces
    for (size_t j = 1; j <= meen.size(); j++)
    {
      int jceq = -meen[j-1];
      if (jceq < 1) continue;

      int jp = mpmceq[jceq-1];
      Real c0 = ttcc[jp-1];

      for (size_t i = 1; i <= meen.size(); i++)
	if (meen[i-1] <= 0)
	{
	  eS.resize(eK.rows());
	  eS(i) = -c0*eK(i,j);
	}
    }

    if (!eS.empty())
      this->assembleReactions(*reactionForces,eS,iel);
  }

  return sysK.assemble(eK,*this,sysRHS,iel);
}


bool SAM::assembleSystem (SystemMatrix& sysM,
			  const Matrix& eM, int iel) const
{
  return sysM.assemble(eM,*this,iel);
}


bool SAM::assembleSystem (SystemVector& sysRHS,
			  const Matrix& eK, int iel,
			  Vector* reactionForces) const
{
  IntVec meen;
  if (!this->getElmEqns(meen,iel,eK.rows()))
    return false;

  // Add (appropriately weighted) terms corresponding to constrained
  // (dependent and prescribed) DOFs in eK into sysRHS and reactionForces
  Vector eS;
  Real* sysrhsPtr = sysRHS.getPtr();
  for (size_t j = 1; j <= meen.size(); j++)
  {
    int jceq = -meen[j-1];
    if (jceq < 1) continue;

    int jp = mpmceq[jceq-1];
    Real c0 = ttcc[jp-1];

    for (size_t i = 1; i <= meen.size(); i++)
    {
      int ieq = meen[i-1];
      this->assembleRHS(sysrhsPtr,-c0*eK(i,j),ieq);
      if (reactionForces && ieq <= 0)
      {
        eS.resize(eK.rows());
        eS(i) = -c0*eK(i,j);
      }
    }
  }

  sysRHS.restore(sysrhsPtr);
  if (!eS.empty())
    this->assembleReactions(*reactionForces,eS,iel);

  return true;
}


bool SAM::assembleSystem (SystemVector& sysRHS,
			  const RealArray& eS, int iel,
			  Vector* reactionForces) const
{
  Real* sysrhsPtr = sysRHS.getPtr();
  int ierr = 0;
#ifdef USE_F77SAM
  int* work = new int[eS.size()];
  addev2_(&eS.front(), ttcc, mpar, madof, meqn, mpmnpc, mmnpc, mpmceq, mmceq,
	  iel, eS.size(), 6, sysrhsPtr, work, ierr);
  delete[] work;
#else
  IntVec meen;
  if (!this->getElmEqns(meen,iel,eS.size()))
    ierr = 1;
  else for (size_t i = 0; i < meen.size(); i++)
    this->assembleRHS(sysrhsPtr,eS[i],meen[i]);
#endif

  sysRHS.restore(sysrhsPtr);
  if (reactionForces)
    this->assembleReactions(*reactionForces,eS,iel);

  return ierr == 0;
}


bool SAM::assembleSystem (SystemVector& sysRHS, const Real* nS, int inod,
			  Vector* reactionForces) const
{
  IntVec mnen;
  if (!this->getNodeEqns(mnen,inod))
    return false;

  Real* sysrhsPtr = sysRHS.getPtr();
  for (size_t i = 0; i < mnen.size(); i++)
    this->assembleRHS(sysrhsPtr,nS[i],mnen[i]);
  sysRHS.restore(sysrhsPtr);

  if (reactionForces)
  {
    int ipR, j, k = 0;
    for (j = madof[inod-1]; j < madof[inod]; j++, k++)
      if ((ipR = -msc[j-1]) > 0 && (size_t)ipR <= reactionForces->size())
        (*reactionForces)(ipR) += nS[k];
  }

  return true;
}


bool SAM::assembleSystem (SystemVector& sysRHS, Real S,
                          const std::pair<int,int>& dof) const
{
  if (dof.first < 1 || dof.first > nnod)
  {
    std::cerr <<"SAM::assembleSystem: Node "<< dof.first
              <<" is out of range [1,"<< nnod <<"]"<< std::endl;
    return false;
  }

  int idof = madof[dof.first-1] + dof.second-1;
  if (dof.second < 1 || idof >= madof[dof.first])
  {
    int nndof = madof[dof.first] - madof[dof.first-1];
    std::cerr <<"SAM::assembleSystem: Local dof "<< dof.second
              <<" is out of range [1,"<< nndof <<"]"<< std::endl;
    return false;
  }

  Real* sysrhsPtr = sysRHS.getPtr();
  this->assembleRHS(sysrhsPtr,S,meqn[idof-1]);
  sysRHS.restore(sysrhsPtr);

  return true;
}


void SAM::addToRHS (SystemVector& sysRHS, const RealArray& S) const
{
  if (ndof < 1 || S.empty()) return;

  Real* sysrhsPtr = sysRHS.getPtr();

  for (size_t i = 0; i < S.size() && i < (size_t)ndof; i++)
    this->assembleRHS(sysrhsPtr,S[i],meqn[i]);

  sysRHS.restore(sysrhsPtr);
}


void SAM::assembleRHS (Real* RHS, Real value, int ieq) const
{
  int iceq = -ieq;
  if (ieq > 0)
    RHS[ieq-1] += value;
  else if (iceq > 0)
    for (int ip = mpmceq[iceq-1]; ip < mpmceq[iceq]-1; ip++)
      if (mmceq[ip] > 0)
      {
        ieq = meqn[mmceq[ip]-1];
        RHS[ieq-1] += ttcc[ip]*value;
      }
}


void SAM::assembleReactions (Vector& reac, const RealArray& eS, int iel) const
{
  int i, j, k, ipR, node;
  int ip = mpmnpc[iel-1];
  int nenod = mpmnpc[iel] - ip;
  for (i = k = 0; i < nenod; i++, ip++)
    if ((node = mmnpc[ip-1]) < 0)
      k += madof[-node] - madof[-node-1];
    else if (node > 0)
      for (j = madof[node-1]; j < madof[node]; j++, k++)
        if ((ipR = -msc[j-1]) > 0 && (size_t)ipR <= reac.size())
          reac(ipR) += eS[k];
}


bool SAM::getElmNodes (IntVec& mnpc, int iel) const
{
  mnpc.clear();
  if (iel < 1 || iel > nel)
  {
    std::cerr <<"SAM::getElmNodes: Element "<< iel <<" is out of range [1,"
	      << nel <<"]"<< std::endl;
    return false;
  }

  int ip = mpmnpc[iel-1]-1;
  int jp = mpmnpc[iel]-1;
  if (jp <= ip) return true;

  mnpc.reserve(jp-ip);
  mnpc.insert(mnpc.end(),mmnpc+ip,mmnpc+jp);
  return true;
}


bool SAM::getElmEqns (IntVec& meen, int iel, int nedof) const
{
  meen.clear();
  if (iel < 1 || iel > nel)
  {
    std::cerr <<"SAM::getElmEqns: Element "<< iel <<" is out of range [1,"
	      << nel <<"]"<< std::endl;
    return false;
  }

  int ip = mpmnpc[iel-1];
  int nenod = mpmnpc[iel] - ip;
  if (nenod <= 0) return true;
#ifndef USE_F77SAM
  int oldof = nedof;
#endif
  if (nedof < 1) nedof = nenod*ndof/nnod;

#ifdef USE_F77SAM
  int neldof, neslv, neprd;
  meen.resize(nedof,0);
  elmeq_(madof,mmnpc+ip-1,mpmceq,meqn,nenod,&meen.front(),neldof,neslv,neprd);
  if (neldof == nedof) return true;
#else
  meen.reserve(nedof);
  for (int i = 0; i < nenod; i++, ip++)
  {
    int node = mmnpc[ip-1];
    if (node > 0)
      meen.insert(meen.end(),meqn+madof[node-1]-1,meqn+madof[node]-1);
    else if (node < 0)
      meen.insert(meen.end(),madof[-node]-madof[-node-1],0);
  }
  int neldof = meen.size();
  if (neldof == nedof || oldof < 1) return true;
#endif

  std::cerr <<"SAM::getElmEqns: Invalid element matrix dimension "
	    << nedof <<" (should have been "<< neldof <<")"<< std::endl;
  return false;
}


size_t SAM::getNoElmEqns (int iel) const
{
  size_t result = 0;
  if (iel < 1 || iel > nel)
    std::cerr <<"SAM::getNoElmEqns: Element "<< iel <<" is out of range [1,"
	      << nel <<"]"<< std::endl;

  else for (int ip = mpmnpc[iel-1]; ip < mpmnpc[iel]; ip++)
  {
    int node = abs(mmnpc[ip-1]);
    if (node > 0)
      result += madof[node] - madof[node-1];
  }

  return result;
}


bool SAM::getNodeEqns (IntVec& mnen, int inod) const
{
  mnen.clear();
  if (inod < 1 || inod > nnod)
  {
    std::cerr <<"SAM::getNodeEqns: Node "<< inod <<" is out of range [1,"
	      << nnod <<"]"<< std::endl;
    return false;
  }

  mnen.reserve(madof[inod]-madof[inod-1]);
  mnen.insert(mnen.end(),meqn+madof[inod-1]-1,meqn+madof[inod]-1);

  return true;
}


std::pair<int,int> SAM::getNodeDOFs (int inod) const
{
  if (inod < 1 || inod > nnod) return std::make_pair(0,0);

  return std::make_pair(madof[inod-1],madof[inod]-1);
}


int SAM::getEquation (int inod, int ldof) const
{
  if (inod < 1 || inod > nnod || ldof < 1) return -1;

  int idof = madof[inod-1] + ldof-1;
  if (idof >= madof[inod]) return -2;

  int ieq = meqn[idof-1];
  return ieq > 0 ? ieq : 0;
}


bool SAM::expandSolution (const SystemVector& solVec, Vector& dofVec,
			  Real scaleSD) const
{
  if (solVec.dim() < (size_t)neq) return false;

  return this->expandVector(solVec.getRef(),dofVec,scaleSD);
}


bool SAM::expandVector (const Vector& solVec, Vector& dofVec) const
{
  if (solVec.size() < (size_t)neq) return false;

  return this->expandVector(solVec.ptr(),dofVec,0.0);
}


bool SAM::expandVector (const Real* solVec, Vector& dofVec, Real scaleSD) const
{
  if (!meqn) return false;

  dofVec.resize(ndof,true);
#ifdef USE_F77SAM
  expand_(solVec, ttcc, mpmceq, mmceq, meqn, Real(1), scaleSD, ndof, neq,
	  dofVec.ptr());
#else
  for (int idof = 0; idof < ndof; idof++)
  {
    int ieq = meqn[idof];
    int iceq = -ieq;
    if (ieq > 0)
      dofVec[idof] += solVec[ieq-1];
    else if (iceq > 0)
    {
      int ip = mpmceq[iceq-1];
      dofVec[idof] += scaleSD*ttcc[ip-1];
      for (; ip < mpmceq[iceq]-1; ip++)
	if (mmceq[ip] > 0)
	{
	  ieq = meqn[mmceq[ip]-1];
	  if (ieq > 0 && ieq <= neq)
	    dofVec[idof] += ttcc[ip]*solVec[ieq-1];
	}
    }
  }
#endif

  return true;
}


bool SAM::applyDirichlet (Vector& dofVec) const
{
  if (!meqn) return false;

  for (int idof = 0; idof < ndof; idof++)
  {
    int iceq = -meqn[idof];
    if (iceq > 0)
    {
      int ip = mpmceq[iceq-1];
      dofVec[idof] = ttcc[ip-1];
    } else if (iceq == 0)
      dofVec[idof] = 0.0;
  }

  return true;
}



Real SAM::dot (const Vector& x, const Vector& y, char dofType) const
{
  if (nodeType.empty() || dofType == 'A')
    return x.dot(y); // All nodes are of the same type, or consider all of them

  // Consider only the dofType nodes
  int i, j, n = x.size() < y.size() ? x.size() : y.size();
  Real retVal = Real(0);
  for (i = 0; i < nnod; i++)
    if (nodeType[i] == dofType)
      for (j = madof[i]-1; j < madof[i+1]-1 && j < n; j++)
	retVal += x[j]*y[j];

  return retVal;
}


Real SAM::normL2 (const Vector& x, char dofType) const
{
  if (x.empty())
    return Real(0);
  else if (nodeType.empty()) // All nodes are of the same type
    return x.norm2()/sqrt(x.size());

  // Consider only the dofType nodes
  int count, i, j, n = x.size();
  Real retVal = Real(0);
  for (count = i = 0; i < nnod; i++)
    if (nodeType[i] == dofType)
      for (j = madof[i]-1; j < madof[i+1]-1 && j < n; j++)
      {
	retVal += x[j]*x[j];
	count ++;
      }

  return count > 0 ? sqrt(retVal/count) : retVal;
}


Real SAM::normInf (const Vector& x, size_t& comp, char dofType) const
{
  if (x.empty() || comp < 1)
    return Real(0);
  else if (nodeType.empty())
  {
    // All nodes are of the same type with the same number of nodal DOFs
    int nndof = madof[1] - madof[0];
    if ((int)comp <= nndof)
      return x.normInf(--comp,nndof);
    else
      return Real(0);
  }

  // Consider only the dofType nodes
  int i, k = comp-2;
  Real retVal = Real(0);
  for (i = 0; i < nnod; i++)
    if (nodeType[i] == dofType)
    {
      int idof = madof[i] + k;
      if (idof >= 0 && idof < madof[i+1]-1)
	if (fabs(x[idof]) > retVal)
	{
	  comp = i+1;
	  retVal = fabs(x[idof]);
	}
    }

  return retVal;
}


Real SAM::normReact (const Vector& u, const Vector& rf) const
{
  Real retVal = Real(0);

  for (int i = 0; i < ndof; i++)
    if (msc[i] < 0 && -msc[i] <= (int)rf.size())
      retVal += u[i]*rf(-msc[i]);

  return 0.5*retVal;
}


Real SAM::getReaction (int dir, const Vector& rf) const
{
  Real retVal = Real(0);

  if (dir > 0)
    for (int i = 0; i < nnod; i++)
    {
      int idof = madof[i]+dir-2;
      if (idof < madof[i+1]-1 && msc[idof] < 0 && -msc[idof] <= (int)rf.size())
	retVal += rf(-msc[idof]);
    }

  return retVal;
}


bool SAM::getNodalReactions (int inod, const Vector& rf, Vector& nrf) const
{
  if (inod < 1 || inod > nnod)
  {
    std::cerr <<"SAM::getNodalReactions: Node "<< inod <<" is out of range [1,"
	      << nnod <<"]"<< std::endl;
    return false;
  }

  bool haveRF = false;
  int ip = madof[inod-1]-1;
  nrf.resize(madof[inod]-1-ip,true);
  for (size_t i = 0; i < nrf.size(); i++, ip++)
    if (msc[ip] < 0 && -msc[ip] <= (int)rf.size())
    {
      haveRF = true;
      nrf[i] = rf(-msc[ip]);
    }

  return haveRF;
}


IntSet SAM::getEquations(char dofType, int dof) const
{
  IntSet result;
  for (int inod = 0; inod < nnod; inod++) {
    char type = this->getNodeType(inod+1);
    if (type == dofType || (dofType == 'D' && type == ' ')) {
      if (dof > 0) {
	int idof = madof[inod] + dof-1;
	int ieq = idof < madof[inod+1] ? meqn[idof-1] : 0;
	if (ieq > 0) result.insert(ieq);
      }
      else for (int idof = madof[inod]; idof < madof[inod+1]; idof++) {
	int ieq = meqn[idof-1];
	if (ieq > 0) result.insert(ieq);
      }
    }
  }

  return result;
}
