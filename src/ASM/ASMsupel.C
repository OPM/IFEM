// $Id$
//==============================================================================
//!
//! \file ASMsupel.C
//!
//! \date Mar 30 2021
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of general superelements.
//!
//==============================================================================

#include "ASMsupel.h"
#include "Integrand.h"
#include "GlobalIntegral.h"
#include "ElementBlock.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include <numeric>

#ifdef HAS_FMXREADER
extern "C" {
  void initfmx_();
  int  readfmx_(const char* fname, const int* itype,
                double* data, const int* nval, const int nchar);
  int  readfsm_(const char* fname,
                int* ndof, int* ndof1, int* ndof2, int* nceq, int* nmmceq,
                int* meqn, int* meqn1, int* meqn2,
                int* mmceq, int* mpmceq, double* ttcc,
                const int nchar = 0);
}
#endif


bool ASMsupel::read (std::istream& is)
{
  // Lambda function for reading the supernode coordinates.
  auto&& readCoord = [&is](Vec3Vec& Xsup)
  {
    int numNod = 0;
    is >> numNod;
    Xsup.resize(numNod);
    for (Vec3& Xn : Xsup) is >> Xn;
  };

  // Lambda function for reading a FEDEM superelement matrix file.
  auto&& readFMX = [&is](Matrix& M, int itype) -> int
  {
    int m, n, ierr = -99;
    std::string fileName;
    is >> m >> n >> fileName;
    M.resize(m,n);
#ifdef HAS_FMXREADER
    initfmx_();
    int ndim = m*n;
    ierr = readfmx_(fileName.c_str(),&itype,M.ptr(),&ndim,fileName.size());
#endif
    if (ierr < 0)
      std::cerr <<" *** readFMX: Failed to read matrix file \""
                << fileName <<"\" of type "<< itype
                <<" (ierr = "<< ierr <<")" << std::endl;

    return ierr;
  };

  // Lambda function for reading a SAM data file from FEDEM.
  auto&& readFSM = [&is](MiniSAM& sam) -> int
  {
    int ierr = -99;
    std::string fileName;
    is >> fileName;
#ifdef HAS_FMXREADER
    initfmx_();
    double dummy = 0.0;
    int ndof, ndof1, ndof2, nceq, nmmceq;
    ierr = readfsm_(fileName.c_str(),&ndof,&ndof1,&ndof2,&nceq,&nmmceq,
                    &ierr,&ierr,&ierr,&ierr,&ierr,&dummy,fileName.size());
    if (ierr >= 0)
    {
      sam.meqn.resize(ndof);
      sam.meqn1.resize(ndof1);
      sam.meqn2.resize(ndof2);
      sam.ttcc.resize(nmmceq);
      sam.mmceq.resize(nmmceq);
      sam.mpmceq.resize(nceq+1);
      ierr = readfsm_("",&ndof,&ndof1,&ndof2,&nceq,&nmmceq,
                      sam.meqn.data(),sam.meqn1.data(),sam.meqn2.data(),
                      sam.mmceq.data(),sam.mpmceq.data(),sam.ttcc.data());
    }
#endif
    if (ierr < 0)
      std::cerr <<" *** readFSM: Failed to read SAM arrays from file \""
                << fileName <<"\" (ierr = "<< ierr <<")" << std::endl;
    else
      sam.findDofP2();

    return ierr;
  };

  // Lambda function checking if myElmMat needs to be allocated,
  // and whether we are reading a binary file or not.
  auto&& checkMatrix = [&is,this](bool fmx = false) -> bool
  {
    if (myElmMat.empty())
      myElmMat.resize(2,1,3);
    if (fmx)
    {
      char c = 0;
      if (is.get(c) && c == 'x')
        return true;
      else if (is)
        is.putback(c);
    }
    return false;
  };

  char c = 0;
  Matrix tmpMat;
  std::string comment;
  int readMat = 0, ierr = 0;
  while (readMat < 7 && is.get(c) && ierr >= 0)
    switch (c) {
    case 'K':
    case 'k':
      // Read stiffness matrix file
      if (checkMatrix(true))
        ierr = readFMX(myElmMat.A.front(),1);
      else
      {
        is >> myElmMat.A.front();
        if (c == 'K') // assume stored column-wise
          myElmMat.A.front().transpose();
      }
      readMat |= 1;
      break;

    case 'M':
    case 'm':
      // Read mass matrix file
      if (checkMatrix(true))
        ierr = readFMX(myElmMat.A.back(),2);
      else
      {
        is >> myElmMat.A.back();
        if (c == 'M') // assume stored column-wise
          myElmMat.A.back().transpose();
      }
      readMat |= 1;
      break;

    case 'B':
      // Read recovery matrix file (binary only)
      if (is.get(c) && c == 'x')
        ierr = readFMX(myRecMat,4);
      else if (is)
        is.putback(c);
      break;

    case 'S':
      // Read SAM matrix file (binary only)
      if (is.get(c) && c == 'x')
        ierr = readFSM(mySAM);
      else if (is)
        is.putback(c);
      break;

    case 'g':
      // Read gravity forces file (binary only)
      if (checkMatrix(true))
      {
        if ((ierr = readFMX(tmpMat,3)) >= 0 && !gravity.isZero())
          if (!tmpMat.multiply(gravity.vec(),myElmMat.b.front()))
            return false;
        readMat |= 2;
      }
      break;

    case 'R':
      // Read right-hand-side vector
      checkMatrix();
      is >> myElmMat.b.front();
      readMat |= 2;
      break;

    case 'L':
      // Read load vector.
      // It is assumed stored as a ndof x 1 matrix and not a vector.
      checkMatrix();
      is >> tmpMat;
      myElmMat.b.front() = tmpMat.getColumn(1);
      readMat |= 2;
      break;

    case 'G':
      // Read supernode coordinates
      readCoord(myNodes);
      readMat |= 4;
      break;

    case '#':
      // Comment line - print and ignore
      if (std::getline(is,comment))
        std::cout << comment << std::endl;
      break;

    case '\n':
      // Blank line - skip
      break;

    default:
      is.putback(c);
      if (readMat)
      {
        std::cerr <<" *** ASMsupel::read: Unknown label "<< c
                  <<" ("<< static_cast<int>(c) <<")"<< std::endl;
        return false;
      }
      else
      {
        // Assuming the order G, K, L but without the labels
        checkMatrix();
        readCoord(myNodes);
        is >> myElmMat.A.front() >> myElmMat.b.front();
        readMat = 7;
      }
    }

  // Calculate the total external force
  for (size_t i = 0; i < myElmMat.c.size() && i < 3; i++)
    myElmMat.c[i] = myElmMat.b.front().sum(i,nf);

  return readMat == 7 && is.good() && ierr >= 0;
}


void ASMsupel::MiniSAM::findDofP2 ()
{
  dofP2.clear();
  dofP2.reserve(meqn2.size());

  int ipos = 0;
  for (int ieq : meqn)
    if ((ipos = utl::findIndex(meqn2,ieq)) >= 0)
      dofP2.push_back(ipos+1);
}


void ASMsupel::MiniSAM::print (std::ostream& os) const
{
  if (meqn.empty())
    return;

  if (mpmceq.size() > 1)
  {
    os <<"\nSAM constraints:";
    for (size_t iceq = 1; iceq < mpmceq.size(); iceq++)
    {
      int ip = mpmceq[iceq-1];
      if (ip < mpmceq[iceq])
      {
        os <<"\n"<< iceq <<": "<< mmceq[ip-1];
        if (fabs(ttcc[ip-1]) > 1.0e-15)
          os <<"*("<< ttcc[ip-1] <<")";
        while (++ip < mpmceq[iceq])
          os <<" "<< mmceq[ip-1] <<"*("<< ttcc[ip-1] <<")";
      }
    }
  }

  os <<"\nSAM meqn:";
  for (size_t idof = 0; idof < meqn.size(); idof++)
    os << (idof%6 ? " " : "\n") << meqn[idof];
  os <<"\nSAM meqn1:";
  for (size_t i1 = 0; i1 < meqn1.size(); i1++)
    os << (i1%6 ? " " : "\n") << meqn1[i1];
  os <<"\nSAM meqn2:";
  for (size_t i2 = 0; i2 < meqn2.size(); i2++)
    os << (i2%6 ? " " : "\n") << meqn2[i2];
  os <<"\nSAM dofP2:";
  for (size_t i = 0; i < dofP2.size(); i++)
    os << (i%6 ? " " : "\n") << dofP2[i];
  os << std::endl;
}


bool ASMsupel::write (std::ostream& os, int) const
{
  // Write out the spider as a lagrangian mesh
  os <<"# LAGRANGIAN nodes="<< 1+nnod <<" elements="<< nnod
     <<" type=superelement\n";
  os << this->getCoord(0) <<"\n";
  for (const Vec3& X : myNodes)
    os << X <<"\n";
  for (size_t i = 1; i <= nnod; i++)
    os <<"0 "<< i <<"\n";

  return os.good();
}


bool ASMsupel::generateFEMTopology ()
{
  nnod = myNodes.size();
  nel  = 1;

  myMLGE = { ++gEl };
  myMLGN.resize(nnod,0);
  myMNPC.resize(1);
  myMNPC.front().resize(nnod,0);
  std::iota(myMLGN.begin(),myMLGN.end(),gNod+1);
  std::iota(myMNPC.front().begin(),myMNPC.front().end(),0);
  gNod += nnod;

  if (nodeSets.size() == 1 && nodeSets.front().second.empty())
  {
    // Define the set of supernodes for this patch
    nodeSets.front().second.resize(nnod);
    std::iota(nodeSets.front().second.begin(),nodeSets.front().second.end(),1);
  }

  return true;
}


int ASMsupel::getNodeSetIdx (const std::string& setName) const
{
  int iset = 1;
  for (const ASM::NodeSet& ns : nodeSets)
    if (ns.first == setName)
      return iset;
    else
      ++iset;

  return 0;
}


const IntVec& ASMsupel::getNodeSet (int iset) const
{
  if (iset > 0 && iset <= static_cast<int>(nodeSets.size()))
    return nodeSets[iset-1].second;

  return this->ASMbase::getNodeSet(iset);
}


bool ASMsupel::isInNodeSet (int iset, int inod) const
{
  if (iset < 1 || iset > static_cast<int>(nodeSets.size()))
    return false;

  return utl::findIndex(nodeSets[iset-1].second,inod) >= 0;
}


int ASMsupel::parseNodeSet (const std::string& setName, const char* cset)
{
  int iset = this->getNodeSetIdx(setName)-1;
  if (iset < 0)
  {
    iset = nodeSets.size();
    nodeSets.emplace_back(setName,IntVec());
  }

  utl::parseIntegers(nodeSets[iset].second,cset);

  return 1+iset;
}


void ASMsupel::getBoundaryNodes (int lIndex, IntVec& nodes,
                                 int, int, int, bool local) const
{
  nodes = this->getNodeSet(lIndex);
  if (!local)
    for (int& node : nodes)
      node = this->getNodeID(node);
}


Vec3 ASMsupel::getCoord (size_t inod) const
{
  Vec3 Xn;
  if (inod == 0)
  {
    // Calculate patch center
    for (const Vec3& X : myNodes)
      Xn += X;
    Xn /= nnod;
  }
  else if (inod <= nnod)
    Xn = myNodes[inod-1];

  return Xn;
}


void ASMsupel::getNodalCoordinates (Matrix& X, bool) const
{
  X.resize(3,nnod);
  for (size_t inod = 1; inod <= nnod; inod++)
  {
    X(1,inod) = myNodes[inod-1].x;
    X(2,inod) = myNodes[inod-1].y;
    X(3,inod) = myNodes[inod-1].z;
  }
}


bool ASMsupel::getElementCoordinates (Matrix& X, int iel, bool) const
{
  if (iel != 1)
    return false;

  this->getNodalCoordinates(X);
  return true;
}


bool ASMsupel::integrate (Integrand& integrand, GlobalIntegral& glbInt,
                          const TimeDomain&)
{
  myElmMat.rhsOnly = integrand.getMode() >= SIM::RHS_ONLY;
  return glbInt.assemble(&myElmMat,MLGE.front());
}


/*!
  This method creates a simplified visualization of the superelement
  as lines extending from the centroid to each of the supernodes.
*/

bool ASMsupel::tesselate (ElementBlock& grid, const int*) const
{
  grid.unStructResize(nnod,1+nnod);
  grid.setCoor(0,this->getCoord(0));
  for (size_t i = 1; i <= nnod; i++)
  {
    grid.setCoor(i,myNodes[i-1]);
    grid.setNode(2*i-2,0);
    grid.setNode(2*i-1,i);
    grid.setElmId(i,i);
  }

  return true;
}


bool ASMsupel::evalSolution (Matrix& sField, const Vector& locSol,
                             const int*, int, bool) const
{
  size_t i, j, k, nComp = locSol.size() / nnod;
  sField.resize(nComp,1+nnod,true);
  for (j = k = 1; j <= nnod; j++)
    for (i = 1; i <= nComp; i++)
    {
      sField(i,1+j) = locSol(k++);
      sField(i,1) += sField(i,1+j) / (double)nnod;
    }

  return true;
}


bool ASMsupel::transform (const Matrix& Tlg)
{
  if (Tlg.rows() < 3 || Tlg.cols() < 3)
    return true; // No transformation defined

  // Transform nodal coordinates to global system
  for (Vec3& X : myNodes)
    X = Tlg*X;

#ifdef SP_DEBUG
  std::cout <<"\nGlobal coordinates for superelement "<< idx+1;
  for (const Vec3& X : myNodes)
    std::cout <<"\n" << X;
  std::cout << std::endl;
#endif

  // Transform the element matrices to global system
  for (Matrix& A : myElmMat.A)
    if (!utl::transform(A,Tlg))
      return false;

  // Transform the element force vectors to global system
  for (Vector& b : myElmMat.b)
    if (!utl::transform(b,Tlg))
      return false;

  return true;
}


bool ASMsupel::recoverInternals (const Vector& supSol, Vector& fullSol) const
{
#if SP_DEBUG > 2
  mySAM.print(std::cout);
#endif

  // Reorder supSol to match the expected ordering of the recovery matrix
  Vector ve(mySAM.meqn2.size()), vi;
  for (size_t i = 1; i <= mySAM.dofP2.size(); i++)
    ve(mySAM.dofP2[i-1]) = supSol(i);

  // Calculate the internal DOFs; vi = myRecMat * ve
  if (!myRecMat.multiply(ve,vi))
    return false;

  // Establish the total equation-order solution vector, svec
  int nceq = mySAM.mpmceq.size()-1;
  int ndof = mySAM.meqn.size();
  int neq = *std::max_element(mySAM.meqn.begin(),mySAM.meqn.end());
  Vector sveq(neq);
  for (size_t i1 = 1; i1 <= mySAM.meqn1.size(); i1++)
    sveq(mySAM.meqn1[i1-1]) = vi(i1);
  for (size_t i2 = 1; i2 <= mySAM.meqn2.size(); i2++)
    sveq(mySAM.meqn2[i2-1]) = ve(i2);

#if SP_DEBUG > 2
  std::cout <<"ve:"<< std::scientific << ve;
  std::cout.unsetf(std::ios_base::floatfield);
  std::cout <<"vi:"<< vi <<"sveq:"<< sveq;
#endif

  // Expand to DOF-order while accounting for constraint equations
  fullSol.resize(ndof,true);
  for (int idof = 1; idof <= ndof; idof++)
  {
    int ieq  = mySAM.meqn[idof-1];
    int iceq = -ieq;
    int jdof = 0;
    if (ieq > 0 && ieq <= neq)
      fullSol(idof) = sveq(ieq);
    else if (iceq > 0 && iceq <= nceq)
      for (int ip = mySAM.mpmceq[iceq-1]; ip < mySAM.mpmceq[iceq]-1; ip++)
        if ((jdof = mySAM.mmceq[ip]) > 0 && jdof <= ndof)
          if ((ieq = mySAM.meqn[jdof-1]) > 0 && ieq <= neq)
            fullSol(idof) += mySAM.ttcc[ip]*sveq[ieq-1];
  }

#ifdef SP_DEBUG
  std::cout <<"\nExpanded solution vector for substructure "<< idx+1
            << fullSol;
#endif
  return true;
}
