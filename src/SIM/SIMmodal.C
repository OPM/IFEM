// $Id$
//==============================================================================
//!
//! \file SIMmodal.C
//!
//! \date Aug 30 2019
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of modal FEM system.
//!
//==============================================================================

#include "SIMmodal.h"
#include "AlgEqSystem.h"
#include "NewmarkMats.h"
#include "SAM.h"
#include <numeric>
#include <algorithm>
#ifdef HAS_CEREAL
#include <cereal/cereal.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#endif


/*!
  \brief A simple SAM class for diagonal systems.
*/

class SAMmodal : public SAM
{
public:
  //! \brief The constructor initializes the arrays for a diagonal system.
  SAMmodal(int n)
  {
    nmmnpc = nel = nnod = ndof = neq = n;
    mmnpc  = new int[n];
    mpmnpc = new int[n+1];
    madof  = new int[n+1];
    msc    = new int[n];
    std::iota(mpmnpc,mpmnpc+n+1,1);
    std::iota(mmnpc ,mmnpc +n  ,1);
    std::iota(madof ,madof +n+1,1);
    std::fill(msc   ,msc   +n  ,1);
    this->initSystemEquations();
  }

  //! \brief Empty destructor.
  virtual ~SAMmodal() {}
};


SIMmodal::SIMmodal (std::vector<Mode>& modes) : myModes(modes)
{
  modalSys = nullptr;
  modalSam = nullptr;
  myElmMat = nullptr;
}


SIMmodal::~SIMmodal ()
{
  delete modalSys;
  delete modalSam;
  delete myElmMat;
}


bool SIMmodal::swapSystem (AlgEqSystem*& sys, SAM*& sam)
{
  if (!modalSys || !sam)
    return false;

  std::swap(modalSys,sys);
  std::swap(modalSam,sam);
  return true;
}


bool SIMmodal::expandSolution (const Vectors& mSol, Vectors& pSol) const
{
  pSol.resize(mSol.size());
  for (size_t i = 0; i < pSol.size(); i++)
    if (mSol[i].size() != myModes.size())
    {
      std::cerr <<" *** SIMmodal::expandSolution: Invalid dimension"
                <<" on modal solution vector "<< i+1 <<": "<< mSol[i].size()
                <<" != "<< myModes.size() << std::endl;
      return false;
    }
    else
    {
      pSol[i].resize(myModes.front().eigVec.size(),true);
      for (size_t j = 0; j < myModes.size(); j++)
        pSol[i].add(myModes[j].eigVec,mSol[i][j]);
    }

  return true;
}


bool SIMmodal::assembleModalSystem (const TimeDomain& time,
                                    const Vectors& pSol, const Vector& Rhs,
                                    double beta, double gamma)
{
  if (!modalSam)
    // Create a diagonal SAM object
    modalSam = new SAMmodal(myModes.size());

  if (!modalSys)
  {
    // Create a diagonal equation system
    modalSys = new AlgEqSystem(*modalSam);
    if (!modalSys->init(LinAlg::DIAG))
      return false;
  }

  if (!myElmMat)
  {
    // Create a single-dof element matrix object,
    // the "elements" now being the eigenmodes
    myElmMat = new NewmarkMats(0.0,0.0,beta,gamma);
    myElmMat->resize(3,1);
    myElmMat->redim(1);
  }

  myElmMat->setStepSize(time.dt,time.it);
  myElmMat->vec.resize(pSol.size(),Vector(1));

  int iu = pSol.size() - 3; // index to modal displacement (u)
  if (iu < 0)
  {
    std::cerr <<" SIMmodal::assembleModalSystem: No solutions."<< std::endl;
    return false;
  }

#if SP_DEBUG > 1
  if (time.first && time.it == 0)
    for (const Mode& m : myModes)
      std::cout <<"\nMode shape "<< m.eigNo
                <<": omega = "<< 2.0*M_PI*m.eigVal << m.eigVec;
#endif

  // Assemble the modal system
  bool ok = true;
  modalSys->initialize(true);
  for (size_t m = 0; m < myModes.size() && ok; m++)
  {
    double omega = 2.0*M_PI*myModes[m].eigVal; // Angular eigenfrequency

#if SP_DEBUG > 2
    std::cout <<"\nProcessing mode shape "<< myModes[m].eigNo
              <<": omega = "<< omega << std::endl;
#endif

    // Extract current solution associated with the m'th mode shape
    for (size_t i = 0; i < pSol.size(); i++)
      myElmMat->vec[i].fill(pSol[i][m]);

    myElmMat->A[1].fill(1.0);                        // Modal mass
    myElmMat->A[2].fill(omega*omega);                // Modal stiffness
    myElmMat->b[0].fill(myModes[m].eigVec.dot(Rhs)); // Modal load
    myElmMat->b[0].add({-omega*omega*pSol[iu][m]});  // Modal residual
    ok = modalSys->assemble(myElmMat,myModes[m].eigNo);
  }
  ok &= modalSys->finalize(true);

  return ok;
}


#ifdef HAS_CEREAL
//! \brief Serializes Mode data to/from the \a archive.
template<class T> void doSerialize (T& archive, Mode& mode)
{
  archive(mode.eigNo);
  archive(mode.eigVal);
  archive(mode.eigVec);
}
#endif


bool SIMmodal::saveModes (std::map<std::string,std::string>& data) const
{
#ifdef HAS_CEREAL
  std::ostringstream str;
  {
    cereal::BinaryOutputArchive archive(str);
    archive(myModes.size());
    for (const Mode& mode : myModes)
      doSerialize(archive,const_cast<Mode&>(mode));
  }
  data["ModalSIM"] = str.str();
  return true;
#else
  return false;
#endif
}


bool SIMmodal::restoreModes (const std::map<std::string,std::string>& data)
{
#ifdef HAS_CEREAL
  std::map<std::string,std::string>::const_iterator sit = data.find("ModalSIM");
  if (sit != data.end())
  {
    std::stringstream str(sit->second);
    cereal::BinaryInputArchive archive(str);
    size_t size;
    archive(size);
    myModes.resize(size);
    for (Mode& mode : myModes)
      doSerialize(archive,mode);
    return true;
  }
#endif
  return false;
}
