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
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml2.h"
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
  parsed = false;
  alpha1 = alpha2 = 0.0;

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


bool SIMmodal::parseParams (const tinyxml2::XMLElement* elem)
{
  if (parsed)
    IFEM::cout <<"\t(skipped)"<< std::endl;
  else if (elem && !strcasecmp(elem->Value(),"newmarksolver"))
  {
    utl::getAttribute(elem,"alpha1",alpha1);
    utl::getAttribute(elem,"alpha2",alpha2);
  }
  else
    return false;

  return true;
}


bool SIMmodal::swapSystem (AlgEqSystem*& sys, SAM*& sam)
{
  if (!modalSys || !sam)
    return false;

  std::swap(modalSys,sys);
  std::swap(modalSam,sam);
  return true;
}


const Vectors& SIMmodal::expandSolution (const Vectors& mSol)
{
  if (sol.empty())
    sol.resize(mSol.size());
  else if (sol.size() != mSol.size())
  {
    std::cerr <<" *** SIMmodal::expandSolution: Invalid modal solution, "
              <<" size(mSol) = "<< mSol.size() <<" ("<< sol.size()
              <<" expected)."<< std::endl;
    sol.clear();
  }

  for (size_t i = 0; i < sol.size(); i++)
    if (mSol[i].size() != myModes.size())
    {
      std::cerr <<" *** SIMmodal::expandSolution: Invalid dimension"
                <<" on modal solution vector "<< i+1 <<": "<< mSol[i].size()
                <<" != "<< myModes.size() << std::endl;
      sol.clear();
    }
    else
    {
      if (sol[i].empty())
        sol[i].resize(myModes.front().eigVec.size());
      else if (sol[i].size() == myModes.front().eigVec.size())
        sol[i].fill(0.0);
      else
      {
        std::cerr <<" *** SIMmodal::expandSolution: Logic error, "
                  <<" size(sol["<< i <<"]) = "<< sol[i].size()
                  <<" != size(eigVec) = "<< myModes.front().eigVec.size()
                  << std::endl;
        sol.clear();
        break;
      }

      for (size_t j = 0; j < myModes.size(); j++)
        sol[i].add(myModes[j].eigVec,mSol[i][j]);
    }

  return sol;
}


const Vector& SIMmodal::expandedSolution (int idx) const
{
  if (idx >= 0 && idx < static_cast<int>(sol.size()))
    return sol[idx];

  static Vector empty;
  return empty;
}


/*!
  This method assembles the diagonal modal equation system of the
  dynamic problem, integrated by the Newmark HHT-method.
  If \a beta is set to zero, a quasi-static solution is calculated instead,
  in which the modal mass and damping is ignored and only the modal stiffness
  (i.e., the angular eigenfrequencies squared) enter the equation system.
*/

bool SIMmodal::assembleModalSystem (const TimeDomain& time, const Vectors& mSol,
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
    if (beta > 0.0)
    {
      myElmMat = new NewmarkMats(0.0,0.0,beta,gamma);
      myElmMat->resize(alpha1 > 0.0 || alpha2 > 0.0 ? 4 : 3, 1);
    }
    else // quasi-static simulation
    {
      myElmMat = new ElmMats();
      myElmMat->resize(1,1);
    }
    myElmMat->redim(1);
  }

  myElmMat->vec.resize(mSol.size(),Vector(1));
  myElmMat->setStepSize(time.dt,time.it);

  int iu = mSol.size() - (beta > 0.0 ? 3 : 1); // index to modal displacement
  if (iu < 0)
  {
    std::cerr <<" *** SIMmodal::assembleModalSystem: No solutions."<< std::endl;
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
  double omega, damping;
  modalSys->initialize(true);
  for (size_t m = 0; m < myModes.size() && ok; m++)
  {
    omega = 2.0*M_PI*myModes[m].eigVal; // Angular eigenfrequency
    if (myModes[m].damping > 0.0)
      damping = myModes[m].damping;
    else // Pure Rayleigh damping
      damping = alpha1 + alpha2*omega*omega;

#if SP_DEBUG > 2
    std::cout <<"\nProcessing mode shape "<< myModes[m].eigNo
              <<": omega = "<< omega
              <<", damping = "<< damping << std::endl;
#endif

    // Extract current solution associated with the m'th mode shape
    for (size_t i = 0; i < mSol.size(); i++)
      myElmMat->vec[i].fill(mSol[i][m]);

    if (beta > 0.0)
    {
      myElmMat->A[1].fill(1.0);         // Modal mass
      myElmMat->A[2].fill(omega*omega); // Modal stiffness
      if (myElmMat->A.size() > 3)
        myElmMat->A[3].fill(damping);   // Modal damping
    }
    else // quasi-static simulation
      myElmMat->A[0].fill(omega*omega); // Modal stiffness

    myElmMat->b[0].fill(myModes[m].eigVec.dot(Rhs));  // Modal load
    if (beta > 0.0)
      myElmMat->b[0].add({-omega*omega*mSol[iu][m]}); // Modal residual

    ok = modalSys->assemble(myElmMat,myModes[m].eigNo);
  }
  ok &= modalSys->finalize(true);

  return ok;
}


bool SIMmodal::saveModes (std::map<std::string,std::string>& data) const
{
#ifdef HAS_CEREAL
  std::ostringstream str;
  {
    cereal::BinaryOutputArchive archive(str);
    archive(myModes.size());
    for (const Mode& mode : myModes)
    {
      archive(mode.eigNo);
      archive(mode.eigVal);
      archive(mode.eigVec);
    }
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
    {
      archive(mode.eigNo);
      archive(mode.eigVal);
      archive(mode.eigVec);
    }
    return true;
  }
#endif
  return false;
}
