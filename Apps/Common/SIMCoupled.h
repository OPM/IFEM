//==============================================================================
//!
//! \file SIMCoupled.h
//!
//! \date Oct 12 2012
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Coupled SIM class template
//==============================================================================
#ifndef SIM_COUPLED_H_
#define SIM_COUPLED_H_

  template<class T1, class T2>
class SIMCoupled
{
public:
  //! \brief Constructor 
  //! \param[in] s1 Reference to spline FE model for the first simulator
  //! \param[in] s2 Reference to spline FE model for the second simulator
  SIMCoupled(T1& s1, T2& s2) : S1(s1), S2(s2)
  {
  }

  //! \brief Destructor
  virtual ~SIMCoupled()
  {
  }

  //! \brief Sets up field dependencies.
  virtual void setupDependencies() = 0;

  //! \brief Advances the time step one step forward.
  virtual bool advanceStep(TimeStep& tp)
  {
    return S2.advanceStep(tp) && S1.advanceStep(tp);
  }

  //! \brief Computes the solution for the current time step.
  virtual bool solveStep(TimeStep& tp)
  {
    return S1.solveStep(tp) && S2.solveStep(tp);
  }

  //! \brief Saves the converged results to VTF-file of a given time step.
  //! \param[in] tp Time step identifier
  virtual bool saveStep(const TimeStep& tp, int& nBlock)
  {
    return S2.saveStep(tp, nBlock) && 
           S1.saveStep(tp, nBlock);
  }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  virtual bool saveModel(char* fileName, int& nBlock)
  {
    if (!S1.saveModel(fileName, nBlock))
      return false;

    S2.setVTF(S1.getVTF());

    return true;
  }

  //! \brief Returns a const reference to the time stepping parameters.
  const TimeStep& getTimePrm() const
  {
    return S1.getTimePrm();
  }

  //! \brief Registers a dependency on a field from another SIM object.
  //! \param[in] sim The SIM object holding the field we depend on
  //! \param[in] name Name of field we depend on
  //! \param[in] nvc Number of components in field
  //! \param[in] patches The geometry the field is defined over
  //! \param[in] diffBasis Different basis for the SIM class and the field
  virtual void registerDependency(SIMdependency* sim, const std::string& name,
                                  short int nvc,
                                  const SIMdependency::PatchVec& patches,
                                  bool diffBasis = false)
  {
    S1.registerDependency(sim, name, nvc, patches, diffBasis);
    S2.registerDependency(sim, name, nvc, patches, diffBasis);
  }
  //! \brief Registers a dependency on a field from another SIM object.
  //! \param[in] sim The SIM object holding the field we depend on
  //! \param[in] name Name of field we depend on
  //! \param[in] nvc Number of components in field
  void registerDependency(SIMdependency* sim, const std::string& name,
                          short int nvc = 1)
  {
    S1.registerDependency(sim, name, nvc);
    S2.registerDependency(sim, name, nvc);
  }

protected:
  T1& S1; //!< First substep
  T2& S2; //!< Second substep
};

#endif
