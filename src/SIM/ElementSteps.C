// $Id$
//==============================================================================
//!
//! \file ElementSteps.C
//!
//! \date Jun 30 2026
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Spatial element step function.
//!
//==============================================================================

#include "ElementSteps.h"
#include "SIMbase.h"
#include "ASMbase.h"
#include "Functions.h"
#include "Vec3Oper.h"
#include "IFEM.h"

#include <sstream>
#include <cstring>


ElementSteps::ElementSteps (const char* input, const SIMbase& sim, int nsd)
{
  if (!input || input[0] == 0)
    return; // avoid segfault on empty string

  size_t nc = strlen(input);
  char* cpy = strdup(input);
  // Replace all '\' and '|' characters in the string by newline '\n'
  for (size_t i = 0; i < nc; i++)
    if (cpy[i] == '\\' || cpy[i] == '|')
      cpy[i] = '\n';

  IFEM::cout <<" ElementSteps";
  std::stringstream str(cpy);
  char temp[512];
  while (str.getline(temp,512))
    if (temp[0] != '#' && temp[0] != 0)
    {
      std::stringstream sline(temp);
      double p[3] = { 0.0, 0.0, 0.0 };
      double value = 0.0;
      int patch = 1;
      for (int i = 0; i < nsd; i++)
        sline >> p[i];
      sline >> value >> patch;

      IFEM::cout <<"\n\t\tElement("<< p[0];
      for (int i = 1; i < nsd; i++)
        IFEM::cout <<", "<< p[i];
      IFEM::cout <<", "<< patch <<") = "<< value;

      ASMbase* pch = sim.getPatch(patch);
      if (!pch)
      {
        std::cerr <<"\n *** ElementSteps: No patch "<< patch << std::endl;
        continue;
      }
      int iel = pch->findElementContaining(p);
      if (iel < 1)
      {
        std::cerr <<"\n *** ElementSteps: Failed to locate element"<< std::endl;
        continue;
      }

      Vec3 X0, X1;
      if (pch->getElementBBox(X0,X1,iel))
      {
        IFEM::cout <<" -> inside([ " << X0 <<"] - ["<< X1 <<"])*"<< value;
        this->add(new StepXYZFunc(value,X0,X1,0.001));
      }
    }
  IFEM::cout << std::endl;

  free(cpy);
}
