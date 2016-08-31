#include "FiniteElement.h"
#include <iostream>

std::ostream& operator<<(std::ostream& os, const FiniteElement& obj)
{
  // write obj to stream
  os << "iGP: "    << obj.iGP << std::endl;
  os << "u: "      << obj.u << std::endl;
  os << "v: "      << obj.v << std::endl;
  os << "w: "      << obj.w << std::endl;
  os << "xi: "     << obj.xi << std::endl;
  os << "eta: "    << obj.eta << std::endl;
  os << "zeta: "   << obj.zeta << std::endl;
  os << "detJxW: " << obj.detJxW << std::endl;
  os << "N: "      << obj.N << std::endl;
  os << "dNdX: "   << obj.dNdX << std::endl;
  os << "d2NdX2: " << obj.d2NdX2 << std::endl;
  os << "G: "      << obj.G << std::endl;

  // Element quantities
  os << "iel: "  << obj.iel << std::endl;
  os << "p: "    << obj.p << std::endl;
  os << "Navg: " << obj.Navg << std::endl;
  os << "Xn: "   << obj.Xn << std::endl;
  os << "Te: "   << obj.Te << std::endl;
  return os;
}

