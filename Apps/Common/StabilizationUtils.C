#include "StabilizationUtils.h"
#include "Vec3Oper.h"

namespace StabilizationUtils {

double getElementSize(const std::vector<Vec3>& XC, int nsd)
{
  double h;
  if (nsd == 2) {
    h = (XC[0]-XC[1]).length();
    h = std::min(h, (XC[1]-XC[3]).length());
    h = std::min(h, (XC[3]-XC[2]).length());
    h = std::min(h, (XC[2]-XC[0]).length());
  } else {
    static const int comps[][2] = {{0,1}, {1,5}, {5,4}, {4,0},
                                   {2,3}, {3,7}, {7,6}, {6,2},
                                   {0,4}, {4,6},        {2,0},
                                          {5,7},        {3,1}};
    h = 1e8;
    for (size_t i=0; i < sizeof(comps)/(2*sizeof(int));++i)
      h = std::min(h, (XC[comps[i][0]]-XC[comps[i][1]]).length());
  }

  return h;
}


double getTauPt(double dt, double mu, const Vector& U, const Matrix& G, const double Ct, const double Cl)
{
  double Gnorm2 = 0.0;
  for (size_t i = 1;i <= G.rows();i++)
    for (size_t j = 1;j <= G.cols();j++)  
      Gnorm2 += G(i,j)*G(i,j);

  return 1.0/sqrt( Ct/pow(dt,2) + U.dot(G*U) + Cl*mu*mu*Gnorm2);
}


bool getTauNSPt(double dt, double mu, const Vector& U, const Matrix& G, double & tauM, double& tauC, const double Ct, const double Cl)
{
  tauM = getTauPt(dt,mu,U,G,Ct,Cl);

  double Gtrace = 0.0;
  for (size_t i = 1;i <= G.rows();i++)
    Gtrace += G(i,i);

  tauC = 1.0/(tauM*Gtrace);

  return true;
}


bool getTauNSALEPt(double dt, double mu, const Vector& U, const Matrix& G, double & tauM, double& tauC, const double Ct, const double Cl)
{
  tauM = getTauPt(dt,mu,U,G,Ct,Cl);

  tauC = tauM*(U.dot(U));

  return true;
}


bool getTauPtJac(const Vector& U, const Matrix& G, const double tauM, Vector& tauMjac)
{
  tauMjac = -pow(tauM,3)*G*U;

  return true;
}


bool getTauNSPtJac(const Vector& U, const Matrix& G, const double tauM, const double& tauC, Vector& tauMjac, Vector& tauCjac)
{
  if (!getTauPtJac(U, G, tauM, tauMjac))
    return false;

  tauCjac = (-tauC/tauM) * tauMjac;

  return true;
}

}
