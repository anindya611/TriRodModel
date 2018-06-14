#include <TriRodsProperties.h>

namespace TriRods
{
  // Compute the strain energy density
  double Material::ComputeStrainEnergyDensity(const double Epsilon[2][2], const double Rho[2][2], const double Delta[2]) const
  {
    const double E = props.E;
    const double h = props.h;
    const double nu = props.nu;
    const double k0 = (E*h)/(1-nu*nu);
    const double k1 = (E*h*h*h)/(12*(1-nu*nu));
    const double k2 = (E*h)/(2*(1+nu));
    const double traceEps = Epsilon[1][1] + Epsilon[2][2];
    const double traceRho = Rho[1][1] + Rho[2][2];
    double Identity[2][2];
    for(int i=0; i<2; ++i)
      for(int j=0; j<2; ++j)
	Identity[i][j] = 0.;
    Identity[1][1] += 1.;
    Identity[2][2] += 1.;
    double n[2][2], m[2][2], q[2];
    double Energy = 0.;

    for(int i=0; i<2; ++i)
      {
	q[i] = k2*Delta[i];
	Energy += 0.5*q[i]*Delta[i];
	for(int j=0; j<2; ++j)
	  {
	    n[i][j] = k0*(nu*Identity[i][j]*traceEps + 0.5*(1-nu)*(Epsilon[i][j] + Epsilon[j][i]));
	    m[i][j] = k1*(nu*Identity[i][j]*traceRho + 0.5*(1-nu)*(Rho[i][j] + Rho[j][i]));
	    Energy += 0.5*(n[i][j]*Epsilon[i][j] + m[i][j]*Rho[i][j]);
	  }
      }
    return Energy; // Strain energy density
  }


  // Compute the stress/moment resultant
  void Material::ComputeStressResultants(const double Epsilon[2][2], const double Rho[2][2], const double Delta[2], double n[2][2], double m[2][2], double q[2]) const
  {
    const double E = props.E;
    const double h = props.h;
    const double nu = props.nu;
    const double k0 = (E*h)/(1-nu*nu);
    const double k1 = (E*h*h*h)/(12*(1-nu*nu));
    const double k2 = (E*h)/(2*(1+nu));
    const double traceEps = Epsilon[1][1] + Epsilon[2][2];
    const double traceRho = Rho[1][1] + Rho[2][2];
    double Identity[2][2];
    for(int i=0; i<2; ++i)
      for(int j=0; j<2; ++j)
	Identity[i][j] = 0.;
    Identity[1][1] += 1.;
    Identity[2][2] += 1.;

    for(int i=0; i<2; ++i)
      {
	q[i] = k2*Delta[i];
	for(int j=0; j<2; ++j)
	  {
	    n[i][j] = k0*(nu*Identity[i][j]*traceEps + 0.5*(1-nu)*(Epsilon[i][j] + Epsilon[j][i]));
	    m[i][j] = k1*(nu*Identity[i][j]*traceRho + 0.5*(1-nu)*(Rho[i][j] + Rho[j][i]));
	  }
      }
  }


  // Compute 1st variations of stress/moment resultants
  void Material::ComputeVarResultants(const double vEpsilon[2][2], const double vRho[2][2], const double vDelta[2], double vn[2][2], double vm[2][2], double vq[2]) const
  {
    const double E = props.E;
    const double h = props.h;
    const double nu = props.nu;
    const double k0 = (E*h)/(1-nu*nu);
    const double k1 = (E*h*h*h)/(12*(1-nu*nu));
    const double k2 = (E*h)/(2*(1+nu));
    const double traceVEps = vEpsilon[1][1] + vEpsilon[2][2];
    const double traceVRho = vRho[1][1] + vRho[2][2];
    double Identity[2][2];
    for(int i=0; i<2; ++i)
      for(int j=0; j<2; ++j)
	Identity[i][j] = 0.;
    Identity[1][1] += 1.;
    Identity[2][2] += 1.;

    for(int i=0; i<2; ++i)
      {
	vq[i] = k2*vDelta[i];
	for(int j=0; j<2; ++j)
	  {
	    vn[i][j] = k0*(nu*Identity[i][j]*traceVEps + 0.5*(1-nu)*(vEpsilon[i][j] + vEpsilon[j][i]));
	    vm[i][j] = k1*(nu*Identity[i][j]*traceVRho + 0.5*(1-nu)*(vRho[i][j] + vRho[j][i]));
	  }
      }
  }
}
