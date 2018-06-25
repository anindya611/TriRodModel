#ifndef TRI_RODS_PROPERTIES_H
#define TRI_RODS_PROPERTIES_H

#include <cassert>
#include <TriRodsPointKinematicsNew.h>

namespace TriRods
{
  // Encapsulation of ribbon information
  class MatProperties
  {
  public:
    double E; // young's modulus
    double h; // thickness of the ribbon 
    double nu; // poisson's ratio
    double w; // width of the ribbon/rod  ////

    // Constructor
    inline MatProperties(double v1, double v2, double v3, double v4)
      :E(v1), h(v2), nu(v3), w(v4) {}  ////

    // Copy constructor
    inline MatProperties(const MatProperties& Obj)
      :E(Obj.E), h(Obj.h), nu(Obj.nu), w(Obj.w) {}

    // Destructor
    inline virtual ~MatProperties() {}

    // Assignment operator
    inline MatProperties& operator=(const MatProperties& rhs)
      {
	if(&rhs==this) return *this;
	this->E = rhs.E;
	this->h = rhs.h;
	this->nu = rhs.nu;
	this->w = rhs.w; ////
	return *this;
      }

    ////
    // Returns the width of the ribbon
    double GetWidth()
    { return w; }
    ////
  };


  // Material for the trirods model
  class Material
  {
  public:
    // Constructor with provided material properties
    inline Material(const MatProperties mc)
      :props(mc) {}

    // Destructor, does nothing
    inline virtual ~Material() {}

    // Disable assignment operator
    Material& operator=(const Material) = delete;

    // Returns the set of material properties
    inline MatProperties GetMaterialProperties() const
    { return props; }

    // Compute the energy density
    double ComputeStrainEnergyDensity(const double Epsilon[2][2], const double Rho[2][2], const double Delta[2]) const;

    // Compute the stress/moment resultants
    void ComputeStressResultants(const double Epsilon[2][2], const double Rho[2][2], const double Delta[2], double n[2][2], double m[2][2], double q[2]) const;

    // Compute 1st variations of stress/moment resultants
    void ComputeVarResultants(const double vEpsilon[2][2], const double vRho[2][2], const double vDelta[2], double vn[2][2], double vm[2][2], double vq[2]) const;

    //

  private:
    const MatProperties props;
  };
}

#endif
