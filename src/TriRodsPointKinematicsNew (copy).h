#ifndef TRI_RODS_POINT_KINEMATICS_NEW_H
#define TRI_RODS_POINT_KINEMATICS_NEW_H

#include<cassert>

namespace TriRods
{
  // Encapsulation of kinematical data at a quadrature point
  class PointKinematics
  {
  public:
    double phi0[3], phi1[3], phi2[3], dirSpatial[3];
    double d1Phi0[3], d1Phi1[3], d1Phi2[3], d1DirSpatial[3];
    double Lambda[3][3];
    double a1[3], a2[3];
    double n[2][2], m[2][2], q[2];
    double epsilon[2][2], rho[2][2], delta[2];

    //! Default constructor
    inline PointKinematics() {}

    // Copy constructor
    inline PointKinematics(const PointKinematics& Obj)
    {
      for(int i=0; i<3; ++i)
	{
	  phi0[i] = Obj.phi0[i];
	  phi1[i] = Obj.phi1[i];
	  phi2[i] = Obj.phi2[i];
	  dirSpatial[i] = Obj.dirSpatial[i];
	  d1Phi0[i] = Obj.d1Phi0[i];
	  d1Phi1[i] = Obj.d1Phi1[i];
	  d1Phi2[i] = Obj.d1Phi2[i];
	  d1DirSpatial[i] = Obj.d1DirSpatial[i];
	  a1[i] = Obj.a1[i];
	  a2[i] = Obj.a2[i];
	  for(int j=0; j<3; ++j)
	    {
	      Lambda[i][j] = Obj.Lambda[i][j];
	    }
	}
      for(int i=0; i<2; ++i)
	{
	  q[i] = Obj.q[i];
	  delta[i] = Obj.delta[i];
	  for(int j=0; j<2; ++j)
	    {
	      n[i][j] = Obj.n[i][j];
	      m[i][j] = Obj.m[i][j];
	      epsilon[i][j] = Obj.epsilon[i][j];
	      rho[i][j] = Obj.rho[i][j];
	    }
	}
    }

    // Destructor, does nothing
    inline virtual ~PointKinematics() {}

    // Assignment operator
    inline PointKinematics& operator=(const PointKinematics& rhs)
      {
	if(&rhs==this) return *this;

	for(int i=0; i<3; ++i)
	  {
	    phi0[i] = rhs.phi0[i];
	    phi1[i] = rhs.phi1[i];
	    phi2[i] = rhs.phi2[i];
	    dirSpatial[i] = rhs.dirSpatial[i];
	    d1Phi0[i] = rhs.d1Phi0[i];
	    d1Phi1[i] = rhs.d1Phi1[i];
	    d1Phi2[i] = rhs.d1Phi2[i];
	    d1DirSpatial[i] = rhs.d1DirSpatial[i];
	    a1[i] = rhs.a1[i];
	    a2[i] = rhs.a2[i];
	    for(int j=0; j<3; ++j)
	      {
		Lambda[i][j] = rhs.Lambda[i][j];
	      }
	  }
	for(int i=0; i<2; ++i)
	  {
	    q[i] = rhs.q[i];
	    delta[i] = rhs.delta[i];
	    for(int j=0; j<2; ++j)
	      {
		n[i][j] = rhs.n[i][j];
		m[i][j] = rhs.m[i][j];
		epsilon[i][j] = rhs.epsilon[i][j];
		rho[i][j] = rhs.rho[i][j];
	      }
	  }
	return *this;
      }

    // Compute lambda // doesn't work when t = -E3
    inline void ComputeRotationMatrix()
    {
      double E3[3] = {0.,0.,1.};
      double tempvec[3];
      tempvec[0] = E3[1]*dirSpatial[2] - E3[2]*dirSpatial[1];
      tempvec[1] = E3[2]*dirSpatial[0] - E3[0]*dirSpatial[2];
      tempvec[2] = E3[0]*dirSpatial[1] - E3[1]*dirSpatial[0];
      double tempmat1[3][3];
      for(int i=0; i<3; ++i)
	for(int j=0; j<3; ++j)
	  tempmat1[i][j] = 0.;
      tempmat1[0][1] += -tempvec[2];
      tempmat1[0][2] += tempvec[1];
      tempmat1[1][2] += -tempvec[0];
      tempmat1[1][0] += tempvec[2];
      tempmat1[2][0] += -tempvec[1];
      tempmat1[2][1] += tempvec[0];
      double cross[3];
      for(int i=0; i<3; ++i)
	{
	  cross[i] = 0.;
	  for(int j=0; j<3; ++j)
	    {
	      cross[i] += tempmat1[i][j]*tempvec[j];
	    }
	}
      for(int i=0; i<3; ++i)
	assert(cross[i] < 1.e-10);
      double tempmat2[3][3];
      for(int i=0; i<3; ++i)
	for(int j=0; j<3; ++j)
	  tempmat2[i][j] = tempvec[i]*tempvec[j];
      double temp = dirSpatial[2]; // temp = E3.dirSpatial
      double identity[3][3] = {{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};

      // Calculate Lambda
      for(int i=0; i<3; ++i)
	for(int j=0; j<3; ++j)
	  Lambda[i][j] = temp*identity[i][j] + tempmat1[i][j] + (1/(1+temp))*tempmat2[i][j];
    }

    // Compute the strain measures
    inline void ComputeStrains(const double y, const double w)
    {
      double a01[3], a02[3];
      double alpha[2][2], alpha0[2][2], kappa[2][2], kappa0[2][2];
      double gamma[2], gamma0[2];
      for(int i=0; i<3; ++i)
	{
	  a1[i] = d1Phi0[i] + ((d1Phi2[i] - d1Phi1[i])/w)*y + 2*((d1Phi1[i] + d1Phi2[i] - 2*d1Phi0[i])/(w*w))*(y*y);
	  a2[i] = ((phi2[i] - phi1[i])/w) + 4*((phi2[i] + phi1[i] - 2*phi0[i])/(w*w))*y;
	  a01[i] = 0.;
	  a02[i] = 0.;
	}
      a01[0] += 1.;
      a02[1] += 1.;

      alpha[0][0] = 0.;
      alpha[0][1] = 0.;
      alpha[1][0] = 0.;
      alpha[1][1] = 0.;

      alpha0[0][0] = 0.;
      alpha0[0][1] = 0.;
      alpha0[1][0] = 0.;
      alpha0[1][1] = 0.;

      kappa[0][0] = 0.;
      kappa[0][1] = 0.;
      kappa[1][0] = 0.;
      kappa[1][1] = 0.;

      kappa0[0][0] = 0.;
      kappa0[0][1] = 0.;
      kappa0[1][0] = 0.;
      kappa0[1][1] = 0.;

      gamma[0] = 0.;
      gamma[1] = 0.;

      gamma0[0] = 0.;
      gamma0[1] = 0.;

      for(int i=0; i<3; ++i)
	{
	  alpha[0][0] += a1[i]*a1[i];
	  alpha[0][1] += a1[i]*a2[i];
	  alpha[1][0] += a2[i]*a1[i];
	  alpha[1][1] += a2[i]*a2[i];

	  alpha0[0][0] += a01[i]*a01[i];
	  alpha0[0][1] += a01[i]*a02[i];
	  alpha0[1][0] += a02[i]*a01[i];
	  alpha0[1][1] += a02[i]*a02[i];

	  kappa[0][0] += a1[i]*d1DirSpatial[i];
	  //kappa[0][1] = 0.;
	  kappa[1][0] += a2[i]*d1DirSpatial[i];
	  //kappa[1][1] = 0.;

	  //kappa0[0][0] = 0.;
	  //kappa0[0][1] = 0.;
	  //kappa0[1][0] = 0.;
	  //kappa0[1][1] = 0.;

	  gamma[0] += a1[i]*dirSpatial[i];
	  gamma[1] += a2[i]*dirSpatial[i];

	  //gamma0[0] = 0.;
	  //gamma0[1] = 0.;
	}

      for(int i=0; i<2; ++i)
	{
	  delta[i] = gamma[i] - gamma0[i];
	  for(int j=0; j<2; ++j)
	    {
	      epsilon[i][j] = 0.5*(alpha[i][j] - alpha0[i][j]);
	      rho[i][j] = kappa[i][j] - kappa0[i][j];
	    }
	}
      return;
    }
  };



  // First variation of kinematic data at a quadrature point
  class VarPointKinematics
  {
  public:
    const PointKinematics* PD;
    double eta0[3], eta1[3], eta2[3], delDirMaterial[2]; // delDirSpatial[3]; ////
    double d1Eta0[3], d1Eta1[3], d1Eta2[3], d1DelDirSpatial[3]; 
    double va1[3], va2[3];
    double vn[2][2], vm[2][2], vq[2]; // 1st variation of stress/moment resultants
    double vEpsilon[2][2], vRho[2][2], vDelta[2]; // 1st variation of strain measures

    // Default constructor
    inline VarPointKinematics()
      :PD(nullptr)
    { SetZero(); }

    // Copy Constructor
    inline VarPointKinematics(const VarPointKinematics& Obj)
      :PD(Obj.PD)
    {
      for(int i=0; i<3; ++i)
	{
	  eta0[i] = Obj.eta0[i];
	  eta1[i] = Obj.eta1[i];
	  eta2[i] = Obj.eta2[i];
	  //delDirSpatial[i] = Obj.delDirSpatial[i]; ////
	  d1Eta0[i] = Obj.d1Eta0[i];
	  d1Eta1[i] = Obj.d1Eta1[i];
	  d1Eta2[i] = Obj.d1Eta2[i];
	  d1DelDirSpatial[i] = Obj.d1DelDirSpatial[i]; 
	  va1[i] = Obj.va1[i]; 
	  va2[i] = Obj.va2[i]; 
	}
      for(int i=0; i<2; ++i)
	{
	  delDirMaterial[i] = Obj.delDirMaterial[i]; ////
	  vq[i] = Obj.vq[i];
	  vDelta[i] = Obj.vDelta[i];
	  for(int j=0; j<2; ++j)
	    {
	      vn[i][j] = Obj.vn[i][j];
	      vm[i][j] = Obj.vm[i][j];
	      vEpsilon[i][j] = Obj.vEpsilon[i][j];
	      vRho[i][j] = Obj.vRho[i][j];
	    }
	}
    }

    // Destructor, does nothing
    inline virtual ~VarPointKinematics() {}

    //  Assignment operator
    inline VarPointKinematics& operator=(const VarPointKinematics& rhs)
      { assert(false && "Should not need to assign variation of point data"); }
    
    // Set all the values to zero
    inline void SetZero()
    {
      for(int i=0; i<3; ++i)
	{
	  eta0[i] = 0.;
	  eta1[i] = 0.;
	  eta2[i] = 0.;
	  //delDirSpatial[i] = 0.; ////
	  d1Eta0[i] = 0.;
	  d1Eta1[i] = 0.;
	  d1Eta2[i] = 0.;
	  d1DelDirSpatial[i] = 0.;
	  va1[i] = 0; 
	  va2[i] = 0; 
	}
      for(int i=0; i<2; ++i)
	{
	  delDirMaterial[i] = 0.; ////
	  vq[i] = 0.;
	  vDelta[i] = 0.;
	  for(int j=0; j<2; ++j)
	    {
	      vn[i][j] = 0.;
	      vm[i][j] = 0.;
	      vEpsilon[i][j] = 0.;
	      vRho[i][j] = 0.;
	    }
	}
    }

    // Compute the 1st variations of strain measures
    inline void ComputeVarStrains(const double y, const double w)
    {
      const auto* a1 = PD->a1;
      const auto* a2 = PD->a2;
      const auto* dirSpatial = PD->dirSpatial;
      const auto* d1DirSpatial = PD->d1DirSpatial;
      const auto* Lambda = PD->Lambda; ////

      double delDirMaterial3Comps[3], delDirSpatial[3]; ////
      for(int i=0; i<2; ++i) ////
	delDirMaterial3Comps[i] = delDirMaterial[i]; ////
      delDirMaterial3Comps[2] = 0.; ////

      ////
      for(int i=0; i<3; ++i) 
	{
	  delDirSpatial[i] = 0.; 
	  for(int j=0; j<3; ++j) 
	    {
	      delDirSpatial[i] += Lambda[i][j]*delDirMaterial3Comps[j]; 
	    }
	}
      ////

      for(int i=0; i<3; ++i)
	{
	  va1[i] = d1Eta0[i] + ((d1Eta2[i] - d1Eta1[i])/w)*y + (2*y*y)*((d1Eta2[i] + d1Eta1[i] - 2*d1Eta0[i])/(w*w));
	  va2[i] = (eta2[i] - eta1[i])/w + 4*((eta2[i] + eta1[i] - 2*eta0[i])/(w*w))*y;
	}
      
      vEpsilon[0][0] = 0.;
      vEpsilon[0][1] = 0.;
      vEpsilon[1][0] = 0.;
      vEpsilon[1][1] = 0.;

      vRho[0][0] = 0.;
      vRho[0][1] = 0.;
      vRho[1][0] = 0.;
      vRho[1][1] = 0.;

      vDelta[0] = 0.;
      vDelta[1] = 0.;
	  
      for(int i=0; i<3; ++i)
	{
	  vEpsilon[0][0] += va1[i]*a1[i];
	  vEpsilon[0][1] += 0.5*(va1[i]*a2[i] + a1[i]*va2[i]);
	  vEpsilon[1][0] += 0.5*(va2[i]*a1[i] + a2[i]*va1[i]);
	  vEpsilon[1][1] += va2[i]*a2[i];

	  vRho[0][0] += va1[i]*d1DirSpatial[i] + a1[i]*d1DelDirSpatial[i];
	  //vRho[0][1] = 0.;
	  vRho[1][0] += va2[i]*d1DirSpatial[i] + a2[i]*d1DelDirSpatial[i];
	  //vRho[1][1] = 0.;

	  vDelta[0] += va1[i]*dirSpatial[i] + a1[i]*delDirSpatial[i];
	  vDelta[1] += va2[i]*dirSpatial[i] + a2[i]*delDirSpatial[i];
	}
      return;
    }
  };



  // 2nd variation of kinematic data at a quadrature point
  class VarVarPointKinematics
  {
  public:
    const PointKinematics* PD; // kinematic data at this point
    const VarPointKinematics* deltaPD; // Variation 1
    const VarPointKinematics* DELTAPD; // Variation 2
    double vvEpsilon[2][2], vvRho[2][2], vvDelta[2]; // 2nd variations of strain measures

    // Default constructor
    inline VarVarPointKinematics()
      :PD(nullptr), deltaPD(nullptr), DELTAPD(nullptr)
    {
      for(int i=0; i<2; ++i)
	{
	  vvDelta[i] = 0.;
	  for(int j=0; j<2; ++j)
	    {
	      vvEpsilon[i][j] = 0.;
	      vvRho[i][j] = 0.;
	    }
	}
    }

    // Copy constructor
    inline VarVarPointKinematics(const VarVarPointKinematics& Obj)
      :PD(Obj.PD), deltaPD(Obj.deltaPD), DELTAPD(Obj.DELTAPD)
      {
	for(int i=0; i<2; ++i)
	{
	  vvDelta[i] = Obj.vvDelta[i];
	  for(int j=0; j<2; ++j)
	    {
	      vvEpsilon[i][j] = Obj.vvEpsilon[i][j];
	      vvRho[i][j] = Obj.vvRho[i][j];
	    }
	}
      }

    // Destructor, does nothing
    inline virtual ~VarVarPointKinematics() {}

    // Assignment operator: disable
    inline VarVarPointKinematics* operator=(const VarVarPointKinematics& rhs)
      { assert(false && "Assignment should not be required."); }
    
    // Compute the 2nd variations of strain measures
    inline void ComputeVarVarStrains()
    {
      // Aliases from kinematics
      const auto* a1 = PD->a1;
      const auto* a2 = PD->a2;
      const auto* t = PD->dirSpatial;
      const auto* d1t = PD->d1DirSpatial;
      const auto* Lambda = PD->Lambda; ////

      // Aliases for the 1st variation
      const auto* va1 = deltaPD->va1;
      const auto* va2 = deltaPD->va2;
      //const auto* delt = deltaPD->delDirSpatial; ////
      const auto* delT = deltaPD->delDirMaterial; ////
      const auto* d1delt = deltaPD->d1DelDirSpatial;

      ////
      double delT3[3], delt[3];
      for(int i=0; i<2; ++i)
	delT3[i] = delT[i];
      delT3[2] = 0.;
      for(int i=0; i<3; ++i)
	{
	  delt[i] = 0.;
	  for(int j=0; j<3; ++j)
	    delt[i] += Lambda[i][j]*delT3[j];
	}
      ////

      // Aliases for another 1st variation
      const auto* Va1 = DELTAPD->va1;
      const auto* Va2 = DELTAPD->va2;
      //const auto* DELt = DELTAPD->delDirSpatial; ////
      const auto* DELT = DELTAPD->delDirMaterial; ////
      const auto* d1DELt = DELTAPD->d1DelDirSpatial;

      ////
      double DELT3[3], DELt[3];
      for(int i=0; i<2; ++i)
	DELT3[i] = DELT[i];
      DELT3[2] = 0.;
      for(int i=0; i<3; ++i)
	{
	  DELt[i] = 0.;
	  for(int j=0; j<3; ++j)
	    DELt[i] += Lambda[i][j]*DELT3[j];
	}
      ////

      double tempvec1[3];
      double tempvec2[3];
      
      for(int i=0; i<3; ++i)
	{
	  tempvec1[i] = 0.;
	  tempvec2[i] = 0.;
	  for(int j=0; j<3; ++j)
	    {
	      tempvec1[i] += d1t[i]*delt[j]*DELt[j] + t[i]*d1delt[j]*DELt[j] + t[i]*delt[j]*d1DELt[j];
	      tempvec2[i] += t[i]*delt[j]*DELt[j];
	    }
	}
      
      // (d1t[i]*delt[j]*DELt[j] + t[i]*d1delt[j]*DELt[j] + t[i]*delt[j]*d1DELt[j])
      // (t[i]*delt[j]*DELt[j])

      vvEpsilon[0][0] = 0.;
      vvEpsilon[0][1] = 0.;
      vvEpsilon[1][0] = 0.;
      vvEpsilon[1][1] = 0.;
      
      vvRho[0][0] = 0.;
      vvRho[0][1] = 0.;
      vvRho[1][0] = 0.;
      vvRho[1][1] = 0.;
      
      vvDelta[0] = 0.;
      vvDelta[1] = 0.;

      // 2nd variation of Epsilon
      for(int i=0; i<3; ++i)
	{
	  vvEpsilon[0][0] += va1[i]*Va1[i];
	  vvEpsilon[0][1] += 0.5*(va1[i]*Va2[i] + Va1[i]*va2[i]);
	  vvEpsilon[1][0] += 0.5*(va2[i]*Va1[i] + Va2[i]*va1[i]);
	  vvEpsilon[1][1] += va2[i]*Va2[i];
	  
	  vvRho[0][0] += va1[i]*d1DELt[i] + Va1[i]*d1delt[i] - a1[i]*tempvec1[i]; 
	  //vvRho[0][1] = 0.;
	  vvRho[1][0] += va2[i]*d1DELt[i] + Va2[i]*d1delt[i] - a2[i]*tempvec1[i]; 
	  //vvRho[1][1] = 0.;

	  vvDelta[0] += va1[i]*DELt[i] + Va1[i]*delt[i] - a1[i]*tempvec2[i]; 
	  vvDelta[1] += va2[i]*DELt[i] + Va2[i]*delt[i] - a2[i]*tempvec2[i]; 
	}
      return;
    }
  };
}
  
#endif
