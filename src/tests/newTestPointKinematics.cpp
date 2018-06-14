#include <iostream>
#include <TriRodsPointKinematicsNew.h>
#include <random>
#include <cmath>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_bessel.h>

// Random double
std::random_device rd; // will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis(-1.,1.);

using namespace TriRods;

// Initialize values
void Initialize(PointKinematics& PD);

// Update state
void Update(PointKinematics& PD, const VarPointKinematics& var);

// Compute numerical values for first variations
void ComputeNumericalVariations(const PointKinematics& PD,
				const VarPointKinematics& var,
				const double EPS,
				double vEpsilon[2][2],
				double vRho[2][2],
				double vDelta[2]);

// Compute numerical values for second variations
void ComputeNumerical2ndVariations(const PointKinematics& PD,
				   const VarPointKinematics& var,
				   const VarPointKinematics& VAR,
				   const double EPS,
				   double vvEpsilon[2][2],
				   double vvRho[2][2],
				   double vvDelta[2]);

// Generate d1DelDirSpatia
void GenerateRandomD1DelDirSpatial(const PointKinematics& PD,
		                   const VarPointKinematics& var,
				   double& k,
				   double d1DelDirSpatial[3]);

// Toggle to print results
const bool PRINT = false;

// main function
int main()
{
  // Initialize point data
  PointKinematics PD;
  Initialize(PD);

  //////////
  std::cout<<'\n'<<'\n'<<"for PD=========="<<'\n';
  std::cout<<std::abs(PD.dirSpatial[0]*PD.d1DirSpatial[0] + PD.dirSpatial[1]*PD.d1DirSpatial[1] + PD.dirSpatial[2]*PD.d1DirSpatial[2])<<'\n';
  assert(std::abs(PD.dirSpatial[0]*PD.d1DirSpatial[0] + PD.dirSpatial[1]*PD.d1DirSpatial[1] + PD.dirSpatial[2]*PD.d1DirSpatial[2]) < 1.e-10);
  /////////
  
  const double y = 0.1;
  const double w = 0.2;

  // Compute the strains
  PD.ComputeStrains(y,w);

  // Set up two independent variations
  VarPointKinematics varA;
  varA.PD = &PD;
  for(int i=0; i<3; ++i)
    {
      varA.eta0[i] = dis(gen);
      varA.eta1[i] = dis(gen);
      varA.eta2[i] = dis(gen);
      //varA.delDirSpatial[i] = dis(gen); ////
      varA.d1Eta0[i] = dis(gen);
      varA.d1Eta1[i] = dis(gen);
      varA.d1Eta2[i] = dis(gen);
      //varA.d1DelDirSpatial[i] = dis(gen); ////
    }

  ////
  for(int i=0; i<2; ++i)
    varA.delDirMaterial[i] = dis(gen);
  ////

  varA.ComputeDelDirSpatial();

  ////
  double kA = 0.;
  GenerateRandomD1DelDirSpatial(PD, varA, kA, varA.d1DelDirSpatial);
  ////

  //////////
  std::cout<<'\n'<<'\n'<<"for varA=========="<<'\n';
  std::cout<<std::abs(kA+(PD.dirSpatial[0]*varA.d1DelDirSpatial[0]+PD.dirSpatial[1]*varA.d1DelDirSpatial[1]+PD.dirSpatial[2]*varA.d1DelDirSpatial[2]))<<'\n'<<'\n';
  assert(std::abs(kA+(PD.dirSpatial[0]*varA.d1DelDirSpatial[0]+PD.dirSpatial[1]*varA.d1DelDirSpatial[1]+PD.dirSpatial[2]*varA.d1DelDirSpatial[2]))<1.e-10);
  /////////
  
  
  varA.ComputeVarStrains(y,w);
 
  VarPointKinematics varB;
  varB.PD = &PD;
  for(int i=0; i<3; ++i)
    {
      varB.eta0[i] = dis(gen);
      varB.eta1[i] = dis(gen);
      varB.eta2[i] = dis(gen);
      //varB.delDirSpatial[i] = dis(gen); ////
      varB.d1Eta0[i] = dis(gen);
      varB.d1Eta1[i] = dis(gen);
      varB.d1Eta2[i] = dis(gen);
      //varB.d1DelDirSpatial[i] = dis(gen); ////
    }

  ////
  for(int i=0; i<2; ++i)
    varB.delDirMaterial[i] = dis(gen);
  ////

  varB.ComputeDelDirSpatial();

  ////
  double kB = 0.;
  GenerateRandomD1DelDirSpatial(PD, varB, kB, varB.d1DelDirSpatial);
  ////

  //////////
  std::cout<<'\n'<<'\n'<<"for varB============="<<'\n';
  std::cout<<std::abs(kB+(PD.dirSpatial[0]*varB.d1DelDirSpatial[0]+PD.dirSpatial[1]*varB.d1DelDirSpatial[1]+PD.dirSpatial[2]*varB.d1DelDirSpatial[2]))<<'\n'<<'\n';
  assert(std::abs(kB+(PD.dirSpatial[0]*varB.d1DelDirSpatial[0]+PD.dirSpatial[1]*varB.d1DelDirSpatial[1]+PD.dirSpatial[2]*varB.d1DelDirSpatial[2]))<1.e-10);
  /////////
  
  varB.ComputeVarStrains(y,w);

  // Compute the 1st varition numerically
  const double EPS = 1.e-5;

  // Variation A
  double varAEpsilon[2][2], varARho[2][2], varADelta[2];
  ComputeNumericalVariations(PD, varA, EPS, varAEpsilon, varARho, varADelta);
  
  if(PRINT)
    {
      std::cout<<"\nVariation of Epsilon: ";
      for(int i=0; i<2; ++i)
	for(int j=0; j<2; ++j)
	  std::cout<<"\n"<<varA.vEpsilon[i][j]<<" should be "<<varAEpsilon[i][j];
      std::cout<<"\nVariation of Rho: ";
      for(int i=0; i<2; ++i)
	for(int j=0; j<2; ++j)
	  std::cout<<"\n"<<varA.vRho[i][j]<<" should be "<<varARho[i][j];
      std::cout<<"\nVariation of Delta: ";
      for(int i=0; i<2; ++i)
	std::cout<<"\n"<<varA.vDelta[i]<<" should be "<<varADelta[i];
      std::cout<<"\n\n"; std::fflush( stdout );
    }
  for(int i=0; i<2; ++i)
    for(int j=0; j<2; ++j)
      assert((std::abs(varA.vEpsilon[i][j]-varAEpsilon[i][j]) + std::abs(varA.vRho[i][j]-varARho[i][j]) + std::abs(varA.vDelta[i]-varADelta[i]))<1.e-4); 

  // Variation B
  double varBEpsilon[2][2], varBRho[2][2], varBDelta[2];
  ComputeNumericalVariations(PD, varB, EPS, varBEpsilon, varBRho, varBDelta);

  // cout //
  std::cout<<"\n=======";
  std::cout<<"\n1st variations are OK \n";
  // // // //

  if(PRINT)
    {
      std::cout<<"\nVariation of Epsilon: ";
      for(int i=0; i<2; ++i)
	for(int j=0; j<2; ++j)
	  std::cout<<"\n"<<varB.vEpsilon[i][j]<<" should be "<<varBEpsilon[i][j];
      std::cout<<"\nVariation of Rho: ";
      for(int i=0; i<2; ++i)
	for(int j=0; j<2; ++j)
	  std::cout<<"\n"<<varB.vRho[i][j]<<" should be "<<varBRho[i][j];
      std::cout<<"\nVariation of Delta: ";
      for(int i=0; i<2; ++i)
	std::cout<<"\n"<<varB.vDelta[i]<<" should be "<<varBDelta[i];
      std::cout<<"\n\n"; std::fflush( stdout );
    }
  for(int i=0; i<2; ++i)
    for(int j=0; j<2; ++j)
      assert((std::abs(varB.vEpsilon[i][j]-varBEpsilon[i][j]) + std::abs(varB.vRho[i][j]-varBRho[i][j]) + std::abs(varB.vDelta[i]-varBDelta[i]))<1.e-4); 
  
  // Second variation
  VarVarPointKinematics varAB;
  varAB.PD = &PD;
  varAB.deltaPD = &varA;  
  varAB.DELTAPD = &varB;
  varAB.ComputeVarVarStrains();

  // Check consistency of 2nd variations
  double vvEpsilon[2][2], vvRho[2][2], vvDelta[2];
  ComputeNumerical2ndVariations(PD, varA, varB, EPS, vvEpsilon, vvRho, vvDelta);
  if(true) // PRINT
    {
      std::cout<<"\n2nd variation of Epsilon: ";
      for(int i=0; i<2; ++i)
	for(int j=0; j<2; ++j)
	  std::cout<<"\n"<<varAB.vvEpsilon[i][j]<<" should be "<<vvEpsilon[i][j];
      std::cout<<"\n2nd variation of Rho: ";
      for(int i=0; i<2; ++i)
	for(int j=0; j<2; ++j)
	  std::cout<<"\n"<<varAB.vvRho[i][j]<<" should be "<<vvRho[i][j];
      std::cout<<"\n2nd variation of Delta: ";
      for(int i=0; i<2; ++i)
	std::cout<<"\n"<<varAB.vvDelta[i]<<" should be "<<vvDelta[i];
      std::cout<<"\n\n"; std::fflush( stdout );
    }
  for(int i=0; i<2; ++i)
    for(int j=0; j<2; ++j)
      assert((std::abs(varAB.vvEpsilon[i][j]-vvEpsilon[i][j]) + std::abs(varAB.vvRho[i][j]-vvRho[i][j]) + std::abs(varAB.vvDelta[i]-vvDelta[i]))<1.e-4); 

}


// Returns the axial vector of a skew symmetric matrix
void AxialVector(const double E[][3], double* e)
{ e[0] = E[2][1]; e[1] = E[0][2]; e[2] = E[1][0]; }


// Returns the hodge star of a vector
void HodgeStar(const double* e, double E[][3])
{
  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      E[i][j] = 0.;
  E[0][1] = -e[2]; E[1][0] = e[2];
  E[0][2] =  e[1]; E[2][0] = -e[1];
  E[2][1] =  e[0]; E[1][2] = -e[0];

  // Check that the axial vector of E is e.
  double axial[3];
  AxialVector(E, axial);
  assert(std::abs(axial[0]-e[0])+std::abs(axial[1]-e[1])+std::abs(axial[2]-e[2])<1.e-6);
  return;
}



// Initialize values
void Initialize(PointKinematics& PD)
{
  // Set random values for phi0, phi1, phi2, t, phi0', phi1', phi2', t'
  for(int i=0; i<3; ++i)
    {
      PD.phi0[i] = dis(gen);
      PD.phi1[i] = dis(gen);
      PD.phi2[i] = dis(gen);
      PD.d1Phi0[i] = dis(gen);
      PD.d1Phi1[i] = dis(gen);
      PD.d1Phi2[i] = dis(gen);
      PD.dirSpatial[i] = dis(gen);
      //PD.d1DirSpatial[i] = dis(gen); ////
    }
  double normDirSpatial = sqrt(PD.dirSpatial[0]*PD.dirSpatial[0] + PD.dirSpatial[1]*PD.dirSpatial[1] + PD.dirSpatial[2]*PD.dirSpatial[2]);

  for(int i=0; i<3; ++i)
    PD.dirSpatial[i] /= normDirSpatial;
  
  ////
  double arbvec[3];
  for(int i=0; i<3; ++i)
    arbvec[i] = dis(gen);
  PD.d1DirSpatial[0] = PD.dirSpatial[1]*arbvec[2] - PD.dirSpatial[2]*arbvec[1];
  PD.d1DirSpatial[1] = PD.dirSpatial[2]*arbvec[0] - PD.dirSpatial[0]*arbvec[2];
  PD.d1DirSpatial[2] = PD.dirSpatial[0]*arbvec[1] - PD.dirSpatial[1]*arbvec[0];
  ////
  
  PD.ComputeRotationMatrix(); ////
  
  /*
  // Generate a random rotation tensor using the axis-angle representation
  double theta = dis(gen);
  double e[3];
  double norm = 0.;
  for(int i=0; i<3; ++i)
    {
      e[i] = dis(gen);
      norm +=e[i]*e[i];
    }
  norm = sqrt(norm);
  for(int i=0; i<3; ++i)
    e[i] /= norm;
  
  // Create the rotation matrix
  double StarE[3][3];
  HodgeStar(e, StarE);
  double delta[3][3] = {{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      PD.Lambda[i][j] =
	e[i]*e[j] + (delta[i][j]-e[i]*e[j])*(std::cos(theta)) + StarE[i][j]*(std::sin(theta));
  */
  
  // Check that the rotation matrix has determinant 1
  double det = 0.;
  for(int i=0; i<3; ++i)
    det += PD.Lambda[0][i]*
      (PD.Lambda[1][(i+1)%3]*PD.Lambda[2][(i+2)%3]-PD.Lambda[1][(i+2)%3]*PD.Lambda[2][(i+1)%3]);
  
  assert(std::abs(det-1.)<1.e-6);
}


////
void GenerateRandomD1DelDirSpatial(const PointKinematics& PD,
		                   const VarPointKinematics& var,
				   double& k,
				   double d1DelDirSpatial[3])
{
  /*
  double delDirMaterial3Comps[3], delDirSpatial[3];
  for(int i=0; i<2; ++i)
    delDirMaterial3Comps[i] = var.delDirMaterial[i];
  delDirMaterial3Comps[2] = 0.;
  for(int i=0; i<3; ++i)
    {
      delDirSpatial[i] = 0.;
      for(int j=0; j<3; ++j)
	delDirSpatial[i] += PD.Lambda[i][j]*delDirMaterial3Comps[j];
    }
  */
  k = var.delDirSpatial[0]*PD.d1DirSpatial[0] + var.delDirSpatial[1]*PD.d1DirSpatial[1] + var.delDirSpatial[2]*PD.d1DirSpatial[2];
  d1DelDirSpatial[0] = dis(gen);
  d1DelDirSpatial[1] = dis(gen);
  d1DelDirSpatial[2] = (-k - PD.dirSpatial[0]*d1DelDirSpatial[0] - PD.dirSpatial[1]*d1DelDirSpatial[1])/(PD.dirSpatial[2]);
}
////


// Exponential map
void ExpSO3(const double* theta, double Mat[][3])
{
  double Skw[3][3];
  HodgeStar(theta, Skw);
  double angle = sqrt(theta[0]*theta[0]+theta[1]*theta[1]+theta[2]*theta[2]);
  double sincx = gsl_sf_sinc(angle);
  double sincxby2 = gsl_sf_sinc(angle/2.);
  const double delta[3][3] = {{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};

  // I + sinc(x) Skw + (1-cos(x))/x^2 Skw^2
  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      {
	Mat[i][j] = 0.;
	for(int k=0; k<3; ++k)
	  Mat[i][j] += Skw[i][k]*Skw[k][j];
	Mat[i][j] *= 0.5*sincxby2*sincxby2;
	Mat[i][j] += (delta[i][j] + sincx*Skw[i][j]);
      }

  // Check that ExpSO3 has determinant = 1
  double det = 0.;
  for(int i=0; i<3; ++i)
    det += Mat[0][i]*(Mat[1][(i+1)%3]*Mat[2][(i+2)%3]-Mat[1][(i+2)%3]*Mat[2][(i+1)%3]);

  // cout //
  std::cout<<"\n========  Determinant:  "; 
  std::cout<<det;
  std::cout<<"\n";
  // //
  
  assert(std::abs(det-1.)<1.e-6);
  
  return;
}


// Update state
void Update(PointKinematics& PD, const VarPointKinematics& var)
{
  // Aliases
  auto* phi0 = PD.phi0;
  auto* phi1 = PD.phi1;
  auto* phi2 = PD.phi2;
  auto* d1Phi0 = PD.d1Phi0;
  auto* d1Phi1 = PD.d1Phi1;
  auto* d1Phi2 = PD.d1Phi2;
  auto* dirSpatial = PD.dirSpatial;
  auto* d1DirSpatial = PD.d1DirSpatial;
  auto* Lambda = PD.Lambda;

  const auto* eta0 = var.eta0;
  const auto* eta1 = var.eta1;
  const auto* eta2 = var.eta2;
  const auto* d1Eta0 = var.d1Eta0;
  const auto* d1Eta1 = var.d1Eta1;
  const auto* d1Eta2 = var.d1Eta2;
  const auto* delDirSpatial = var.delDirSpatial; //// ////
  //const auto* delDirMaterial = var.delDirMaterial; //// ////
  const auto* d1DelDirSpatial = var.d1DelDirSpatial;

  /*
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
  */
  
  // Update rotations
  double delTheta[3];
  delTheta[0] = dirSpatial[1]*delDirSpatial[2] - dirSpatial[2]*delDirSpatial[1];
  delTheta[1] = dirSpatial[2]*delDirSpatial[0] - dirSpatial[0]*delDirSpatial[2];
  delTheta[2] = dirSpatial[0]*delDirSpatial[1] - dirSpatial[1]*delDirSpatial[0];
  double normDelDirSpatial = sqrt(delDirSpatial[0]*delDirSpatial[0] + delDirSpatial[1]*delDirSpatial[1] + delDirSpatial[2]*delDirSpatial[2]);
  double sinc = gsl_sf_bessel_j0(normDelDirSpatial); // sin(x)/x
  double d1NormDelDirSpatial = 0.;
  for(int i=0; i<3; ++i)
    {
      d1NormDelDirSpatial += (delDirSpatial[i]*d1DelDirSpatial[i])/(normDelDirSpatial);
    }
  double exp[3][3];
  ExpSO3(delTheta, exp);
  double LambdaNew[3][3];
  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      {
	LambdaNew[i][j] = 0.;
	for(int k=0; k<3; ++k)
	  LambdaNew[i][j] += exp[i][k]*Lambda[k][j];
      }

  // Update kinematics at point
  for(int i=0; i<3; ++i)
    {
      phi0[i] += eta0[i];
      phi1[i] += eta1[i];
      phi2[i] += eta2[i];
      d1Phi0[i] += d1Eta0[i];
      d1Phi1[i] += d1Eta1[i];
      d1Phi2[i] += d1Eta2[i];
      dirSpatial[i] *= std::cos(normDelDirSpatial);
      dirSpatial[i] += sinc*delDirSpatial[i];
      d1DirSpatial[i] *= std::cos(normDelDirSpatial);
      d1DirSpatial[i] -= dirSpatial[i]*(std::sin(normDelDirSpatial))*d1NormDelDirSpatial;
      d1DirSpatial[i] += sinc*d1DelDirSpatial[i];
      d1DirSpatial[i] += delDirSpatial[i]*((normDelDirSpatial*(std::cos(normDelDirSpatial))*d1NormDelDirSpatial - (std::sin(normDelDirSpatial))*d1NormDelDirSpatial)/(normDelDirSpatial*normDelDirSpatial));

      for(int j=0; j<3; ++j)
	{
	  Lambda[i][j] = LambdaNew[i][j];
	}
    }
}

//// ////
// Update the 1st variations with 2nd variations
void UpdateVar(const PointKinematics& PD1,
	       const PointKinematics& PD2,
	       VarPointKinematics& var1,
	       const VarPointKinematics& var2,
	       const eps)
{
  for(int i=0; i<3; ++i)
    {
      var1.eta0[i] += eps*0.;
      var1.eta1[i] += eps*0.;
      var1.eta2[i] += eps*0.;
      var1.d1Eta0[i] += eps*0.;
      var1.d1Eta1[i] += eps*0.;
      var1.d1Eta2[i] += eps*0.;

      // aliases for 2nd variation calculation
      auto* t = PD1.dirSpatial;
      auto* tPrime = PD1.d1DirSpatial;
      auto* delt = var1.delDirSpatial;
      auto* deltPrime = var1.d1DelDirSpatial;
      auto* tEps = PD2.dirSpatial;
      auto* tPrimeEps = PD2.d1DirSpatial;

      double delTheta[3];
      delTheta[0] = t[1]*delt[2] - t[2]*delt[1];
      delTheta[1] = t[2]*delt[0] - t[0]*delt[2];
      delTheta[2] = t[0]*delt[1] - t[1]*delt[0];

      var1.delt[0] = delTheta[1]*tEps[2] - delTheta[2]*tEps[1];
      var1.delt[1] = delTheta[2]*tEps[0] - delTheta[0]*tEps[2];
      var1.delt[2] = delTheta[0]*tEps[1] - delTheta[1]*tEps[0];

      double tempvec1[3], tempvec2[3], tempvec3[3];

      tempvec1[0] = tPrime[1]*delt[2] - tPrime[2]*delt[1];
      tempvec1[1] = tPrime[2]*delt[0] - tPrime[0]*delt[2];
      tempvec1[2] = tPrime[0]*delt[1] - tPrime[1]*delt[0];

      tempvec2[0] = t[1]*deltPrime[2] - t[2]*deltPrime[1];
      tempvec2[1] = t[2]*deltPrime[0] - t[0]*deltPrime[2];
      tempvec2[2] = t[0]*deltPrime[1] - t[1]*deltPrime[0];

      tempvec3[0] = t[1]*delt[2] - t[2]*delt[1];
      tempvec3[1] = t[2]*delt[0] - t[0]*delt[2];
      tempvec3[2] = t[0]*delt[1] - t[1]*delt[0];

      double tempvec4[3], tempvec5[3], tempvec6[3];
      tempvec4[0] = tempvec1[1]*tEps[2] - tempvec1[2]*tEps[1];
      tempvec4[1] = tempvec1[2]*tEps[0] - tempvec1[0]*tEps[2];
      tempvec4[2] = tempvec1[0]*tEps[1] - tempvec1[1]*tEps[0];

      tempvec5[0] = tempvec2[1]*tEps[2] - tempvec2[2]*tEps[1];
      tempvec5[1] = tempvec2[2]*tEps[0] - tempvec2[0]*tEps[2];
      tempvec5[2] = tempvec2[0]*tEps[1] - tempvec2[1]*tEps[0];

      tempvec6[0] = tempvec3[1]*tPrimeEps[2] - tempvec3[2]*tPrimeEps[1];
      tempvec6[1] = tempvec3[2]*tPrimeEps[0] - tempvec3[0]*tPrimeEps[2];
      tempvec6[2] = tempvec3[0]*tPrimeEps[1] - tempvec3[1]*tPrimeEps[0];

      var1.delt[0] = tempvec4[0] + tempvec5[0] + tempvec6[0];
      var1.delt[1] = tempvec4[1] + tempvec5[1] + tempvec6[1];
      var1.delt[2] = tempvec4[2] + tempvec5[2] + tempvec6[2]; 
    }
}
//// ////



// Compute numerical values for 1st variations
void ComputeNumericalVariations(const PointKinematics& PD,
				const VarPointKinematics& var,
				const double eps,
				double varEpsilon[2][2],
				double varRho[2][2],
				double varDelta[2])
{
  // Scale variation
  VarPointKinematics varplus(var);
  VarPointKinematics varminus(var);

  for(int i=0; i<3; ++i)
    {
      varplus.eta0[i] *= eps; varminus.eta0[i] *= -eps;
      varplus.eta1[i] *= eps; varminus.eta1[i] *= -eps;
      varplus.eta2[i] *= eps; varminus.eta2[i] *= -eps;
      varplus.d1Eta0[i] *= eps; varminus.d1Eta0[i] *= -eps;
      varplus.d1Eta1[i] *= eps; varminus.d1Eta1[i] *= -eps;
      varplus.d1Eta2[i] *= eps; varminus.d1Eta2[i] *= -eps;
      //varplus.delDirSpatial[i] *= eps; varminus.delDirSpatial[i] *= -eps; ////
      varplus.d1DelDirSpatial[i] *= eps; varminus.d1DelDirSpatial[i] *= -eps; 
    }
  
  ////
  for(int i=0; i<2; ++i)
    { varplus.delDirMaterial[i] *= eps;
      varminus.delDirMaterial[i] *= -eps; }
  ////

  //// ////
  varplus.ComputeDelDirSpatial();
  varminus.ComputeDelDirSpatial();
  //// ////
  
  const double y = 0.1;
  const double w = 0.2;

  // Perturbed states
  PointKinematics PDplus(PD);
  Update(PDplus, varplus);
  PDplus.ComputeStrains(y,w);
  
  PointKinematics PDminus(PD);
  Update(PDminus, varminus);
  PDminus.ComputeStrains(y,w);

  // Numerical values of strain variations
  for(int i=0; i<2; ++i)
    {
      varDelta[i] = (PDplus.delta[i]-PDminus.delta[i])/(2.*eps);
      for(int j=0; j<2; ++j)
	{
	  varEpsilon[i][j] = (PDplus.epsilon[i][j]-PDminus.epsilon[i][j])/(2.*eps);
	  varRho[i][j] = (PDplus.rho[i][j]-PDminus.rho[i][j])/(2.*eps);
	}
    }
  return;
}


// Compute numerical values of 2nd variations
void ComputeNumerical2ndVariations(const PointKinematics& PD,
				   const VarPointKinematics& var,
				   const VarPointKinematics& VAR,
				   const double eps,
				   double vvEpsilon[2][2],
				   double vvRho[2][2],
				   double vvDelta[2])
{
  // Update PD along the variation VAR
  VarPointKinematics VARplus(VAR), VARminus(VAR);
  for(int i=0; i<3; ++i)
    {
      VARplus.eta0[i] *= eps; VARminus.eta0[i] *= -eps;
      VARplus.eta1[i] *= eps; VARminus.eta1[i] *= -eps;
      VARplus.eta2[i] *= eps; VARminus.eta2[i] *= -eps;
      VARplus.d1Eta0[i] *= eps; VARminus.d1Eta0[i] *= -eps;
      VARplus.d1Eta1[i] *= eps; VARminus.d1Eta1[i] *= -eps;
      VARplus.d1Eta2[i] *= eps; VARminus.d1Eta2[i] *= -eps;
      //VARplus.delDirSpatial[i] *= eps; VARminus.delDirSpatial[i] *= -eps; ////
      VARplus.d1DelDirSpatial[i] *= eps; VARminus.d1DelDirSpatial[i] *= -eps;
    }

  ////
  for(int i=0; i<2; ++i)
    { VARplus.delDirMaterial[i] *= eps;
      VARminus.delDirMaterial[i] *= -eps; }
  ////

  //// ////
  VARplus.ComputeDelDirSpatial();
  VARminus.ComputeDelDirSpatial();
  //// ////
  
  const double y = 0.1;
  const double w = 0.2;
  
  PointKinematics PDplus(PD);
  Update(PDplus, VARplus);
  PDplus.ComputeStrains(y,w); // for updating a1, a2 of PDplus
  
  PointKinematics PDminus(PD);
  Update(PDminus, VARminus);
  PDminus.ComputeStrains(y,w); // for updating a1, a2 of PDminus

  // Compute 1st variations with perturbed states
  VarPointKinematics varplus(var);
  varplus.PD = &PDplus;
  varplus.ComputeVarStrains(y,w);

  VarPointKinematics varminus(var);
  varminus.PD = &PDminus;
  varminus.ComputeVarStrains(y,w);

  // Compute second variations
  for(int i=0; i<2; ++i)
    {
      vvDelta[i] = (varplus.vDelta[i]-varminus.vDelta[i])/(2.*eps);
      for(int j=0; j<2; ++j)
	{
	  vvEpsilon[i][j] = (varplus.vEpsilon[i][j]-varminus.vEpsilon[i][j])/(2.*eps);
	  vvRho[i][j] = (varplus.vRho[i][j]-varminus.vRho[i][j])/(2.*eps);
	}
    }
  return;
}