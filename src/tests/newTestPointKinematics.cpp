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
				   double vvDelta[2],
				   double delDelt[3],
				   double delDeltPrime[3]);

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

  /*
//////////
std::cout<<'\n'<<'\n'<<"for PD=========="<<'\n';
std::cout<<std::abs(PD.dirSpatial[0]*PD.d1DirSpatial[0] + PD.dirSpatial[1]*PD.d1DirSpatial[1] + PD.dirSpatial[2]*PD.d1DirSpatial[2])<<'\n';
assert(std::abs(PD.dirSpatial[0]*PD.d1DirSpatial[0] + PD.dirSpatial[1]*PD.d1DirSpatial[1] + PD.dirSpatial[2]*PD.d1DirSpatial[2]) < 1.e-10);
/////////
*/
  
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
      varA.d1Eta0[i] = dis(gen);
      varA.d1Eta1[i] = dis(gen);
      varA.d1Eta2[i] = dis(gen);
    }
  for(int i=0; i<2; ++i)
    varA.delDirMaterial[i] = dis(gen);
  varA.ComputeDelDirSpatial();
  double kA = 0.;
  GenerateRandomD1DelDirSpatial(PD, varA, kA, varA.d1DelDirSpatial);

  /*
//////////
std::cout<<'\n'<<'\n'<<"for varA=========="<<'\n';
std::cout<<std::abs(kA+(PD.dirSpatial[0]*varA.d1DelDirSpatial[0]+PD.dirSpatial[1]*varA.d1DelDirSpatial[1]+PD.dirSpatial[2]*varA.d1DelDirSpatial[2]))<<'\n'<<'\n';
assert(std::abs(kA+(PD.dirSpatial[0]*varA.d1DelDirSpatial[0]+PD.dirSpatial[1]*varA.d1DelDirSpatial[1]+PD.dirSpatial[2]*varA.d1DelDirSpatial[2]))<1.e-10);
/////////
*/
   
  varA.ComputeVarStrains(y,w);
 
  VarPointKinematics varB;
  varB.PD = &PD;
  for(int i=0; i<3; ++i)
    {
      varB.eta0[i] = dis(gen);
      varB.eta1[i] = dis(gen);
      varB.eta2[i] = dis(gen);
      varB.d1Eta0[i] = dis(gen);
      varB.d1Eta1[i] = dis(gen);
      varB.d1Eta2[i] = dis(gen);
    }
  for(int i=0; i<2; ++i)
    varB.delDirMaterial[i] = dis(gen);
  varB.ComputeDelDirSpatial();
  double kB = 0.;
  GenerateRandomD1DelDirSpatial(PD, varB, kB, varB.d1DelDirSpatial);

  /*
//////////
std::cout<<'\n'<<'\n'<<"for varB============="<<'\n';
std::cout<<std::abs(kB+(PD.dirSpatial[0]*varB.d1DelDirSpatial[0]+PD.dirSpatial[1]*varB.d1DelDirSpatial[1]+PD.dirSpatial[2]*varB.d1DelDirSpatial[2]))<<'\n'<<'\n';
assert(std::abs(kB+(PD.dirSpatial[0]*varB.d1DelDirSpatial[0]+PD.dirSpatial[1]*varB.d1DelDirSpatial[1]+PD.dirSpatial[2]*varB.d1DelDirSpatial[2]))<1.e-10);
/////////
*/
  
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

  double delDELt[3], delDELtPrime[3];
  double temp10 = 0.;
  double temp11 = 0.;
  double temp12 = 0.;
  for(int i=0; i<3; ++i)
    {
      temp10 += varA.delDirSpatial[i]*varB.delDirSpatial[i];
      temp11 += varA.d1DelDirSpatial[i]*varB.delDirSpatial[i];
      temp12 += varA.delDirSpatial[i]*varB.d1DelDirSpatial[i];
    }
  for(int i=0; i<3; ++i)
    {
      delDELt[i] = -temp10*PD.dirSpatial[i];
      delDELtPrime[i] = -temp11*PD.dirSpatial[i] - temp12*PD.dirSpatial[i] - temp10*PD.d1DirSpatial[i]; // 
    }

  /*
////
double temp20 = 0.;
double temp21 = 0.;
double temp22 = 0.;

for(int i=0; i<3; ++i)
{
temp20 += PD.d1DirSpatial[i]*varB.delDirSpatial[i];
temp21 += PD.dirSpatial[i]*varB.d1DelDirSpatial[i];
temp22 += PD.dirSpatial[i]*varB.delDirSpatial[i];
}
std::cout<<'\n'<<'\n'<<"temp20 is: "<<temp20<<'\n'<<"temp21 is: "<<temp21<<'\n'<<"temp20+temp21 is: "<<(temp20+temp21)<<'\n'<<"temp22 is: "<<temp22<<'\n'<<'\n';
assert((temp20+temp21)<1e-10);
////
*/

  /*
////
std::cout<<"\n\nPD.t is: ";
for(int i=0; i<3; ++i)
{
std::cout<<PD.dirSpatial[i]<<'\t';
}
std::cout<<"\n\n";

std::cout<<"\n\nPD.t' is: ";
for(int i=0; i<3; ++i)
{
std::cout<<PD.d1DirSpatial[i]<<'\t';
}
std::cout<<"\n\n";

std::cout<<"\n\nvarA.delt is: ";
for(int i=0; i<3; ++i)
{
std::cout<<varA.delDirSpatial[i]<<'\t';
}
std::cout<<"\n\n";

std::cout<<"\n\nvarA.delt' is: ";
for(int i=0; i<3; ++i)
{
std::cout<<varA.d1DelDirSpatial[i]<<'\t';
}
std::cout<<"\n\n";

std::cout<<"\n\nvarB.delt is: ";
for(int i=0; i<3; ++i)
{
std::cout<<varB.delDirSpatial[i]<<'\t';
}
std::cout<<"\n\n";

std::cout<<"\n\nvarB.delt' is: ";
for(int i=0; i<3; ++i)
{
std::cout<<varB.d1DelDirSpatial[i]<<'\t';
}
std::cout<<"\n\n";
////

////
double a[3], b[3], c[3];
double d=0.;
double e=0.;
double f=0.;
for(int i=0; i<3; ++i)
{
d += varA.delDirSpatial[i]*varB.delDirSpatial[i];
e += varA.d1DelDirSpatial[i]*varB.delDirSpatial[i];
f += varA.delDirSpatial[i]*varB.d1DelDirSpatial[i]; 
}
for(int i=0; i<3; ++i)
{
a[i] = -PD.d1DirSpatial[i]*d;
b[i] = -PD.dirSpatial[i]*e;
c[i] = -PD.dirSpatial[i]*f;
}

std::cout<<"\n\na is: \n";
for(int i=0; i<3; ++i)
std::cout<<a[i]<<'\t';
std::cout<<"\n\n";

std::cout<<"\n\nb is: \n";
for(int i=0; i<3; ++i)
std::cout<<b[i]<<'\t';
std::cout<<"\n\n";

std::cout<<"\n\nc is: \n";
for(int i=0; i<3; ++i)
std::cout<<c[i]<<'\t';
std::cout<<"\n\n";
////
*/

  // Check consistency of 2nd variations
  double vvEpsilon[2][2], vvRho[2][2], vvDelta[2];
  double delDELtNum[3], delDELtPrimeNum[3];
  ComputeNumerical2ndVariations(PD, varA, varB, EPS, vvEpsilon, vvRho, vvDelta, delDELtNum, delDELtPrimeNum); 

  if(PRINT)
    {
      std::cout<<"\n2nd variation of delt : ";
      for(int i=0; i<3; ++i)
	std::cout<<"\n"<<delDELt[i]<<" should be "<<delDELtNum[i];
      std::cout<<'\n'<<'\n';
      std::cout<<"\n2nd variation of delt' : ";
      for(int i=0; i<3; ++i)
	std::cout<<"\n"<<delDELtPrime[i]<<" should be "<<delDELtPrimeNum[i];
      std::cout<<"\n\n"; std::fflush( stdout );
    }
  for(int i=0; i<3; ++i)
    {
      assert((std::abs(delDELt[i]-delDELtNum[i]))<1.e-4);
      assert((std::abs(delDELtPrime[i]-delDELtPrimeNum[i]))<1.e-4);
    }


  if(PRINT) // true
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
    }
  double normDirSpatial = sqrt(PD.dirSpatial[0]*PD.dirSpatial[0] + PD.dirSpatial[1]*PD.dirSpatial[1] + PD.dirSpatial[2]*PD.dirSpatial[2]);

  for(int i=0; i<3; ++i)
    PD.dirSpatial[i] /= normDirSpatial;
  
  double arbvec[3];
  for(int i=0; i<3; ++i)
    arbvec[i] = dis(gen);
  PD.d1DirSpatial[0] = PD.dirSpatial[1]*arbvec[2] - PD.dirSpatial[2]*arbvec[1];
  PD.d1DirSpatial[1] = PD.dirSpatial[2]*arbvec[0] - PD.dirSpatial[0]*arbvec[2];
  PD.d1DirSpatial[2] = PD.dirSpatial[0]*arbvec[1] - PD.dirSpatial[1]*arbvec[0];
  
  PD.ComputeRotationMatrix(); 
 
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

  /*
  // cout //
  std::cout<<"\n========  Determinant:  "; 
  std::cout<<det;
  std::cout<<"\n";
  // //
  */
  
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
  const auto* delDirSpatial = var.delDirSpatial; 
  const auto* d1DelDirSpatial = var.d1DelDirSpatial;
  
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

// Compute cross product of three vectors
void CrossProduct3(const double a[3],
		   const double b[3],
		   const double c[3],
		   double d[3])
{
  double k1=0.;
  double k2=0.;
  for(int i=0; i<3; ++i)
    {
      k1 += a[i]*c[i];
      k2 += b[i]*c[i];
    }
  for(int i=0; i<3; ++i)
    d[i] = k1*b[i] - k2*a[i];
}

// Update the 1st variations with 2nd variations
void UpdateVar(VarPointKinematics& var1,
	       const PointKinematics& PD1,
	       const PointKinematics& PD2)
{
  for(int i=0; i<3; ++i)
    {
      //var1.eta0[i] += eps*0.;
      //var1.eta1[i] += eps*0.;
      //var1.eta2[i] += eps*0.;
      //var1.d1Eta0[i] += eps*0.;
      //var1.d1Eta1[i] += eps*0.;
      //var1.d1Eta2[i] += eps*0.;

      // aliases for 2nd variation calculation
      const auto* t = PD1.dirSpatial;
      const auto* tPrime = PD1.d1DirSpatial;
      auto* delt = var1.delDirSpatial;
      auto* deltPrime = var1.d1DelDirSpatial;
      const auto* tEps = PD2.dirSpatial;
      const auto* tPrimeEps = PD2.d1DirSpatial;

      double tempvec1[3], tempvec2[3], tempvec3[3], tempvec4[3];
      
      CrossProduct3(t, delt, tEps, tempvec1);
      CrossProduct3(tPrime, delt, tEps, tempvec2);
      CrossProduct3(t, deltPrime, tEps, tempvec3);
      CrossProduct3(t, delt, tPrimeEps, tempvec4);

      for(int i=0; i<3; ++i)
	{
	  delt[i] = tempvec1[i];
	  deltPrime[i] = tempvec2[i] + tempvec3[i] + tempvec4[i];
	}
      ////

      /*
      double delTheta[3];
      delTheta[0] = t[1]*delt[2] - t[2]*delt[1];
      delTheta[1] = t[2]*delt[0] - t[0]*delt[2];
      delTheta[2] = t[0]*delt[1] - t[1]*delt[0];

      delt[0] = delTheta[1]*tEps[2] - delTheta[2]*tEps[1];
      delt[1] = delTheta[2]*tEps[0] - delTheta[0]*tEps[2];
      delt[2] = delTheta[0]*tEps[1] - delTheta[1]*tEps[0];

      // delt is being UPDATED, now we should not use that updated value, but we were using those updated values of delt, that's why it was not working

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

      for(int i=0; i<3; ++i)
	deltPrime[i] = tempvec4[i] + tempvec5[i] + tempvec6[i];
      */
    }
}

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
      varplus.d1DelDirSpatial[i] *= eps; varminus.d1DelDirSpatial[i] *= -eps; 
    }
  for(int i=0; i<2; ++i)
    { varplus.delDirMaterial[i] *= eps;
      varminus.delDirMaterial[i] *= -eps; }

  varplus.ComputeDelDirSpatial();
  varminus.ComputeDelDirSpatial();

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
				   double vvDelta[2],
				   double delDelt[3],      
				   double delDeltPrime[3]) 
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
      VARplus.delDirSpatial[i] *= eps; VARminus.delDirSpatial[i] *= -eps; 
      VARplus.d1DelDirSpatial[i] *= eps; VARminus.d1DelDirSpatial[i] *= -eps;
    }

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

  /*
  double tempvec20[3], tempvec21[3], tempvec22[3], tempvec23[3];
  CrossProduct3(PD.d1DirSpatial, varplus.delDirSpatial, PDplus.dirSpatial, tempvec20);
  CrossProduct3(PD.dirSpatial, varplus.d1DelDirSpatial, PDplus.dirSpatial, tempvec21);
  CrossProduct3(PD.dirSpatial, varplus.delDirSpatial, PDplus.d1DirSpatial, tempvec22);
  for(int i=0; i<3; ++i)
    tempvec23[i] = tempvec20[i] + tempvec21[i] + tempvec22[i];
  */

  varplus.PD = &PDplus;
  UpdateVar(varplus, PD, PDplus);
  varplus.ComputeVarStrains(y,w);

  VarPointKinematics varminus(var);

  /*
  double tempvec30[3], tempvec31[3], tempvec32[3], tempvec33[3];
  CrossProduct3(PD.d1DirSpatial, varminus.delDirSpatial, PDplus.dirSpatial, tempvec30);
  CrossProduct3(PD.dirSpatial, varminus.d1DelDirSpatial, PDplus.dirSpatial, tempvec31);
  CrossProduct3(PD.dirSpatial, varminus.delDirSpatial, PDplus.d1DirSpatial, tempvec32);
  for(int i=0; i<3; ++i)
    tempvec33[i] = tempvec30[i] + tempvec31[i] + tempvec32[i];
  */
  
  varminus.PD = &PDminus;
  UpdateVar(varminus, PD, PDminus);
  varminus.ComputeVarStrains(y,w);

  /*
  ////
  std::cout<<"\n\nCalculated varplus.d1DirSpatial is: \n";
  for(int i=0; i<3; ++i)
    std::cout<<tempvec23[i]<<'\t';
  std::cout<<"\n\n";
  
  std::cout<<"\n\nvarplus.d1DelDirSpatial is: \n";
  for(int i=0; i<3; ++i)
    std::cout<<varplus.d1DelDirSpatial[i]<<'\t';
  std::cout<<"\n\n";

  std::cout<<"\n\nCalculated varminus.d1DirSpatial is: \n";
  for(int i=0; i<3; ++i)
    std::cout<<tempvec33[i]<<'\t';
  std::cout<<"\n\n";
  
  std::cout<<"\n\nvarminus.d1DelDirSpatial is: \n";
  for(int i=0; i<3; ++i)
    std::cout<<varminus.d1DelDirSpatial[i]<<'\t';
  std::cout<<"\n\n";
  ////
  */

  ////
  for(int i=0; i<3; ++i)
    {
      delDelt[i] = (varplus.delDirSpatial[i] - varminus.delDirSpatial[i])/(2*eps);
      delDeltPrime[i] = (varplus.d1DelDirSpatial[i] - varminus.d1DelDirSpatial[i])/(2*eps);
    }
  
  //varminus.ComputeVarStrains(y,w);

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
