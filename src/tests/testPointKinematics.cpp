#include <iostream>
#include <TriRodsPointKinematics.h>
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

// // Compute numerical values for second variations
void ComputeNumerical2ndVariations(const PointKinematics& PD,
				   const VarPointKinematics& var,
				   const VarPointKinematics& VAR,
				   const double EPS,
				   double vvEpsilon[2][2],
				   double vvRho[2][2],
				   double vvDelta[2]);

// Toggle to print results
const bool PRINT = false;

// main function
int main()
{
  // Initialize point data
  PointKinematics PD;
  Initialize(PD);
  
  const double y = 0.1;
  const double w = 0.2;

  // Compute the strains
  PD.ComputeStrains(y,w);

  // cout // // //
  /*
  std::cout<<"phi0 is : \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PD.phi0[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"phi1 is : \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PD.phi1[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"phi2 is : \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PD.phi2[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"d1Phi0 is : \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PD.d1Phi0[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"d1Phi1 is : \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PD.d1Phi1[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"d1Phi2 is : \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PD.d1Phi2[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"dirSpatial is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PD.dirSpatial[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"d1DirSpatial is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PD.d1DirSpatial[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"a1 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PD.a1[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"a2 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PD.a2[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';
  // // // //
  */
  

  // Set up two independent variations
  VarPointKinematics varA;
  varA.PD = &PD;
  for(int i=0; i<3; ++i)
    {
      varA.eta0[i] = dis(gen);
      varA.eta1[i] = dis(gen);
      varA.eta2[i] = dis(gen);
      varA.delDirSpatial[i] = dis(gen);
      varA.d1Eta0[i] = dis(gen);
      varA.d1Eta1[i] = dis(gen);
      varA.d1Eta2[i] = dis(gen);
      varA.d1DelDirSpatial[i] = dis(gen);
    }

  /*
  /////
  for(int i=0; i<2; ++i)
    varA.delDirMaterial[i] = dis(gen);
  /////
  */
  
  varA.ComputeVarStrains(y,w);
  
  /*
  // cout // // //
  
  std::cout<<"va1 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varA.va1[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"va2 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varA.va2[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';
  // // // //
  */
  

  VarPointKinematics varB;
  varB.PD = &PD;
  for(int i=0; i<3; ++i)
    {
      varB.eta0[i] = dis(gen);
      varB.eta1[i] = dis(gen);
      varB.eta2[i] = dis(gen);
      varB.delDirSpatial[i] = dis(gen);
      varB.d1Eta0[i] = dis(gen);
      varB.d1Eta1[i] = dis(gen);
      varB.d1Eta2[i] = dis(gen);
      varB.d1DelDirSpatial[i] = dis(gen);
    }

  /*
  /////
  for(int i=0; i<2; ++i)
    varB.delDirMaterial[i] = dis(gen);
  /////
  */
  
  varB.ComputeVarStrains(y,w);

  // Compute the 1st varition numerically
  const double EPS = 1.e-5;

  // Variation A
  double varAEpsilon[2][2], varARho[2][2], varADelta[2];
  ComputeNumericalVariations(PD, varA, EPS, varAEpsilon, varARho, varADelta);

  /*
  // cout // // //
  std::cout<<"varAEpsilon is: \n";
  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<2; ++j)
	{
	  std::cout<<varAEpsilon[i][j]<<'\t';
	}
      std::cout<<'\n'<<" "<<'\n';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varA.vEpsilon is: \n";
  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<2; ++j)
	{
	  std::cout<<varA.vEpsilon[i][j]<<'\t';
	}
      std::cout<<'\n'<<" "<<'\n';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varARho is: \n";
  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<2; ++j)
	{
	  std::cout<<varARho[i][j]<<'\t';
	}
      std::cout<<'\n'<<" "<<'\n';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varA.vRho is: \n";
  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<2; ++j)
	{
	  std::cout<<varA.vRho[i][j]<<'\t';
	}
      std::cout<<'\n'<<" "<<'\n';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varADelta is: \n";
  for(int i=0; i<2; ++i)
    {
      std::cout<<varADelta[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varA.delta is: \n";
  for(int i=0; i<2; ++i)
    {
      std::cout<<varA.vDelta[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  // // // //
  */
  
  
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

  /*
  // cout // // //

  std::cout<<"varBEpsilon is: \n";
  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<2; ++j)
	{
	  std::cout<<varBEpsilon[i][j]<<'\t';
	}
      std::cout<<'\n'<<" "<<'\n';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varB.vEpsilon is: \n";
  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<2; ++j)
	{
	  std::cout<<varB.vEpsilon[i][j]<<'\t';
	}
      std::cout<<'\n'<<" "<<'\n';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varBRho is: \n";
  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<2; ++j)
	{
	  std::cout<<varBRho[i][j]<<'\t';
	}
      std::cout<<'\n'<<" "<<'\n';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varB.vRho is: \n";
  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<2; ++j)
	{
	  std::cout<<varB.vRho[i][j]<<'\t';
	}
      std::cout<<'\n'<<" "<<'\n';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varBDelta is: \n";
  for(int i=0; i<2; ++i)
    {
      std::cout<<varBDelta[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varB.delta is: \n";
  for(int i=0; i<2; ++i)
    {
      std::cout<<varB.vDelta[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  // // // //
  */
  
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
      PD.d1DirSpatial[i] = dis(gen);
    }
  double normDirSpatial = sqrt(PD.dirSpatial[0]*PD.dirSpatial[0] + PD.dirSpatial[1]*PD.dirSpatial[1] + PD.dirSpatial[2]*PD.dirSpatial[2]);

  for(int i=0; i<3; ++i)
    PD.dirSpatial[i] /= normDirSpatial;
  
  /*for(int i=0; i<2; ++i)
    {
      PD.q[i] = 0.;
      PD.delta[i] = 0.;
      for(int j=0; j<2; ++j)
	{
	  PD.n[i][j] = 0.;
	  PD.m[i][j] =0.;
	  PD.epsilon[i][j] = 0.;
	  PD.rho[i][j] = 0.;
	}
	}*/
  
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
  
  /*while(true)
    {
    double norm = 0.;
    for(int i=0; i<3; ++i)
    {
    e[i] = dis(gen);
    norm +=e[i]*e[i];
    }
    norm = sqrt(norm);
    if(norm>1.e-2)
    {
    for(int i=0; i<3; ++i)
    e[i] /= norm;
    break;
    }
    }*/
  
  // Create the rotation matrix
  double StarE[3][3];
  HodgeStar(e, StarE);
  double delta[3][3] = {{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      PD.Lambda[i][j] =
	e[i]*e[j] + (delta[i][j]-e[i]*e[j])*(std::cos(theta)) + StarE[i][j]*(std::sin(theta));
  
  // Check that the rotation matrix has determinant 1
  double det = 0.;
  for(int i=0; i<3; ++i)
    det += PD.Lambda[0][i]*
      (PD.Lambda[1][(i+1)%3]*PD.Lambda[2][(i+2)%3]-PD.Lambda[1][(i+2)%3]*PD.Lambda[2][(i+1)%3]);

  /*
  // cout //
  std::cout<<"\n========";
  std::cout<<"initialization of PD part is being executed \n";
  */
  
  assert(std::abs(det-1.)<1.e-6);
}


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
      varplus.delDirSpatial[i] *= eps; varminus.delDirSpatial[i] *= -eps;
      varplus.d1DelDirSpatial[i] *= eps; varminus.d1DelDirSpatial[i] *= -eps; //
    }
  const double y = 0.1;
  const double w = 0.2;

  // Perturbed states
  PointKinematics PDplus(PD);
  Update(PDplus, varplus);
  PDplus.ComputeStrains(y,w);

  /*
  // cout // // //

  std::cout<<"PDplus.phi0 is : \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PDplus.phi0[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDplus.phi1 is : \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PDplus.phi1[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDplus.phi2 is : \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PDplus.phi2[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDplus.d1Phi0 is : \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PDplus.d1Phi0[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDplus.d1Phi1 is : \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PDplus.d1Phi1[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDplus.d1Phi2 is : \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PDplus.d1Phi2[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDplus.dirSpatial is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PDplus.dirSpatial[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDplus.d1DirSpatial is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PDplus.d1DirSpatial[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDplus.a1 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PDplus.a1[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDplus.a2 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PDplus.a2[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varplus.eta0 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varplus.eta0[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varplus.eta1 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varplus.eta1[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varplus.eta2 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varplus.eta2[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varplus.d1Eta0 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varplus.d1Eta0[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varplus.d1Eta1 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varplus.d1Eta1[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varplus.d1Eta2 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varplus.d1Eta2[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  // // // //
  */

  
  PointKinematics PDminus(PD);
  Update(PDminus, varminus);
  PDminus.ComputeStrains(y,w);

  /*
  // cout // // //

  std::cout<<"PDminus.phi0 is : \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PDminus.phi0[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDminus.phi1 is : \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PDminus.phi1[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDminus.phi2 is : \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PDminus.phi2[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDminus.d1Phi0 is : \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PDminus.d1Phi0[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDminus.d1Phi1 is : \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PDminus.d1Phi1[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDminus.d1Phi2 is : \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PDminus.d1Phi2[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDminus.dirSpatial is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PDminus.dirSpatial[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDminus.d1DirSpatial is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PDminus.d1DirSpatial[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDminus.a1 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PDminus.a1[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDminus.a2 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<PDminus.a2[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varplus.eta0 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varplus.eta0[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varplus.eta1 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varplus.eta1[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varminus.eta2 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varminus.eta2[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varminus.d1Eta0 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varminus.d1Eta0[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varminus.d1Eta1 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varminus.d1Eta1[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"varminus.d1Eta2 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varminus.d1Eta2[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  // // // //


  // cout // // //

  std::cout<<"PDplus.epsilon is: \n";
  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<2; ++j)
	{
	  std::cout<<PDplus.epsilon[i][j]<<'\t';
	}
      std::cout<<'\n'<<" "<<'\n';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDminus.epsilon is: \n";
  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<2; ++j)
	{
	  std::cout<<PDminus.epsilon[i][j]<<'\t';
	}
      std::cout<<'\n'<<" "<<'\n';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDplus.rho is: \n";
  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<2; ++j)
	{
	  std::cout<<PDplus.rho[i][j]<<'\t';
	}
      std::cout<<'\n'<<" "<<'\n';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDminus.rho is: \n";
  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<2; ++j)
	{
	  std::cout<<PDminus.rho[i][j]<<'\t';
	}
      std::cout<<'\n'<<" "<<'\n';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDplus.delta is: \n";
  for(int i=0; i<2; ++i)
    {
      std::cout<<PDplus.delta[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  std::cout<<"PDminus.delta is: \n";
  for(int i=0; i<2; ++i)
    {
      std::cout<<PDminus.delta[i]<<'\t';
    }
  std::cout<<'\n'<<" "<<'\n';

  // // // //
  */

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

  // cout //
  std::cout<<"\n=======";
  std::cout<<"PDplus & PDminus are updated \n";
  // //
  
  // Compute 1st variations with perturbed states
  VarPointKinematics varplus(var);
  varplus.PD = &PDplus;
  varplus.ComputeVarStrains(y,w);

  // cout //
  std::cout<<"\n=======";
  std::cout<<"varplus.PD->t is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varplus.PD->dirSpatial[i]<<'\t';
    }
  std::cout<<'\n'<<"  "<<'\n';
  // //

  // cout //
  std::cout<<"\n=======";
  std::cout<<"varplus.PD->t' is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varplus.PD->d1DirSpatial[i]<<'\t';
    }
  std::cout<<'\n'<<"  "<<'\n';
  // //

  // cout //
  std::cout<<"\n=======";
  std::cout<<"varplus.PD->t.varplus->t' is: \n";
  double tempvec = 0.;
  for(int i=0; i<3; ++i)
    {
      tempvec += (varplus.PD->dirSpatial[i])*(varplus.PD->d1DirSpatial[i]);
    }
  std::cout<<tempvec<<'\n'<<"  "<<'\n';
  // //

  // cout //
  std::cout<<"\n=======";
  std::cout<<"varplus.PD->a1 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varplus.PD->a1[i]<<'\t';
    }
  std::cout<<'\n'<<"  "<<'\n';
  // //

  // cout //
  std::cout<<"\n=======";
  std::cout<<"varplus.PD->a2 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varplus.PD->a2[i]<<'\t';
    }
  std::cout<<'\n'<<"  "<<'\n';
  // //

  // cout //
  std::cout<<"\n=======";
  std::cout<<"varplus.va1 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varplus.va1[i]<<'\t';
    }
  std::cout<<'\n'<<"  "<<'\n';
  // //

  // cout //
  std::cout<<"\n=======";
  std::cout<<"varplus.va2 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varplus.va2[i]<<'\t';
    }
  std::cout<<'\n'<<"  "<<'\n';
  // //

  // cout //
  std::cout<<"\n=======";
  std::cout<<"varplus.vt is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varplus.delDirSpatial[i]<<'\t';
    }
  std::cout<<'\n'<<"  "<<'\n';
  // //

  // cout //
  std::cout<<"\n=======";
  std::cout<<"varplus.vt' is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varplus.d1DelDirSpatial[i]<<'\t';
    }
  std::cout<<'\n'<<"  "<<'\n';
  // //

  // cout //
  std::cout<<"\n=======";
  std::cout<<"varplus.vEpsilon is: \n";
  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<2; ++j)
	{
	  std::cout<<varplus.vEpsilon[i][j]<<'\t';
	}
      std::cout<<'\n'<<"  "<<'\n';
    }
  // //

  
  VarPointKinematics varminus(var);
  varminus.PD = &PDminus;
  varminus.ComputeVarStrains(y,w);

  // cout //
  std::cout<<"\n=======";
  std::cout<<"varminus.PD->t is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varminus.PD->dirSpatial[i]<<'\t';
    }
  std::cout<<'\n'<<"  "<<'\n';
  // //

  // cout //
  std::cout<<"\n=======";
  std::cout<<"varminus.PD->t' is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varminus.PD->d1DirSpatial[i]<<'\t';
    }
  std::cout<<'\n'<<"  "<<'\n';
  // //

  // cout //
  std::cout<<"\n=======";
  std::cout<<"varminus.PD->a1 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varminus.PD->a1[i]<<'\t';
    }
  std::cout<<'\n'<<"  "<<'\n';
  // //

  // cout //
  std::cout<<"\n=======";
  std::cout<<"varminus.PD->a2 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varminus.PD->a2[i]<<'\t';
    }
  std::cout<<'\n'<<"  "<<'\n';
  // //

  // cout //
  std::cout<<"\n=======";
  std::cout<<"varminus.va1 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varminus.va1[i]<<'\t';
    }
  std::cout<<'\n'<<"  "<<'\n';
  // //

  // cout //
  std::cout<<"\n=======";
  std::cout<<"varminus.va2 is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varminus.va2[i]<<'\t';
    }
  std::cout<<'\n'<<"  "<<'\n';
  // //

  // cout //
  std::cout<<"\n=======";
  std::cout<<"varminus.vt is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varminus.delDirSpatial[i]<<'\t';
    }
  std::cout<<'\n'<<"  "<<'\n';
  // //

  // cout //
  std::cout<<"\n=======";
  std::cout<<"varminus.vt' is: \n";
  for(int i=0; i<3; ++i)
    {
      std::cout<<varminus.d1DelDirSpatial[i]<<'\t';
    }
  std::cout<<'\n'<<"  "<<'\n';
  // //

  // cout //
  std::cout<<"\n=======";
  std::cout<<"varminus.vEpsilon is: \n";
  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<2; ++j)
	{
	  std::cout<<varminus.vEpsilon[i][j]<<'\t';
	}
      std::cout<<'\n'<<"  "<<'\n';
    }
  // //

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
