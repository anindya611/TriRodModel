#include <TriRodsOps.h>
#include <TriRodsPointKinematicsNew.h>
#include <TriRodsUtils.h>
#include <cmath>
#include <cassert>
#include <iostream>

namespace TriRods
{
  // Computes the energy functional
  double Ops::GetEnergy(const void* arg) const
  {
    // Cast configuration
    const TriRodsConfiguration* Config = static_cast<const TriRodsConfiguration*>(arg);
    assert(Config!=nullptr && "TriRods::Ops::GetEnergy() - Invalid pointer to TriRodsConfiguration");

    // Check that this configuration is "set"
    assert(Config->IsInitialized() && "Ops::GetEnergy()- TriRodsConfiguration not initialized");

    // Access workspace for quadrature point-level calculatons
    auto& PD = WrkSpc->PD;
    
    // Number of quadrature points
    const auto& nQuadsX = ConfigStr.nQuadsPerElmX;
    assert(nQuadsX==static_cast<int>(Elm->GetIntegrationWeights(fields[0]).size()));

    // Node number for this element
    const auto& nodesX = ConfigStr.nodesX;
    assert(static_cast<int>(nodesX.size())==GetFieldDof(fields[0]));
    const int ndofs = static_cast<int>(nodesX.size());
    
    // This element number
    const auto& e = ConfigStr.elmX;

    // Quadrature weights along x direction
    const auto& QwtsX = Elm->GetIntegrationWeights(fields[0]);

    // Return the properties of the Ribbon
    MatProperties props = Mat->GetMaterialProperties();

    // Return the width of the ribbon
    double w = props.GetWidth();

    // Define the 5 point quadrature point coordinates and the corresponding weights
    const double qPtsY[5] = {-0.90618,
			     -0.538469,
			      0.,
			      0.538469,
			      0.90618}; // Quadrature points from 5 point gauss quadrature rule, they lie in between -1 and 1
    
    const double qWtsY[5] = {0.236927,
			     0.478629,
			     0.568889,
			     0.478629,
			     0.236927}; // Quadratute weights

    double qPtCoordinatesY[5];
    for(int i=0; i<5; ++i)
      qPtCoordinatesY[i] = (w*qPtsY[i])/2.;
    
    double qPtWtsY[5];
    for(int i=0; i<5; ++i)
      qPtWtsY[i] = (w*qWtsY[i])/2.;

    ////
    // Rotation matrix array
    std::vector<Mat3> Lambda1;
    for(int a=0; a<ndofs; ++a)
      Lambda1[a] = Config->GetRotations(nodesX[a]);
    ////

    // Energy to be computed
    double Energy = 0.;
    
    // Integrate
    for(int qx=0; qx<nQuadsX; ++qx)
      {
	// Centerline, edge_1 and edge_2 coordinates at a quadrature point
	for(int f=0; f<3; ++f)
	  {
	    PD.phi0[f] = 0.;
	    PD.phi1[f] = 0.;
	    PD.phi2[f] = 0.;
	    PD.d1Phi0[f] = 0.;
	    PD.d1Phi1[f] = 0.;
	    PD.d1Phi2[f] = 0.;
	    for(int a=0; a<ndofs; ++a)
	      {
		const auto& x0 = Config->GetLine0(nodesX[a]);
		const auto& x1 = Config->GetLine1(nodesX[a]);
		const auto& x2 = Config->GetLine2(nodesX[a]);
		PD.phi0[f] += x0[f]*Elm->GetShape(fields[f],qx,a);
		PD.phi1[f] += x1[f]*Elm->GetShape(fields[f],qx,a);
		PD.phi2[f] += x2[f]*Elm->GetShape(fields[f],qx,a);
		PD.d1Phi0[f] += x0[f]*Elm->GetDShape(fields[f],qx,a,0);
		PD.d1Phi1[f] += x1[f]*Elm->GetDShape(fields[f],qx,a,0);
		PD.d1Phi2[f] += x2[f]*Elm->GetDShape(fields[f],qx,a,0);
	      }
	  }

	// Get Spatial_Director at a quadrature point
	const auto& dirSpat = Config->GetDirSpatial(e,qx);
	for(int i=0; i<3; ++i)
	  PD.dirSpatial[i] = dirSpat[i];

	// Derivative of Spatial_director with respect to xi1 at a quadrature point
	const auto& d1DirSpat = Config->GetD1DirSpatial(e,qx); ////
	for(int i=0; i<3; ++i)
	  PD.d1DirSpatial[i] = d1DirSpat[i];

	double EnergyY = 0.;
	for(int i=0; i<5; ++i)
	  {
	    PD.ComputeStrains(qPtCoordinatesY[i],w);

	    // Update the energy at a particular x
	    EnergyY += qPtWtsY[i]*Mat->ComputeStrainEnergyDensity(PD.epsilon, PD.rho, PD.delta);
	  }

	// Update the total energy
	Energy += QwtsX[qx]*EnergyY;
      }

    return Energy;
  }




  // Residual calculation
  void Ops::GetVal(const void* arg,
		   std::vector<std::vector<double>>* resval) const
  { GetDVal(arg, resval, nullptr); }

  // Residual and stiffness calculation
  void Ops::GetDVal(const void* arg,
		    std::vector<std::vector<double>>* resval,
		    std::vector<std::vector<std::vector<std::vector<double>>>>* dresval) const
  {
    // Cast configuration variable
    const TriRodsConfiguration* Config = static_cast<const TriRodsConfiguration*>(arg);
    assert(Config!=nullptr && "TriRods::Ops::GetDVal- Invalid configuration pointer");

    // Check that this provided configuration has been set
    assert(Config->IsInitialized() &&
	   "TriRods::Ops::GetDResidue()- TriRodsConfiguration not initialized.");

    // Initialize residual and stiffness to zero
    SetZero(resval, dresval);

    // Assume that at least one of the residual/stiffness is required
    assert((resval!=nullptr || dresval!=nullptr) &&
	   "TriRods::Ops::GetDResidue()- Residual is assumed to be required");

    // Access the workspace 
    auto& PD = WrkSpc->PD;
    auto& varPD = WrkSpc->varPD;
    auto& VARvar = WrkSpc->varvarPD;
    //auto& tempval = WrkSpc->tempval;

    // Number of quadrature points
    const auto& nQuadsX = ConfigStr.nQuadsPerElmX;

    // Quadrature weights for x direction
    const auto& QwtsX = Elm->GetIntegrationWeights(fields[0]);
    assert(static_cast<int>(QwtsX.size())==nQuadsX);

    // Access details for configuration
    const auto& e = ConfigStr.elmX;
    const auto& nodesX = ConfigStr.nodesX;
    const int ndofs = static_cast<int>(nodesX.size());
    assert(ndofs==GetFieldDof(fields[0]) && "TriRods::Ops::GetDResidue- Inconsistent number of dofs");

    // Return the properties of the Ribbon
    MatProperties props = Mat->GetMaterialProperties();

    // Return the width of the ribbon
    double w = props.GetWidth();

    // Define the 5 point quadrature point coordinates and the corresponding weights
    const double qPtsY[5] = {-0.90618,
			     -0.538469,
			      0.,
			      0.538469,
			      0.90618}; // Quadrature points from 5 point gauss quadrature rule, they lie in between -1 and 1
    
    const double qWtsY[5] = {0.236927,
			     0.478629,
			     0.568889,
			     0.478629,
			     0.236927}; // Quadratute weights

    double qPtCoordinatesY[5];
    for(int i=0; i<5; ++i)
      qPtCoordinatesY[i] = (w*qPtsY[i])/2.;
    double qPtWtsY[5];
    for(int i=0; i<5; ++i)
      qPtWtsY[i] = (w*qWtsY[i])/2.;

    ////
    // Rotation matrix array
    std::vector<Mat3> Lambda1;
    for(int a=0; a<ndofs; ++a)
      Lambda1[a] = Config->GetRotations(nodesX[a]);
    ////
    
    // Integrate
    for(int qx=0; qx<nQuadsX; ++qx)
      {
	for(int z=0; z<5; ++z)
	  {
	    // Set the configuration(centerline, edge_1, edge_2) at this point
	    for(int f=0; f<3; ++f)
	      {
		PD.phi0[f] = 0.;
		PD.phi1[f] = 0.;
		PD.phi2[f] = 0.;
		PD.d1Phi0[f] = 0.;
		PD.d1Phi1[f] = 0.;
		PD.d1Phi2[f] = 0.;
		for(int a=0; a<ndofs; ++a)
		  {
		    const auto& x0 = Config->GetLine0(nodesX[a]);
		    const auto& x1 = Config->GetLine1(nodesX[a]);
		    const auto& x2 = Config->GetLine2(nodesX[a]);
		    PD.phi0[f] += x0[f]*Elm->GetShape(fields[f],qx,a);
		    PD.phi1[f] += x1[f]*Elm->GetShape(fields[f+3],qx,a);
		    PD.phi2[f] += x2[f]*Elm->GetShape(fields[f+6],qx,a);
		    PD.d1Phi0[f] += x0[f]*Elm->GetDShape(fields[f],qx,a,0);
		    PD.d1Phi1[f] += x1[f]*Elm->GetDShape(fields[f+3],qx,a,0);
		    PD.d1Phi2[f] += x2[f]*Elm->GetDShape(fields[f+6],qx,a,0);
		  }
	      }

	    // Rotation and Spatial_Director at a quadrature point
	    const auto& dirSpat = Config->GetDirSpatial(e,qx);
	    for(int i=0; i<3; ++i)
	      PD.dirSpatial[i] = dirSpat[i];

	    // Derivative of Spatial_Director at a quadrature point
	    const auto& d1DirSpat = Config->GetD1DirSpatial(e,qx);
	    for(int i=0; i<3; ++i)
	      PD.d1DirSpatial[i] = d1DirSpat[i];

	    /*
	    // Get the rotation matrix
	    const auto& Lambda = Config->GetRotations(e,qx);
	    for(int i=0; i<3; ++i)
	      for(int j=0; j<3; ++j)
		PD.Lambda[i][j] = Lambda[i][j];
	    */

	    ////// Compute strains of PD, then compute sress-resultants
	    PD.ComputeStrains(qPtCoordinatesY[z],w);
	    Mat->ComputeStressResultants(PD.epsilon, PD.rho, PD.delta,PD.n, PD.m, PD.q);

	    // Set the variations of centerline, edge_1, edge_2
	    for(int f=0; f<3; ++f)
	      for(int a=0; a<ndofs; ++a)
		{
		  auto& var = varPD[f][a];
		  var.SetZero();
		  var.PD = &PD;
		  var.eta0[f] = Elm->GetShape(fields[f],qx,a);
		  var.d1Eta0[f] = Elm->GetDShape(fields[f],qx,a,0);

		  // Compute the variations of strains
		  var.ComputeVarStrains(qPtCoordinatesY[z],w);

		  // Compute variations of stress/moment resultants
		  Mat->ComputeVarResultants(var.vEpsilon, var.vRho, var.vDelta, var.vn, var.vm, var.vq);
		}
	
	    for(int f=0; f<3; ++f)
	      for(int a=0; a<ndofs; ++a)
		{
		  auto& var = varPD[f+3][a];
		  var.SetZero(); 
		  var.PD = &PD; 
		  var.eta1[f] = Elm->GetShape(fields[f+3],qx,a);
		  var.d1Eta1[f] = Elm->GetDShape(fields[f+3],qx,a,0);
 
		  // Compute the variations of strains
		  var.ComputeVarStrains(qPtCoordinatesY[z],w);

		  // Compute variations of stress/moment resultants
		  Mat->ComputeVarResultants(var.vEpsilon, var.vRho, var.vDelta, var.vn, var.vm, var.vq);
		}
	
	    for(int f=0; f<3; ++f)
	      for(int a=0; a<ndofs; ++a)
		{
		  auto& var = varPD[f+6][a];
		  var.SetZero();
		  var.PD = &PD;
		  var.eta2[f] = Elm->GetShape(fields[f+6],qx,a);
		  var.d1Eta2[f] = Elm->GetDShape(fields[f+6],qx,a,0);
		  
		  // Compute the variations of strains
		  var.ComputeVarStrains(qPtCoordinatesY[z],w);

		  // Compute variations of stress/moment resultants
		  Mat->ComputeVarResultants(var.vEpsilon, var.vRho, var.vDelta, var.vn, var.vm, var.vq);
		}
	
	    for(int a=0; a<ndofs; ++a)
	      {
		auto& var = varPD[9][a];
		var.SetZero();
		var.PD = &PD;
		var.delDirSpatial[0] = Elm->GetShape(fields[0],qx,a)*Lambda1[a][0][0];
		var.delDirSpatial[1] = Elm->GetShape(fields[0],qx,a)*Lambda1[a][1][0];
		var.delDirSpatial[2] = Elm->GetShape(fields[0],qx,a)*Lambda1[a][2][0];
		var.d1DelDirSpatial[0] = Elm->GetDShape(fields[0],qx,a,0)*Lambda1[a][0][0];
		var.d1DelDirSpatial[1] = Elm->GetDShape(fields[0],qx,a,0)*Lambda1[a][1][0];
		var.d1DelDirSpatial[2] = Elm->GetDShape(fields[0],qx,a,0)*Lambda1[a][2][0];
		
		// Compute the variations of strains
		var.ComputeVarStrains(qPtCoordinatesY[z],w);

		// Compute variations of stress/moment resultants
		Mat->ComputeVarResultants(var.vEpsilon, var.vRho, var.vDelta, var.vn, var.vm, var.vq);
	      }
	
	    for(int a=0; a<ndofs; ++a)
	      {
		auto& var = varPD[10][a];
		var.SetZero();
		var.PD = &PD;
		var.delDirSpatial[0] = Elm->GetShape(fields[0],qx,a)*Lambda1[a][0][1];
		var.delDirSpatial[1] = Elm->GetShape(fields[0],qx,a)*Lambda1[a][1][1];
		var.delDirSpatial[2] = Elm->GetShape(fields[0],qx,a)*Lambda1[a][2][1];
		var.d1DelDirSpatial[0] = Elm->GetDShape(fields[0],qx,a,0)*Lambda1[a][0][1];
		var.d1DelDirSpatial[1] = Elm->GetDShape(fields[0],qx,a,0)*Lambda1[a][1][1];
		var.d1DelDirSpatial[2] = Elm->GetDShape(fields[0],qx,a,0)*Lambda1[a][2][1];
		
		// Compute the variations of strains
		var.ComputeVarStrains(qPtCoordinatesY[z],w);

		// Compute variations of stress/moment resultants
		Mat->ComputeVarResultants(var.vEpsilon, var.vRho, var.vDelta, var.vn, var.vm, var.vq);
	      }

	    // Compute the residue
	    if(resval!=nullptr)
	      for(int f=0; f<11; ++f)
		for(int a=0; a<ndofs; ++a)
		  {
		    auto& var = varPD[f][a];

		    // n.delta(epsilon) + m.delta(rho) + q.delta(delta)
		    double tempval1 = 0.;
		    for(int i=0; i<2; ++i)
		      {
			tempval1 += PD.q[i]*var.vDelta[i];
			for(int j=0; j<2; ++j)
			  tempval1 += PD.n[i][j]*var.vEpsilon[i][j] + PD.m[i][j]*var.vRho[i][j];
		      }

		    // Update the residue
		    (*resval)[fields[f]][a] += qPtWtsY[z]*QwtsX[qx]*tempval1;
		  }

	    // Compute the stiffness
	    if(dresval!=nullptr)
	      for(int f=0; f<11; ++f)
		for(int a=0; a<ndofs; ++a)
		  {
		    const auto& var = varPD[f][a];
		    for(int g=0; g<11; ++g)
		      for(int b=0; b<ndofs; ++b)
			{
			  const auto& VAR = varPD[g][b];

			  // 2nd variation
			  VARvar.PD = &PD;
			  VARvar.deltaPD = &var;
			  VARvar.DELTAPD = &VAR;

			  // Compute the second variations of strains
			  VARvar.ComputeVarVarStrains();

			  // VAR( n.delta(epsilon) + m.delta(rho) + q.delta(delta) )
			  double tempval2 = 0.;
			  for(int i=0; i<2; ++i)
			    {
			      tempval2 += VAR.vq[i]*var.vDelta[i] + PD.q[i]*VARvar.vvDelta[i];
			      for(int j=0; j<2; ++j)
				tempval2 += VAR.vn[i][j]*var.vEpsilon[i][j] + VAR.vm[i][j]*var.vRho[i][j] + PD.n[i][j]*VARvar.vvEpsilon[i][j] + PD.m[i][j]*VARvar.vvRho[i][j];
			    }

			  // Update the stiffness
			  (*dresval)[fields[f]][a][fields[g]][b] += qPtWtsY[z]*QwtsX[qx]*tempval2;
			}
		  }
	  }
      }
    return;
  }



  // Consistency test for this operation at the given configuration
  bool Ops::ConsistencyTest(const void* arg,
			    const double pertEPS, const double tolEPS) const
  {
    // Cast configuration pointer
    const TriRodsConfiguration* ptr = static_cast<const TriRodsConfiguration*>(arg);
    assert(ptr!=nullptr && "TriRods::Ops::ConsistencyTest- invalid configuration pointer");
    const TriRodsConfiguration& Config = *ptr;
    assert(Config.IsInitialized() && "TriRods::Ops::ConsistencyTest- configuration not initialized");

    // Compute the residue and the derivative
    const int ndofs = GetFieldDof(0);
    const int nfields = static_cast<int>(fields.size());
    assert(static_cast<int>(fields.size())==11);
    std::vector<std::vector<double>> resvals(nfields, std::vector<double>(ndofs));
    std::vector<std::vector<std::vector<std::vector<double>>>>dresvals(nfields);
    for(int f=0; f<nfields; ++f)
      {
	dresvals[f].resize(ndofs);
	for(int a=0; a<ndofs; ++a)
	  {
	    dresvals[f][a].resize(nfields);
	    for(int g=0; g<nfields; ++g)
	      dresvals[f][a][g].resize(ndofs);
	  }
      }
    GetDVal(&Config, &resvals, &dresvals);
    
    ////
    // Rotation matrix array
    std::vector<Mat3> Lambda1;
    for(int a=0; a<ndofs; ++a)
      Lambda1[a] = Config.GetRotations(ConfigStr.nodesX[a]);
    ////

    // Compute the residue and dresidue numerically
    auto nresvals = resvals;
    auto ndresvals = dresvals;
    SetZero(&nresvals, &ndresvals);

    // Local residual calculations at perturbed configs
    auto pres = resvals;
    auto mres = resvals;

    // Compute the residual/dresidual numerically
    // Fields for centerline
    for(int f=0; f<3; ++f) 
      for(int a=0; a<ndofs; ++a)
	{
	  // Create perturbed configurations here
	  TriRodsConfiguration pConfig(Config), mConfig(Config);

	  // Centerline update
	  Vec3 incPhi0;
	  for(int i=0; i<3; ++i)
	    incPhi0[i] = 0.;
	  incPhi0[f] = pertEPS;
	  pConfig.UpdateLine0(ConfigStr.nodesX[a], incPhi0);
	  incPhi0[f] = -pertEPS;
	  mConfig.UpdateLine0(ConfigStr.nodesX[a], incPhi0);
	  
	  // No rotation/vorticity updates
	  pConfig.SetInitialized();
	  mConfig.SetInitialized();

	  // Compute perturbed energies
	  auto Eplus = GetEnergy(&pConfig);
	  auto Eminus = GetEnergy(&mConfig);
	  nresvals[f][a] = (Eplus-Eminus)/(2.*pertEPS);

	  // Compute the perturbed residual
	  GetVal(&mConfig, &mres);
	  GetVal(&pConfig, &pres);

	  // Update the numerical stiffnesses
	  for(int g=0; g<nfields; ++g)
	    for(int b=0; b<ndofs; ++b)
	      ndresvals[g][b][f][a] = (pres[g][b]-mres[g][b])/(2.*pertEPS);
	}

    // Fields for edge_1
    for(int f=0; f<3; ++f) 
      for(int a=0; a<ndofs; ++a)
	{
	  // Create perturbed configurations here
	  TriRodsConfiguration pConfig(Config), mConfig(Config);

	  // Centerline update
	  Vec3 incPhi1;
	  for(int i=0; i<3; ++i)
	    incPhi1[i] = 0.;
	  incPhi1[f] = pertEPS;
	  pConfig.UpdateLine0(ConfigStr.nodesX[a], incPhi1);
	  incPhi1[f] = -pertEPS;
	  mConfig.UpdateLine0(ConfigStr.nodesX[a], incPhi1);
	  
	  // No rotation/vorticity updates
	  pConfig.SetInitialized();
	  mConfig.SetInitialized();

	  // Compute perturbed energies
	  auto Eplus = GetEnergy(&pConfig);
	  auto Eminus = GetEnergy(&mConfig);
	  nresvals[f+3][a] = (Eplus-Eminus)/(2.*pertEPS);

	  // Compute the perturbed residual
	  GetVal(&mConfig, &mres);
	  GetVal(&pConfig, &pres);

	  // Update the numerical stiffnesses
	  for(int g=0; g<nfields; ++g)
	    for(int b=0; b<ndofs; ++b)
	      ndresvals[g][b][f+3][a] = (pres[g][b]-mres[g][b])/(2.*pertEPS);
	}

    // Fields for edge_2
    for(int f=0; f<3; ++f) 
      for(int a=0; a<ndofs; ++a)
	{
	  // Create perturbed configurations here
	  TriRodsConfiguration pConfig(Config), mConfig(Config);

	  // Centerline update
	  Vec3 incPhi2;
	  for(int i=0; i<3; ++i)
	    incPhi2[i] = 0.;
	  incPhi2[f] = pertEPS;
	  pConfig.UpdateLine0(ConfigStr.nodesX[a], incPhi2);
	  incPhi2[f] = -pertEPS;
	  mConfig.UpdateLine0(ConfigStr.nodesX[a], incPhi2);
	  
	  // No rotation/vorticity updates
	  pConfig.SetInitialized();
	  mConfig.SetInitialized();

	  // Compute perturbed energies
	  auto Eplus = GetEnergy(&pConfig);
	  auto Eminus = GetEnergy(&mConfig);
	  nresvals[f+6][a] = (Eplus-Eminus)/(2.*pertEPS);

	  // Compute the perturbed residual
	  GetVal(&mConfig, &mres);
	  GetVal(&pConfig, &pres);

	  // Update the numerical stiffnesses
	  for(int g=0; g<nfields; ++g)
	    for(int b=0; b<ndofs; ++b)
	      ndresvals[g][b][f+6][a] = (pres[g][b]-mres[g][b])/(2.*pertEPS);
	}

    // Fields for director
    for(int f=0; f<2; ++f)
      for(int a=0; a<ndofs; ++a)
	{
	  // Create perturbed configurations
	  TriRodsConfiguration pConfig(Config), mConfig(Config);

	  // incDirMat at a node
	  Vec2 incT;
	  for(int i=0; i<2; ++i)
	    incT[i] = 0.;
	  incT[f] = pertEPS;

	  // incDirSpat at a node
	  Vec3 incT3, inct;
	  for(int i=0; i<2; ++i)
	    incT3[i] = incT[i];
	  incT3[2] = 0.;

	  for(int i=0; i<3; ++i)
	    {
	      inct[i] = 0.;
	      for(int j=0; j<3; ++j)
		inct[i] += Lambda1[a][i][j]*incT3[j];
	    }

	  // Interpolate incDirSpat at the quadrature points and update t and tPrime at the quadrature points
	  for(int qx=0; qx<ConfigStr.nQuadsPerElmX; ++qx)
	    {
	      Vec3 incDirSpat, incD1DirSpat;
	      for(int i=0; i<3; ++i)
		{
		  incDirSpat[i] = inct[i]*Elm->GetShape(fields[0], qx, a);
		  incD1DirSpat[i] = inct[i]*Elm->GetDShape(fields[0], qx, a,0);
		}
	      pConfig.UpdateDirSpatial(ConfigStr.elmX, qx, incDirSpat);
	      pConfig.UpdateD1DirSpatial(ConfigStr.elmX, qx, incDirSpat, incD1DirSpat);
	      for(int i=0; i<3; ++i)
		{
		  incDirSpat[i] *= -1.;
		  incD1DirSpat[i] *= -1.;
		}
	      mConfig.UpdateDirSpatial(ConfigStr.elmX, qx, incDirSpat);
	      mConfig.UpdateD1DirSpatial(ConfigStr.elmX, qx, incDirSpat, incD1DirSpat);
	    }

	  // Update rotations at the nodal points  
	  pConfig.UpdateRotations(ConfigStr.nodesX[a], incT);
	  incT[f] *= -1;
	  mConfig.UpdateRotations(ConfigStr.nodesX[a], incT);

	  // No other field updates
	  pConfig.SetInitialized();
	  mConfig.SetInitialized();

	  // Compute perturbed energies
	  auto Eplus = GetEnergy(&pConfig);
	  auto Eminus = GetEnergy(&mConfig);
	  nresvals[f+9][a] = (Eplus-Eminus)/(2.*pertEPS);

	  // Compute the residuals at the perturbed configuration
	  GetVal(&mConfig, &mres);
	  GetVal(&pConfig, &pres);

	  // numerical stiffness values
	   for(int g=0; g<nfields; ++g)
	     for(int b=0; b<ndofs; ++b)
	       ndresvals[g][b][f+9][a] = (pres[g][b]-mres[g][b])/(2.*pertEPS);
	}

    // Check consistency
    for(int f=0; f<11; ++f)
      for(int a=0; a<ndofs; ++a)
	{
	  assert(std::abs(resvals[f][a]-nresvals[f][a])<tolEPS && "TriRods::Ops::Consistency test for residuals failed");
	  for(int g=0; g<11; ++g)
	    for(int b=0; b<ndofs; ++b)
	      assert(std::abs(dresvals[f][a][g][b]-ndresvals[f][a][g][b])<tolEPS && "TriRods::Ops::Consistency test for dresiduals failed");
	}

    // For debugging puposes only
    if(0)
      {
	std::cout<<"\n\nConsistency of residuals: ";
	for(int f=0; f<11; ++f)
	  for(int a=0; a<ndofs; ++a)
	    std::cout<<"\n"<<resvals[f][a]<<" should be "<<nresvals[f][a];
	std::fflush( stdout );

	std::cout<<"\n\nConsistency of dresiduals: ";
	for(int f=0; f<11; ++f)
	  for(int a=0; a<ndofs; ++a)
	    for(int g=0; g<11; ++g)
	      for(int b=0; b<ndofs; ++b)
		std::cout<<"\n"<<dresvals[f][a][g][b]<<" should be "<<ndresvals[f][a][g][b];
	std::fflush( stdout );
      }
    return true;  
  }
}
