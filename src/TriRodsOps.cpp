#include <TriRodsOps.h>
#include <TriRodsPointKinematics.h>
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
    const Configuration* Config = static_cast<const Configuration*>(arg);
    assert(Config!=nullptr && "TriRods::Ops::GetEnergy() - Invalid pointer to configuration");

    // Check that this configuration is "set"
    assert(Config->IsInitialized() && "Ops::GetEnergy()- Configuration not initialized");
    
    // Energy to be computed
    double Energy = 0.;

    // Access workspace for quadrature point-level calculatons
    auto& PD = WrkSpc->PD;
    
    // Number of quadrature points
    const auto& nQuads = ConfigStr.nQuadsPerElmX;
    assert(nQuads==static_cast<int>(Elm->GetIntegrationWeights(fields[0]).size()));

    // Node number for this element
    const auto& nodes = ConfigStr.nodesX;
    assert(static_cast<int>(nodes.size())==GetFieldDof(fields[0]));
    const int ndofs = static_cast<int>(nodes.size());
    
    // This element number
    const auto& e = ConfigStr.elm;

    // Quadrature weights
    const auto& Qwts = Elm->GetIntegrationWeights(fields[0]);

    // Integrate
    for(int )
  }
}
