#ifndef TRI_RODS_UTILS_H
#define TRI_RODS_UTILS_H

#include <array>

using Vec3 = std::array<double,3>;

//Utility functions
namespace TriRods
{
  // Computes the axial vector of a 3x3 skew symmetric tensor
  void AxialVector(const double skw[][3], double* theta);

  // Computes the skew symmetric matrix for a given axial vector
  void HodgeStar(const double* theta, double skw[][3]);

  // Computes the exponential map given an axial vector
  void ExpSO3(const double* vec, double Mat[][3]);

  // Computes the exponential map given a skew symmetric tensor
  void ExpSO3(const double skw[][3], double Mat[][3]);

  // Comutes cross-product of two vectors
  Vec3 CrossProduct(const Vec3 v1, const Vec3 v2);
}

#endif
