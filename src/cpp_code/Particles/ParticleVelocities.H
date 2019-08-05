#ifndef _PARTICLEVELOCITIES_H_
#define _PARTICLEVELOCITIES_H_
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "DBox.H"
#include "RectMDArray.H"
#include "Hockney.H"
#include "ParticleSet.H"
using namespace std;
/// Class for computing the RHS in RK4, corresponding to the template parameter F.
class ParticleVelocities
{
public:
  /// Member function that computes increment to solution.Conforms to the m_f(...) interface in RK4. 
  /** Implements an operator that evaluates a_dt*F(t,a_X+a_k), given the inputs (a_k, a_time, a_dt, a_X ). It first computes a temporary a_X + a_k, then evaluates the right-hand side, scales the result by  a_dt and stores it in  a_k. 
   */

  ParticleVelocities();
  void operator()(ParticleShift& a_k, 
                     const double& a_time, const double& dt, 
                     ParticleSet& a_state);
};
#endif
