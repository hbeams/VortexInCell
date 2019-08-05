#include "ParticleSet.H"

void ParticleShift::init(const ParticleSet& a_particles)
{
  int sizeOfa_particle = a_particles.m_particles.size();
  m_particles.resize(sizeOfa_particle);
  for (int j = 0; j<m_particles.size(); j++)
    {
      m_particles[j] *= 0.0;
    }
}
  /// m_particles[k] += a_rhs.m_particles[k]*a_scale.
void ParticleShift::increment(
                 double a_scale, 
                 const ParticleShift& a_rhs)
{
  for (int j = 0; j<a_rhs.m_particles.size(); j++)
    {
      m_particles[j].increment(a_scale, a_rhs.m_particles[j]);
    }
}
  /// m_particles[k] *= a_scale
void ParticleShift::operator*=(double a_scale)
{
  for (int j = 0; j<m_particles.size(); j++)
    {
      m_particles[j] *= a_scale;
    }
}
  /// reinitializes the values m_particles[k] to zero. Not used in RK4.
void ParticleShift::setToZero()
{
  for (int j = 0; j<m_particles.size(); j++)
    {
      m_particles[j] *= 0.0;
    }
}
