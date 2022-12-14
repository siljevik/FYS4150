#include <cmath>

// fetch class definition
#include "../header/particle.hpp"
#include "../header/penningtrap.hpp"

// Constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
{
  // Assigning member variables
  B0 = B0_in;
  V0 = V0_in;
  d = d_in;
}

// Method to add particle to trap
void PenningTrap::add_particle(Particle p_in)
{
  particles.push_back(p_in);
}

// Method to calculate ext. E-field
arma::mat PenningTrap::external_E_field(arma::mat R)
{
  int N = particles.size();
  arma::mat E = arma::mat(3, N);

  const double A = V0 / (d * d);

  for (int i = 0; i < N; i++)
  {
    arma::vec r = R.col(i);
    E.col(i) = A * arma::vec{r(0), r(1), -2. * r(2)};
  }
  return E;
}

// Method to calculate ext. B-field (magnetic field)
arma::mat PenningTrap::external_B_field(arma::mat V)
{
  int N = particles.size(); // Number of particles
  arma::vec B = arma::vec{0.0, 0.0, B0};
  arma::mat B_tot = arma::mat(3, N);

  for (int i = 0; i < N; i++)
  {
    arma::vec v_i = V.col(i);
    double q_i = particles[i].q;

    B_tot.col(i) = arma::cross(q_i * v_i, B);
  }
  return B_tot;
}

// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(arma::vec r_i, arma::vec r_j,
                                      double q_i, double q_j)
{
  // Electric constant in u(mu*m)^3/(mu*s)^2e^2
  const double k_e = 1.38935333 * pow(10,5);

  arma::vec E = k_e * q_j * (r_i - r_j) / std::pow(arma::norm(r_i - r_j), 3.0);

  return q_i * E;
}

// // Total force on particle_i from other particles
arma::mat PenningTrap::total_force_particles(arma::mat R)
{
  int N = particles.size();

  arma::mat F_tot = arma::mat(3, N);

  for (int i = 0; i < N; i++)
  {
    arma::vec r_i = R.col(i);
    double q_i = particles[i].q;

    for (int j = 0; j < N; j++)
    {
      if (j == i)
      {}

        else
        {
          double q_j = particles[j].q;
          arma::vec r_j = R.col(j);

          F_tot.col(i) += PenningTrap::force_particle(r_i, r_j, q_i, q_j);
        }
      }
    }
    return F_tot;
  }

// Total ext. force
arma::mat PenningTrap::total_force_external(arma::mat R, arma::mat V)
{
  arma::mat extE = PenningTrap::external_E_field(R);
  arma::mat extB = PenningTrap::external_B_field(V);

  arma::mat extF_tot = extE + extB;

  return extF_tot;
}

// Total force on particle_i from external field and other particles
arma::mat PenningTrap::total_force(arma::mat R, arma::mat V, bool particle_interaction)
{
  arma::mat F_tot;
// For particle interaction we consider the total force form external fields and particles
  if (particle_interaction)
  {
    F_tot = PenningTrap::total_force_external(R, V)
    + PenningTrap::total_force_particles(R);
  }
// With no interaction we only consider the external forces
  else
  {
    F_tot = PenningTrap::total_force_external(R, V);
  }

  return F_tot;
}

// Evolve the system one time step using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt, bool particle_interaction)
{

  int N = particles.size(); // Number of particles

  arma::mat R = arma::zeros(3, N); // Position
  arma::mat V = arma::zeros(3, N); // Velocity
  arma::mat a = arma::zeros(3, N); // Acceleration
  double m = particles[0].m;       // Mass of the particle 1

  arma::mat K1_r, K2_r, K3_r, K4_r; // Position
  arma::mat K1_v, K2_v, K3_v, K4_v; // Velocity

  // Looping through each of the N particles to update position and velocity
  for (int i = 0; i < N; i++)
  {
    R.col(i) = particles[i].r;
    V.col(i) = particles[i].v;
  }

  // Acceleration of particle (Newtons 2nd law):
  a = (1/m) * total_force(R, V, particle_interaction);

  //K1
  K1_v = dt * a;
  K1_r = dt * V;

  a = (1/m) * total_force(R + 0.5 * K1_r, V + 0.5 * K1_v, particle_interaction);

  //K2
  K2_v = dt * a;
  K2_r = dt * (V + 0.5 * K1_v);

  a = (1/m) * total_force(R + 0.5 * K2_r, V + 0.5 * K2_v, particle_interaction);

  //K3
  K3_v = dt * a;
  K3_r = dt * (V + .5 * K2_v);

  a = (1/m) * total_force(R + K3_r, V + K3_v, particle_interaction);

  //K4
  K4_v = dt * a;
  K4_r = dt * (V + K3_v);

  // Adding the values to position and vector matrices
  R += (1.0/6) * (K1_r + 2 * K2_r + 2 * K3_r + K4_r);
  V += (1.0/6) * (K1_v + 2 * K2_v + 2 * K3_v + K4_v);

  for (int i = 0; i < N; i++)
  {
    // Updating the old position/velocity to the new one calculated above
    particles[i].r = R.col(i);
    particles[i].v = V.col(i);
  }
}

// Evolve the system one time step using Forward Euler
void PenningTrap::forward_euler(double dt, bool particle_interaction)
{
  int N = particles.size();
  arma::mat R = arma::zeros(3, N);
  arma::mat V = arma::zeros(3, N);
  arma::mat a = arma::zeros(3, N);
  double m = particles[0].m;

  for (int i = 0; i < N; i++)
  {
    R.col(i) = particles[i].r;
    V.col(i) = particles[i].v;
  }

  a = (1/m) * total_force(R, V, particle_interaction);

  R += dt * V;
  V += dt * a;

  for (int i = 0; i < N; i++)
  {
    particles[i].r = R.col(i);
    particles[i].v = V.col(i);
  }
}
