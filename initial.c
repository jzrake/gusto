#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "gusto.h"
#include "quartic.h"



static const char **id_uniform(struct gusto_user *user,
			       struct aux_variables *A, double *X)
{
  if (A == NULL) {
    static const char *help[] = { "-- Uniform initial data --",
				  "density1: uniform density",
				  "pressure1: uniform pressure",
				  NULL };
    return help;
  }

  A->comoving_mass_density = user->density1;
  A->gas_pressure = user->pressure1;

  return NULL;
}



const char **id_cylindrical_shock(struct gusto_user *user,
				  struct aux_variables *A, double *X)
{
  if (A == NULL) {
    static const char *help[] = { "-- Cylindrical shock wave --",
				  NULL };
    return help;
  }

  double x = X[1] - 0.5;
  double y = X[3] - 0.5;
  if (sqrt(x*x + y*y) < 0.125) {
    A->comoving_mass_density = 1.0;
    A->gas_pressure = 1.0;
  }
  else {
    A->comoving_mass_density = 0.1;
    A->gas_pressure = 0.125;
  }

  return NULL;
}



const char **id_sound_wave(struct gusto_user *user,
			   struct aux_variables *A, double *X)
{
  if (A == NULL) {
    static const char *help[] = { "-- Sound wave --",
				  NULL };
    return help;
  }
  double k = 2 * M_PI;
  double z = X[3];

  double gm = gamma_law_index;
  double d0 = 1.0;
  double p0 = 1.0;
  double u0 = 0.0;
  double h0 = 1.0 + (p0/d0) * gm / (gm - 1);
  double cs = sqrt(gm * p0 / d0 / h0);
  double dd = 0.1;
  double dp = cs * cs * dd;
  double du = cs * dd / d0;

  double d = d0 + dd * sin(k * z);
  double p = p0 + dp * sin(k * z);
  double u = u0 + du * sin(k * z);

  A->velocity_four_vector[3] = u;
  A->comoving_mass_density = d;
  A->gas_pressure = p;

  return NULL;
}



const char **id_density_wave(struct gusto_user *user,
			     struct aux_variables *A, double *X)
{
  if (A == NULL) {
    static const char *help[] = { "-- Density wave --",
				  "pressure0: ambient gas pressure",
				  "density0: ambient density",
				  "density1: density perturbation",
				  "fourvel0: four velocity in the R direction",
				  "fourvel1: four velocity in the z direction",
				  NULL };
    return help;
  }
  double k = 2 * M_PI;
  double R = X[1];
  double z = X[3];
  double d = user->density0 + user->density1 * sin(k * z + k * R);

  A->velocity_four_vector[1] = user->fourvel0;
  A->velocity_four_vector[3] = user->fourvel1;
  A->comoving_mass_density = d;
  A->gas_pressure = user->pressure0;

  return NULL;
}



const char **id_abc_ff(struct gusto_user *user,
		       struct aux_variables *A, double *X)
{
  if (A == NULL) {
    static const char *help[] =
      { "-- ABC magnetic field --",
	"fourvel0: four velocity in the R direction",
	"fourvel1: four velocity in the z direction",
	NULL };
    return help;
  }

  double a = user->abc[0];
  double b = user->abc[1];
  double c = user->abc[2];

  double al = 2 * M_PI;
  double B1 = c * cos(X[3] * al) - b * sin(X[2] * al);
  double B2 = a * cos(X[1] * al) - c * sin(X[3] * al);
  double B3 = b * cos(X[2] * al) - a * sin(X[1] * al);

  A->magnetic_four_vector[1] = B1;
  A->magnetic_four_vector[2] = B2;
  A->magnetic_four_vector[3] = B3;

  A->vector_potential = B2 / al;
  A->velocity_four_vector[1] = user->fourvel0;
  A->velocity_four_vector[3] = user->fourvel1;

  return NULL;
}



const char **id_cyl_ff(struct gusto_user *user,
		       struct aux_variables *A, double *X)
{
  if (A == NULL) {
    static const char *help[] =
      { "-- Cylindrical force-free magnetic field --",
	NULL };
    return help;
  }

  double R = X[1];
  double B1 = 0.0;
  double B2 = j1(R);
  double B3 = j0(R);

  A->magnetic_four_vector[1] = B1;
  A->magnetic_four_vector[2] = B2;
  A->magnetic_four_vector[3] = B3;

  return NULL;
}



const char **id_michel69(struct gusto_user *user,
			 struct aux_variables *A, double *X)
{
  if (A == NULL) {
    static const char *help[] =
      {"-- Michel (1969) cold MHD wind --",
       "Use with spherical coordinate source terms and radial mesh if in 1D.",
       "",
       "entropy: uniform entropy (choose something small)",
       "sigma: magnetization, ratio of r0 to light cylinder radius (squared)",
       NULL};
    return help;
  }

  double r = X[1];
  double c = 1.0;
  double Phi = 1.0;
  double sig = user->sigma;
  double eta = pow(sig, 1./3);
  double lam = (pow(1 + eta * eta, 1.5) - 1) / sig;
  double rL = 1.0; /* light cylinder radius */
  double r0 = rL * sqrt(sig);
  double Br = Phi / (r * r);
  double xc2 = lam / (1 + sig * lam); /* critical point, xc^2 */

  double x2 = pow(r / r0, 2);
  double x4 = pow(r / r0, 4);
  double d0 = -sig * x2 + sig * sig * x4;
  double d1 = 2 * sig * x4;
  double d2 = 1 + sig * x2 * (lam * lam - 2) + sig * x4 * (sig - lam * (2 + lam * sig));
  double d3 = -2 * x2 + 2 * sig * x4;
  double d4 = x4;
  double roots[4];

  if (solve_quartic_equation(d4, d3, d2, d1, d0, roots) != 4) {
    printf("[gusto] WARNING: michel69 solution has imaginary roots\n");
  }

  double ur = x2 < xc2 ? roots[1] : roots[2];
  double uf = sqrt(sig * x2) * (1 - lam * ur) / (1 - x2 * (ur + sig));
  double Bf = sqrt(sig * x2) * Br * (x2 * (1 + lam * sig) - lam) / (1 - x2 * (ur + sig));
  double u0 = sqrt(1.0 + ur*ur + uf*uf);
  double b0 = Br*ur + Bf*uf;
  double br = (Br + b0 * ur) / u0;
  double bf = (Bf + b0 * uf) / u0;
  double f = Phi * Phi / (4 * M_PI * c * r0 * r0); /* mass rate per steradian */
  double d = f / (r * r * ur);
  double s = user->entropy; /* log(p / rho^Gamma) */
  double p = exp(s) * pow(d, gamma_law_index);

  A->velocity_four_vector[1] = ur;
  A->velocity_four_vector[2] = uf;
  A->magnetic_four_vector[1] = br / sqrt(4 * M_PI); /* code units of B field */
  A->magnetic_four_vector[2] = bf / sqrt(4 * M_PI);
  A->comoving_mass_density = d;
  A->gas_pressure = p;

  return NULL;
}



const char **id_michel73(struct gusto_user *user,
			 struct aux_variables *A, double *X)
{
  if (A == NULL) {
    static const char *help[] =
      {"-- Michel (1973) split-monopole pulsar magnetosphere --",
       NULL};
    return help;
  }

  double R = X[1];
  double z = X[3];
  double r = sqrt(R*R + z*z);
  double omega = 0.5 * (X[0] < 1 ? X[0] : 1); /* ramps up */

  double Bz = z / pow(r, 3);
  double BR = R / pow(r, 3);
  double Bf = -omega * R / pow(r, 2);

  double Y = 1 - z / r;

  double gradY[4] = {0, R*Bz, 0, -R*BR};
  double B[4] = {0, BR, Bf, Bz};
  double E[4] = {0, -omega * gradY[1], 0, -omega * gradY[3]};
  double S[4] = VEC4_CROSS(E, B);
  double U = VEC4_DOT(B,B);
  double vm[4] = { 0, S[1]/U, S[2]/U, S[3]/U };
  double u0 = pow(1.0 - VEC4_DOT(vm, vm), -0.5);
  double uR = vm[1] * u0;
  double uf = vm[2] * u0;
  double uz = vm[3] * u0;
  double b0 = uR*BR + uf*Bf + uz*Bz;
  double Bp = sqrt(BR*BR + Bz*Bz);
  double up = sqrt(uR*uR + uz*uz);
  double dg = user->density0;// / user->sigma * Bp / up; if (dg > 1) dg = 1;
  double s0 = user->entropy; /* log(p / rho^Gamma) */
  double pg = exp(s0) * pow(dg, gamma_law_index);

  /* printf("B^2/d=%f\n", Bp*Bp/dg); */
  /* double dg = user->density0; */
  /* double pg = user->pressure0; */

  A->velocity_four_vector[1] = uR;
  A->velocity_four_vector[2] = uf;
  A->velocity_four_vector[3] = uz;
  A->magnetic_four_vector[1] = (BR + b0 * uR) / u0;
  A->magnetic_four_vector[2] = (Bf + b0 * uf) / u0;
  A->magnetic_four_vector[3] = (Bz + b0 * uz) / u0;
  A->comoving_mass_density = dg;
  A->gas_pressure = pg;
  A->vector_potential = Y / R;


  /*
   * [x] constant poloidal flux
   * [x] constant mass flux
   * [x] constant energy flux --- need to use v_min
   * [x] constant L flux
   */

  /* double nhat[4] = {0, B[1]/Bp, 0, B[3]/Bp}; */
  /* double SFlux = VEC4_DOT(S, nhat); */
  /* double dA = 1.0 / Bp; */
  /* double F[8]; */
  /* gusto_complete_aux(A); */
  /* gusto_fluxes(A, nhat, F); */
  /* printf("S.dA=%f F[TAU].dA=%f\n", SFlux * dA, (F[TAU] + F[DDD]) * dA); */

  return NULL;
}



const char **id_narayan07(struct gusto_user *user,
			  struct aux_variables *A, double *X)
{
  if (A == NULL) {
    static const char *help[] =
      {"-- Narayan et al. (2007) self-similar force-free disk winds --",
       NULL};
    return help;
  }

  double M = 0.2;
  double s = 2.0;
  double R = X[1];
  double z = X[3];
  double r = sqrt(R*R + z*z);

  double Bz = 1.0 / r;
  double BR = (1.0 - z / r) / R;
  double Bf = -s * M / R; /* Equation (45) N07 */
  double Bp = sqrt(BR*BR + Bz*Bz);

  double Y = r - z;       /* the flux function */
  double T = Y / R;
  double omega = M / (T * R);

  double gradY[4] = {0, R*Bz, 0, -R*BR};
  double B[4] = {0, BR, Bf, Bz};
  double E[4] = {0, -omega * gradY[1], 0, -omega * gradY[3]};
  double S[4] = VEC4_CROSS(E, B);
  double U = VEC4_DOT(B,B);
  double vm[4] = { 0, S[1]/U, S[2]/U, S[3]/U };
  double u0 = pow(1.0 - VEC4_DOT(vm, vm), -0.5);
  double uR = vm[1] * u0;
  double uf = vm[2] * u0;
  double uz = vm[3] * u0;
  double up = sqrt(uR*uR + uz*uz);
  double b0 = uR*BR + uf*Bf + uz*Bz;
  double dg = user->density0 * Bp / up;
  double s0 = user->entropy; /* log(p / rho^Gamma) */
  double pg = exp(s0) * pow(dg, gamma_law_index);

  A->velocity_four_vector[1] = uR;
  A->velocity_four_vector[2] = uf;
  A->velocity_four_vector[3] = uz;
  A->magnetic_four_vector[1] = (BR + b0 * uR) / u0;
  A->magnetic_four_vector[2] = (Bf + b0 * uf) / u0;
  A->magnetic_four_vector[3] = (Bz + b0 * uz) / u0;
  A->comoving_mass_density = dg;
  A->gas_pressure = pg;
  A->vector_potential = Y / R;

  /*
   * [x] constant poloidal flux
   * [x] constant mass flux
   * [x] constant energy flux --- need to use v_min
   * [x] constant L flux
   */

  /* double nhat[4] = {0, B[1]/Bp, 0, B[3]/Bp}; */
  /* double SFlux = VEC4_DOT(S, nhat); */
  /* double dA = 1.0 / Bp; */
  /* double F[8]; */
  /* gusto_complete_aux(A); */
  /* gusto_fluxes(A, nhat, F); */
  /* printf("S.dA=%f F[TAU].dA=%f\n", SFlux * dA, (F[TAU] + F[DDD]) * dA); */

  return NULL;
}



OpInitialData gusto_lookup_initial_data(const char *user_key)
{
  const char *keys[] = {
    "uniform",
    "cylindrical_shock",
    "density_wave",
    "abc_ff",
    "cyl_ff",
    "michel69",
    "michel73",
    "narayan07",
    NULL
  } ;
  OpInitialData vals[] = {
    id_uniform,
    id_cylindrical_shock,
    id_density_wave,
    id_abc_ff,
    id_cyl_ff,
    id_michel69,
    id_michel73,
    id_narayan07,
    NULL } ;
  int n = 0;
  const char *key;
  OpInitialData val;
  do {
    key = keys[n];
    val = vals[n];
    if (key == NULL) {
      break;
    }
    else if (strcmp(user_key, key) == 0) {
      return val;
    }
    else {
      n += 1;
    }
  } while (1);
  printf("[gusto] ERROR: no such initial_data=%s\n", user_key);
  return NULL;
}
