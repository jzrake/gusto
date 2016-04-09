#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "utlist.h"
#include "gusto.h"
#include "srmhd_c2p.h"
#include "quartic.h"



const char **initial_data_uniform(struct gusto_user *user,
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



const char **initial_data_cylindrical_shock(struct gusto_user *user,
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



const char **initial_data_sound_wave(struct gusto_user *user,
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



const char **initial_data_density_wave(struct gusto_user *user,
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



const char **initial_data_abc_ff(struct gusto_user *user,
				 struct aux_variables *A, double *X)
{
  if (A == NULL) {
    static const char *help[] = { "-- ABC magnetic field --",
				  NULL };
    return help;
  }

  double a = 1.0;
  double b = 0.0;
  double c = 1.0;

  double al = 2 * M_PI;
  double B1 = c * cos(X[3] * al) - b * sin(X[2] * al);
  double B2 = a * cos(X[1] * al) - c * sin(X[3] * al);
  double B3 = b * cos(X[2] * al) - a * sin(X[1] * al);

  A->magnetic_four_vector[1] = B1;
  A->magnetic_four_vector[2] = B2;
  A->magnetic_four_vector[3] = B3;

  A->vector_potential = B2 / al;

  return NULL;
}



const char **initial_data_cyl_ff(struct gusto_user *user,
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



const char **initial_data_michel69(struct gusto_user *user,
				   struct aux_variables *A, double *X)
{
  if (A == NULL) {
    static const char *help[] =
      {"-- Michel 1969 cold MHD wind --",
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



OpInitialData gusto_lookup_initial_data(const char *user_key)
{
  const char *keys[] = {
    "uniform",
    "cylindrical_shock",
    "density_wave",
    "abc_ff",
    "cyl_ff",
    "michel69",
    NULL
  } ;
  OpInitialData vals[] = {
    initial_data_uniform,
    initial_data_cylindrical_shock,
    initial_data_density_wave,
    initial_data_abc_ff,
    initial_data_cyl_ff,
    initial_data_michel69,
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



void gusto_initial_data(struct gusto_sim *sim)
{
  struct mesh_vert *V;
  struct mesh_cell *C;
  struct aux_variables *A;

  for (int n=0; n<sim->num_rows; ++n) {

    /* Vertex initial data is only used to get the vector potential. We don't
       need a completed set of aux variables there. */
    DL_FOREACH(sim->rows[n].verts, V) {
      A = &V->aux[0];
      gusto_default_aux(A);
      sim->initial_data(&sim->user, A, V->x);
    }

    /* Cell initial data is used to get everything else. The aux variables are
       then completed so we can use them after differencing the vector potential
       to get the magnetic field at the cell centers. */
    DL_FOREACH(sim->rows[n].cells, C) {
      A = &C->aux[0];
      gusto_default_aux(A);
      sim->initial_data(&sim->user, A, C->x);
      gusto_complete_aux(A);
      gusto_to_conserved(A, C->U, C->dA);
    }
  }


  if (sim->user.advance_poloidal_field) {

    gusto_compute_cell_magnetic_field(sim);

    for (int n=0; n<sim->num_rows; ++n) {
      DL_FOREACH(sim->rows[n].cells, C) {
	A = &C->aux[0];
	double b0 = A->magnetic_four_vector[0];
	double u0 = A->velocity_four_vector[0];
	double u1 = A->velocity_four_vector[1];
	double u3 = A->velocity_four_vector[3];
	double B1 = C->U[B11] / C->dA[1];
	double B3 = C->U[B33] / C->dA[3];
	double b1 = (B1 + b0 * u1) / u0;
	double b3 = (B3 + b0 * u3) / u0;

	A->magnetic_four_vector[1] = b1;
	A->magnetic_four_vector[3] = b3;
	gusto_complete_aux(A);
	gusto_to_conserved(A, C->U, C->dA);
      }
    }
  }
}



void gusto_compute_fluxes(struct gusto_sim *sim)
{
  struct mesh_face *F;
  struct mesh_vert *V;

  for (int n=0; n<sim->num_rows; ++n) {
    DL_FOREACH(sim->rows[n].verts, V) {
      V->Efield = 0.0;
      V->num_counts = 0;
    }
  }

  DL_FOREACH(sim->faces, F) {

    struct mesh_cell *CL = F->cells[0];
    struct mesh_cell *CR = F->cells[1];

    double *v0 = F->verts[0]->v;
    double *v1 = F->verts[1]->v;
    double vpar = 0.5 * (VEC4_DOT(v0, F->nhat) + VEC4_DOT(v1, F->nhat));

    if (CL && CR) {
      gusto_riemann(CL->aux, CR->aux, F->nhat, F->Fhat, vpar);
    }
    else if (CL) {
      gusto_riemann(CL->aux, CL->aux, F->nhat, F->Fhat, vpar);
    }
    else if (CR) {
      gusto_riemann(CR->aux, CR->aux, F->nhat, F->Fhat, vpar);
    }

    /*
     * Compute electric field at the face endpoints (x0, x1) from the Godunov
     * flux.
     *
     * E2 = -F3(B1) = F1(B3) = F.n[B.dx] = Fhat[B1] dx1 + Fhat[B3] dx3
     *
     * where dx is the unit vector from x0 -> x1.
     */
    double df[4] = {0, 0, 1, 0};
    double dx[4] = VEC4_CROSS(F->nhat, df); /* unit vector from x0 -> x1 */
    double ef = dx[1] * F->Fhat[B11] + dx[3] * F->Fhat[B33];

    F->verts[0]->Efield += ef;
    F->verts[1]->Efield += ef;
    F->verts[0]->num_counts += 1;
    F->verts[1]->num_counts += 1;
  }

  for (int n=0; n<sim->num_rows; ++n) {
    DL_FOREACH(sim->rows[n].verts, V) {
      V->Efield /= V->num_counts;
      V->num_counts = 0;
    }
  }
}



void gusto_compute_cell_magnetic_field(struct gusto_sim *sim)
/*
 * Evaluate the poloidal (in-plane) magnetic field at the cell centers from the
 * toroidal vector potential A.phi stored at cell vertices. The operation leaves
 * the result in U[B11] and U[B33], which are magnetic fluxes, being multiplied
 * by the nominal cell cross-sections dA[1] and dA[3] to be consistent with the
 * convention used for the toroidal (out-of-plane) field.
 */
{
  struct mesh_cell *C;

  for (int n=0; n<sim->num_rows; ++n) {
    DL_FOREACH(sim->rows[n].cells, C) {

      /*
       *
       *  (+)1-------3(+)
       *     |       |
       *     |   x   |
       *     |       |
       *  (+)0-------2(+)
       *
       *  x : [curl(1,2;3,0)]
       *  + : [curl(0,1;0,2) + curl(1,3;1,0) + curl(2,0;2,3) + curl(3,2;3,1)]/4
       */

      if (sim->user.curl_mode == '+') {
	/* Evaluate curl A using four estimates of curlA, one at each vertex
	   based on the difference of A biased toward the cell center. */
	double crlA0[2];
	double crlA1[2];
	double crlA2[2];
	double crlA3[2];
	gusto_curlA1(C, 0, 1, 2, crlA0);
	gusto_curlA1(C, 1, 3, 0, crlA1);
	gusto_curlA1(C, 2, 0, 3, crlA2);
	gusto_curlA1(C, 3, 2, 1, crlA3);
	C->U[B11] = (crlA0[0] + crlA1[0] + crlA2[0] + crlA3[0]) / 4 * C->dA[1];
	C->U[B33] = (crlA0[1] + crlA1[1] + crlA2[1] + crlA3[1]) / 4 * C->dA[3];
      }
      else if (sim->user.curl_mode == 'x') {
	/* Evaluate curl A using all four vertices at once to get one estimate
	   at (or near) the cell center. */
	double crlA[2];
	gusto_curlA2(C, crlA);
	C->U[B11] = crlA[0] * C->dA[1];
	C->U[B33] = crlA[1] * C->dA[3];
      }
    }
  }
}



void gusto_transmit_fluxes(struct gusto_sim *sim, double dt)
{
  struct mesh_face *F;

  DL_FOREACH(sim->faces, F) {
    struct mesh_cell *CL = F->cells[0];
    struct mesh_cell *CR = F->cells[1];

    for (int q=0; q<5; ++q) {
      if (CL) CL->U[q] -= F->Fhat[q] * F->nhat[0] * dt;
      if (CR) CR->U[q] += F->Fhat[q] * F->nhat[0] * dt;
    }

    if (CL) CL->U[B22] -= F->Fhat[B22] * F->length * dt;
    if (CR) CR->U[B22] += F->Fhat[B22] * F->length * dt;
  }
}



void gusto_advance_vector_potential(struct gusto_sim *sim, double dt)
/*
 * Use the Efield data stored on cell vertices to advance the vector potential
 * data at the same location. Efield is part of the vertex data itself, and was
 * evaluated at the compte_fluxes stage.
 *
 *                           \dot{A} = -E
 *
 * E will be the motional electric field E = -(v - v0) \times B where v0 is the
 * vertex velocity. So if the vertex moves precisely with the fluid velocity, E
 * will be zero.
 */
{
  struct mesh_vert *V;

  for (int n=0; n<sim->num_rows; ++n) {
    DL_FOREACH(sim->rows[n].verts, V) {
      V->aux[0].vector_potential -= V->Efield * dt;
    }
  }
}



void gusto_add_source_terms(struct gusto_sim *sim, double dt)
{
  struct mesh_cell *C;
  double Udot[5] = {0,0,0,0,0};
  for (int n=0; n<sim->num_rows; ++n) {
    DL_FOREACH(sim->rows[n].cells, C) {

      switch (sim->user.coordinates) {
      case 'c': gusto_cylindrical_source_terms(&C->aux[0], Udot); break;
      case 's': gusto_spherical_source_terms(&C->aux[0], Udot); break;
      }

      for (int q=0; q<5; ++q) {
	C->U[q] += Udot[q] * C->dA[0] * dt;
      }
    }
  }
}



void gusto_recover_variables(struct gusto_sim *sim)
{
  struct mesh_cell *C;
  for (int n=0; n<sim->num_rows; ++n) {
    DL_FOREACH(sim->rows[n].cells, C) {
      gusto_from_conserved(C->aux, C->U, C->dA);
    }
  }
}



void gusto_compute_vertex_velocities(struct gusto_sim *sim)
{
  struct mesh_vert *V;
  for (int n=0; n<sim->num_rows; ++n) {
    DL_FOREACH(sim->rows[n].verts, V) {
      double *u = V->aux[0].velocity_four_vector;
      V->v[1] = u[1] / u[0];
      V->v[2] = u[2] / u[0];
      V->v[3] = u[3] / u[0];
    }
  }
}



void gusto_compute_variables_at_vertices(struct gusto_sim *sim)
{
  struct mesh_vert *V;
  struct mesh_cell *C;

  /* IMPROVE: This function needs some work. Presently it sets the velocity of
     any vertex twice (except for the first and last of the row). It also would
     allow, formally, for two vertices at the same location to move in different
     directions. */

  for (int n=0; n<sim->num_rows; ++n) {
    DL_FOREACH(sim->rows[n].cells, C) {

      for (int v=0; v<4; ++v) {
	struct aux_variables *A = &C->verts[v]->aux[0];

	for (int d=1; d<4; ++d) {
	  A->velocity_four_vector[d] = C->aux[0].velocity_four_vector[d];
	  A->magnetic_four_vector[d] = C->aux[0].magnetic_four_vector[d];
	}

	A->comoving_mass_density = C->aux[0].comoving_mass_density;
	A->gas_pressure = C->aux[0].gas_pressure;
      }
    }

    DL_FOREACH(sim->rows[n].verts, V) {
      gusto_complete_aux(&V->aux[0]);
    }
  }
}



void gusto_default_aux(struct aux_variables *A)
{
  A->comoving_mass_density = 1.0;
  A->gas_pressure = 1.0;
  A->vector_potential = 0.0;
  A->velocity_four_vector[1] = 0.0;
  A->velocity_four_vector[2] = 0.0;
  A->velocity_four_vector[3] = 0.0;
  A->magnetic_four_vector[1] = 0.0;
  A->magnetic_four_vector[2] = 0.0;
  A->magnetic_four_vector[3] = 0.0;
}



void gusto_complete_aux(struct aux_variables *A)
/*
 * Fills out all primitive variables provided the following are already valid:
 *
 * u1, u2, u3
 * b1, b2, b3
 * comoving_mass_density, gas_pressure
 */
{
  double *u = A->velocity_four_vector;
  double *b = A->magnetic_four_vector;
  double uu3 = u[1]*u[1] + u[2]*u[2] + u[3]*u[3];
  double ub3 = u[1]*b[1] + u[2]*b[2] + u[3]*b[3];

  u[0] = sqrt(1.0 + uu3);
  b[0] = ub3 / u[0];

  double bb = b[1]*b[1] + b[2]*b[2] + b[3]*b[3] - b[0]*b[0];
  double dg = A->comoving_mass_density;
  double pg = A->gas_pressure;
  double ug = pg / (gamma_law_index - 1.0);
  double pb = 0.5 * bb;
  double ub = 0.5 * bb;
  double H0 = dg + (ug + ub) + (pg + pb);

  A->momentum_density[0] = H0 * u[0] * u[0] - b[0] * b[0];
  A->momentum_density[1] = H0 * u[0] * u[1] - b[0] * b[1];
  A->momentum_density[2] = H0 * u[0] * u[2] - b[0] * b[2];
  A->momentum_density[3] = H0 * u[0] * u[3] - b[0] * b[3];
  A->magnetic_pressure = pb;
  A->enthalpy_density = H0;
}



void gusto_to_conserved(struct aux_variables *A, double U[8], double dA[4])
{
  double *b = A->magnetic_four_vector;
  double *u = A->velocity_four_vector;
  double u0 = u[0];
  double u1 = u[1];
  double u2 = u[2];
  double u3 = u[3];
  double b0 = b[0];
  double b1 = b[1];
  double b2 = b[2];
  double b3 = b[3];
  double dg = A->comoving_mass_density;
  double pg = A->gas_pressure;
  double pb = A->magnetic_pressure;
  double H0 = A->enthalpy_density;

  U[DDD] = dA[0] * (dg * u0);
  U[TAU] = dA[0] * (H0 * u0 * u0 - b0 * b0 - (pg + pb) - dg * u0);
  U[S11] = dA[0] * (H0 * u0 * u1 - b0 * b1);
  U[S22] = dA[0] * (H0 * u0 * u2 - b0 * b2);
  U[S33] = dA[0] * (H0 * u0 * u3 - b0 * b3);
  U[B11] = dA[1] * (b1 * u0 - u1 * b0);
  U[B22] = dA[2] * (b2 * u0 - u2 * b0);
  U[B33] = dA[3] * (b3 * u0 - u3 * b0);

  U[S22] *= A->R; /* R != 1 if in cylindrical coordinates */
}



int gusto_from_conserved(struct aux_variables *A, double U[8], double dA[4])
{
  double Uin[8]; /* conserved densities */
  double Pin[8]; /* primitive variables d, p, v, B */
  srmhd_c2p *c2p = srmhd_c2p_new();

  Uin[DDD] = U[DDD] / dA[0];
  Uin[TAU] = U[TAU] / dA[0];
  Uin[S11] = U[S11] / dA[0];
  Uin[S22] = U[S22] / dA[0];
  Uin[S33] = U[S33] / dA[0];
  Uin[B11] = U[B11] / dA[1];
  Uin[B22] = U[B22] / dA[2]; /* toroidal field (along phi) */
  Uin[B33] = U[B33] / dA[3]; /* poloidal field (along z) */

  Uin[S22] /= A->R;

  srmhd_c2p_set_pressure_floor(c2p, -1.0);
  srmhd_c2p_set_gamma(c2p, gamma_law_index);
  srmhd_c2p_new_state(c2p, Uin);
  srmhd_c2p_estimate_from_cons(c2p);

  int error = srmhd_c2p_solve_noble1dw(c2p, Pin);
  //int error = srmhd_c2p_solve_anton2dzw(c2p, Pin);
  double B1 = Uin[B11];
  double B2 = Uin[B22];
  double B3 = Uin[B33];
  double v1 = Pin[2];
  double v2 = Pin[3];
  double v3 = Pin[4];
  double u0 = 1.0 / sqrt(1.0 - (v1*v1 + v2*v2 + v3*v3));
  double u1 = u0 * v1;
  double u2 = u0 * v2;
  double u3 = u0 * v3;
  double b0 = B1*u1 + B2*u2 + B3*u3;
  double b1 = (B1 + b0 * u1) / u0;
  double b2 = (B2 + b0 * u2) / u0;
  double b3 = (B3 + b0 * u3) / u0;
  double dg = Pin[0];
  double pg = Pin[1];

  A->comoving_mass_density = dg;
  A->gas_pressure = pg;
  A->velocity_four_vector[1] = u1;
  A->velocity_four_vector[2] = u2;
  A->velocity_four_vector[3] = u3;
  A->magnetic_four_vector[1] = b1;
  A->magnetic_four_vector[2] = b2;
  A->magnetic_four_vector[3] = b3;
  gusto_complete_aux(A);

  srmhd_c2p_del(c2p);

  if (error != 0) {
    printf("c2p failed: %s\n", srmhd_c2p_get_error(c2p, error));
    printf("%f %f %f %f %f %f %f %f\n",
	   Uin[0], Uin[1], Uin[2], Uin[3], Uin[4], Uin[5], Uin[6], Uin[7]);
    exit(1);
  }

  return error;
}



int gusto_wavespeeds(struct aux_variables *A, double n[4], double evals[8])
{
  const double n1 = n[1];
  const double n2 = n[2];
  const double n3 = n[3];
  const double u0 = A->velocity_four_vector[0];
  const double v1 = A->velocity_four_vector[1] / u0;
  const double v2 = A->velocity_four_vector[2] / u0;
  const double v3 = A->velocity_four_vector[3] / u0;
  const double b0 = A->magnetic_four_vector[0];
  const double b1 = A->magnetic_four_vector[1];
  const double b2 = A->magnetic_four_vector[2];
  const double b3 = A->magnetic_four_vector[3];
  const double bb = b1*b1 + b2*b2 + b3*b3 - b0*b0;
  const double bn = b1*n1 + b2*n2 + b3*n3;
  const double vn = v1*n1 + v2*n2 + v3*n3;
  const double dg = A->comoving_mass_density;
  const double pg = A->gas_pressure;
  const double ug = pg / (gamma_law_index - 1.0);
  const double Hg = dg + ug + pg; /* gas enthalpy density */
  const double cs2 = gamma_law_index * pg / Hg; /* sound speed squared */
  const double C = Hg + bb; /* constant in Alfven wave expression */

  const double W2 = u0*u0;
  const double W4 = W2*W2;
  const double V2 = vn*vn; /* Vx := vn^x */
  const double V3 = vn*V2;
  const double V4 = vn*V3;
  const double K  =    Hg * (1.0/cs2-1)  * W4;
  const double L  =  -(Hg +   bb/cs2)    * W2;
  const double A4 =    K    - L          -   b0*b0;
  const double A3 = -4*K*vn + L*vn*2     + 2*b0*bn;
  const double A2 =  6*K*V2 + L*(1.0-V2) +   b0*b0 - bn*bn;
  const double A1 = -4*K*V3 - L*vn*2     - 2*b0*bn;
  const double A0 =    K*V4 + L*V2       +   bn*bn;

  double roots[4];
  int nr = solve_quartic_equation(A4, A3, A2, A1, A0, roots);

  if (nr == 4) {
    evals[0] = roots[0];
    evals[1] = (bn - sqrt(C) * vn * u0) / (b0 - sqrt(C) * u0);
    evals[2] = roots[1];
    evals[3] = vn;
    evals[4] = vn;
    evals[5] = roots[2];
    evals[6] = (bn + sqrt(C) * vn * u0) / (b0 + sqrt(C) * u0);
    evals[7] = roots[3];
    return 0;
  }
  else {
    printf("bad wavespeeds\n");
    return 1;
  }
}



int gusto_fluxes(struct aux_variables *A, double n[4], double F[8])
{
  const double u0 = A->velocity_four_vector[0];
  const double u1 = A->velocity_four_vector[1];
  const double u2 = A->velocity_four_vector[2];
  const double u3 = A->velocity_four_vector[3];
  const double b0 = A->magnetic_four_vector[0];
  const double b1 = A->magnetic_four_vector[1];
  const double b2 = A->magnetic_four_vector[2];
  const double b3 = A->magnetic_four_vector[3];
  const double S0 = A->momentum_density[0]; /* T^0_0 = tau + D + pg + pb */
  const double S1 = A->momentum_density[1];
  const double S2 = A->momentum_density[2];
  const double S3 = A->momentum_density[3];
  const double dg = A->comoving_mass_density;
  const double pg = A->gas_pressure;
  const double pb = A->magnetic_pressure;
  const double D0 = dg * u0;
  const double n1 = n[1];
  const double n2 = n[2];
  const double n3 = n[3];
  const double T0 = S0 - (D0 + pg + pb); /* tau */
  const double v1 = u1 / u0;
  const double v2 = u2 / u0;
  const double v3 = u3 / u0;
  const double B1 = b1 * u0 - b0 * u1;
  const double B2 = b2 * u0 - b0 * u2;
  const double B3 = b3 * u0 - b0 * u3;
  const double Bn = B1*n1 + B2*n2 + B3*n3;
  const double vn = v1*n1 + v2*n2 + v3*n3;

  F[DDD] = D0 * vn;
  F[TAU] = T0 * vn + (pg + pb) * vn - b0 * Bn / u0;
  F[S11] = S1 * vn + (pg + pb) * n1 - b1 * Bn / u0;
  F[S22] = S2 * vn + (pg + pb) * n2 - b2 * Bn / u0;
  F[S33] = S3 * vn + (pg + pb) * n3 - b3 * Bn / u0;
  F[B11] = B1 * vn - Bn * v1;
  F[B22] = B2 * vn - Bn * v2;
  F[B33] = B3 * vn - Bn * v3;

  F[S22] *= A->R; /* R != 1 if in cylindrical coordinates */

  return 0;
}



void gusto_cylindrical_source_terms(struct aux_variables *A, double Udot[8])
{
  double *u = A->velocity_four_vector;
  double *b = A->magnetic_four_vector;
  double pg = A->gas_pressure;
  double pb = A->magnetic_pressure;
  double H0 = A->enthalpy_density;
  Udot[S11] = (pg + pb + u[2]*u[2]*H0 - b[2]*b[2]) / A->R;
}



void gusto_spherical_source_terms(struct aux_variables *A, double Udot[8])
{
  double *u = A->velocity_four_vector;
  double *b = A->magnetic_four_vector;
  double pg = A->gas_pressure;
  double pb = A->magnetic_pressure;
  double H0 = A->enthalpy_density;

  /* This is intended to be used for 1D meshes, and is valid only on the
     equatorial plane. */
  Udot[S11] = (2*(pg + pb) + u[2]*u[2]*H0 - b[2]*b[2]) / A->R;
}



void gusto_riemann(struct aux_variables *AL,
		   struct aux_variables *AR,
		   double nhat[4],
		   double Fhat[8], double s)
{
  double lamL[8];
  double lamR[8];
  double UL[8];
  double UR[8];
  double FL[8];
  double FR[8];
  double dA_unit[4] = {1, 1, 1, 1};

  gusto_wavespeeds(AL, nhat, lamL);
  gusto_wavespeeds(AR, nhat, lamR);
  gusto_to_conserved(AL, UL, dA_unit);
  gusto_to_conserved(AR, UR, dA_unit);
  gusto_fluxes(AL, nhat, FL);
  gusto_fluxes(AR, nhat, FR);

  double Sm = gusto_min3(lamL[0], lamR[0], 0.0);
  double Sp = gusto_max3(lamL[7], lamR[7], 0.0);
  double F_hll[8];
  double U_hll[8];

  for (int q=0; q<8; ++q) {
    U_hll[q] = (Sp*UR[q] - Sm*UL[q] +       (FL[q] - FR[q])) / (Sp - Sm);
    F_hll[q] = (Sp*FL[q] - Sm*FR[q] + Sp*Sm*(UR[q] - UL[q])) / (Sp - Sm);
  }

  double F[8], U[8];

  if      (        s<=Sm) for (int q=0; q<8; ++q) U[q] = UL[q];
  else if (Sm<s && s<=Sp) for (int q=0; q<8; ++q) U[q] = U_hll[q];
  else if (Sp<s         ) for (int q=0; q<8; ++q) U[q] = UR[q];

  if      (        s<=Sm) for (int q=0; q<8; ++q) F[q] = FL[q];
  else if (Sm<s && s<=Sp) for (int q=0; q<8; ++q) F[q] = F_hll[q];
  else if (Sp<s         ) for (int q=0; q<8; ++q) F[q] = FR[q];

  for (int q=0; q<8; ++q) {
    Fhat[q] = F[q] - s * U[q];
  }
}



void bc_none(struct gusto_sim *sim)
{

}



void bc_periodic_longitudinal(struct gusto_sim *sim)
{
  /* Assumes we have exactly one guard zone in the longitudinal direction. Does
     nothing to the transverse direction. */
  for (int n=0; n<sim->num_rows; ++n) {
    struct mesh_cell *Cinner_bc = sim->rows[n].cells;
    struct mesh_cell *Couter_bc = sim->rows[n].cells;

    while (Couter_bc->next) Couter_bc = Couter_bc->next;

    struct mesh_cell *Cinner_in = Cinner_bc->next;
    struct mesh_cell *Couter_in = Couter_bc->prev;

    Cinner_bc->aux[0] = Couter_in->aux[0];
    Couter_bc->aux[0] = Cinner_in->aux[0];

    gusto_to_conserved(&Cinner_bc->aux[0], Cinner_bc->U, Cinner_bc->dA);
    gusto_to_conserved(&Couter_bc->aux[0], Couter_bc->U, Couter_bc->dA);
  }
}



void bc_periodic_all(struct gusto_sim *sim)
{
  /* Assumes we have exactly one guard row in the transverse direction. The
     guard rows must have the same length and vertex locations as their interior
     row counterpart on the other side of the domain. */

  struct mesh_cell *CL_bc; /* |+| | | | | | | */
  struct mesh_cell *CL_in; /* | |+| | | | | | */
  struct mesh_cell *CR_in; /* | | | | | |+| | */
  struct mesh_cell *CR_bc; /* | | | | | | |+| */

  int N = sim->num_rows;
  CL_bc = sim->rows[  0].cells;
  CL_in = sim->rows[  1].cells;
  CR_in = sim->rows[N-2].cells;
  CR_bc = sim->rows[N-1].cells;

  while (CL_bc && CL_in && CR_in && CR_bc) {

    CL_bc->aux[0] = CR_in->aux[0];
    CR_bc->aux[0] = CL_in->aux[0];

    gusto_to_conserved(&CL_bc->aux[0], CL_bc->U, CL_bc->dA);
    gusto_to_conserved(&CR_bc->aux[0], CR_bc->U, CR_bc->dA);

    CL_bc = CL_bc->next;
    CL_in = CL_in->next;
    CR_in = CR_in->next;
    CR_bc = CR_bc->next;
  }

  bc_periodic_longitudinal(sim);
}



void bc_inflow(struct gusto_sim *sim)
{
  for (int n=0; n<sim->num_rows; ++n) {
    struct mesh_cell *Cinner0 = sim->rows[n].cells;
    struct mesh_cell *Couter0 = sim->rows[n].cells;

    while (Couter0->next) Couter0 = Couter0->next;

    gusto_default_aux(&Cinner0->aux[0]);
    gusto_default_aux(&Couter0->aux[0]);
    sim->initial_data(&sim->user, &Couter0->aux[0], Couter0->x);
    sim->initial_data(&sim->user, &Cinner0->aux[0], Cinner0->x);

    gusto_complete_aux(&Cinner0->aux[0]);
    gusto_complete_aux(&Couter0->aux[0]);
    gusto_to_conserved(&Cinner0->aux[0], Cinner0->U, Cinner0->dA);
    gusto_to_conserved(&Couter0->aux[0], Couter0->U, Couter0->dA);
  }
}



void bc_inflow_outflow(struct gusto_sim *sim)
{
  for (int n=0; n<sim->num_rows; ++n) {
    struct mesh_cell *Cinner0 = sim->rows[n].cells;
    struct mesh_cell *Couter0 = sim->rows[n].cells;

    while (Couter0->next) Couter0 = Couter0->next;

    struct mesh_cell *Couter1 = Couter0->prev;
    Couter0->aux[0] = Couter1->aux[0];
    Couter0->aux[0].gas_pressure = sim->user.pressure1;
    Cinner0->aux[0].gas_pressure = sim->user.pressure0;
    Cinner0->aux[0].comoving_mass_density = sim->user.density0;
    Cinner0->aux[0].velocity_four_vector[1] = sim->user.fourvel0;

    gusto_complete_aux(&Cinner0->aux[0]);
    gusto_complete_aux(&Couter0->aux[0]);
    gusto_to_conserved(&Cinner0->aux[0], Cinner0->U, Cinner0->dA);
    gusto_to_conserved(&Couter0->aux[0], Couter0->U, Couter0->dA);
  }
}



OpBoundaryCon gusto_lookup_boundary_con(const char *user_key)
{
  const char *keys[] = {
    "none",
    "periodic_longitudinal",
    "periodic_all",
    "inflow",
    "inflow_outflow",
    NULL
  } ;
  OpBoundaryCon vals[] = {
    bc_none,
    bc_periodic_longitudinal,
    bc_periodic_all,
    bc_inflow,
    bc_inflow_outflow,
    NULL } ;
  int n = 0;
  const char *key;
  OpBoundaryCon val;
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
  printf("[gusto] ERROR: no such boundary_con=%s\n", user_key);
  return NULL;
}
