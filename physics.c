#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "utlist.h"
#include "gusto.h"
#include "srmhd_c2p.h"
#include "quartic.h"



void gusto_initial_data(struct gusto_sim *sim)
{
  struct mesh_cell *C;
  struct aux_variables *A;

  for (int n=0; n<sim->num_rows; ++n) {

    DL_FOREACH(sim->rows[n].cells, C) {
      A = &C->aux;
      gusto_default_aux(A);
      sim->initial_data(&sim->user, A, C->x);
      gusto_complete_aux(A);
      gusto_to_conserved(A, C->U, C->dA);
    }

  }
}



void gusto_validate_fluxes(struct gusto_sim *sim)
{
  struct mesh_cell *C;

  for (int n=0; n<sim->num_rows; ++n) {

    DL_FOREACH(sim->rows[n].cells, C) {

      double F[8];
      double nhat[4] = {0, 0, 0, 1};
      double dA = C->geom.area_element[3];

      gusto_fluxes(&C->aux, nhat, F);

      printf("F.dA = %f %f %f %f %f %f (dA=%f)\n",
	     F[0] * dA,
	     F[1] * dA,
	     F[2] * dA,
	     F[3] * dA,
	     F[4] * dA, F[B22] * C->geom.line_element[2], dA);

    }
  }
}



void gusto_compute_fluxes(struct gusto_sim *sim)
{
  struct mesh_face *F;

  for (int n=0; n<sim->num_rows; ++n) {
    DL_FOREACH(sim->rows[n].faces, F) {

      struct mesh_cell *CL = F->cells[0];
      struct mesh_cell *CR = F->cells[1];
      struct aux_variables AL, AR;
      double vpar = 0.0;

      if (CL && CR) {
	AL = CL->aux;
	AR = CR->aux;
      }
      else if (CL) {
	AL = CL->aux;
	AR = CL->aux;
      }
      else if (CR) {
	AL = CR->aux;
	AR = CR->aux;
      }

      double nhat[4] = {0, 0, 0, 1};
      gusto_riemann(&AL, &AR, nhat, F->Fhat, vpar);

    }
  }
}



void gusto_cache_rk(struct gusto_sim *sim)
{
  struct mesh_cell *C;
  struct mesh_face *F;

  for (int n=0; n<sim->num_rows; ++n) {
    DL_FOREACH(sim->rows[n].cells, C) {
      for (int q=0; q<5; ++q) {
	C->U_rk[q] = C->U[q];
      }
    }
    DL_FOREACH(sim->rows[n].faces, F) {
      for (int d=0; d<4; ++d) {
	F->x_rk[d] = F->x[d];
      }
    }
  }
}



void gusto_average_rk(struct gusto_sim *sim, double b)
{
  struct mesh_cell *C;
  struct mesh_face *F;

  for (int n=0; n<sim->num_rows; ++n) {
    DL_FOREACH(sim->rows[n].cells, C) {
      for (int q=0; q<5; ++q) {
	C->U[q] = b * C->U[q] + (1 - b) * C->U_rk[q];
      }
    }
    DL_FOREACH(sim->rows[n].faces, F) {
      for (int d=0; d<4; ++d) {
	F->x[d] = b * F->x[d] + (1 - b) * F->x_rk[d];
      }
    }
  }
}



void gusto_transmit_fluxes(struct gusto_sim *sim, double dt)
{
  struct mesh_face *F;

  for (int n=0; n<sim->num_rows; ++n) {
    DL_FOREACH(sim->rows[n].faces, F) {
      struct mesh_cell *CL = F->cells[0];
      struct mesh_cell *CR = F->cells[1];

      for (int q=0; q<5; ++q) {
	if (CL) CL->U[q] -= F->Fhat[q] * F->dA[0] * dt;
	if (CR) CR->U[q] += F->Fhat[q] * F->dA[0] * dt;
      }

      if (CL) CL->U[B22] -= F->Fhat[B22] * F->dA[2] * dt;
      if (CR) CR->U[B22] += F->Fhat[B22] * F->dA[2] * dt;
    }
  }
}



void gusto_add_source_terms(struct gusto_sim *sim, double dt)
{
  struct mesh_cell *C;

  for (int n=0; n<sim->num_rows; ++n) {
    DL_FOREACH(sim->rows[n].cells, C) {

      double Udot[5] = {0,0,0,0,0};
      gusto_geometry_source_terms(&C->aux, &C->geom, Udot);

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
      if (C->cell_type == 'g') {
	continue;
      }
      if (gusto_from_conserved(&C->aux, C->U, C->dA)) {
	printf("[gusto] at position (%f -> %f)\n", C->faces[0]->x[3], C->faces[1]->x[3]);
	exit(1);
      }
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
  gusto_complete_aux(A);
}



void gusto_complete_aux(struct aux_variables *A)
/*
 * Fill out all primitive variables provided the following are already valid:
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
  Uin[B22] = U[B22] / dA[2];
  Uin[B33] = U[B33] / dA[3];

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

  if (error != 0) {
    printf("[gusto] ERROR: %s\n", srmhd_c2p_get_error(c2p, error));
    printf("[gusto] failed on U = [%f %f %f %f %f %f %f %f]\n",
	   Uin[0], Uin[1], Uin[2], Uin[3],
	   Uin[4], Uin[5], Uin[6], Uin[7]);
  }

  srmhd_c2p_del(c2p);

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



void gusto_electric_field(struct aux_variables *A, double E[4])
{
  double *u = A->velocity_four_vector;
  double *b = A->magnetic_four_vector;
  double B[4] = {0,
		 b[1] * u[0] - b[0] * u[1],
		 b[2] * u[0] - b[0] * u[2],
		 b[3] * u[0] - b[0] * u[3]};
  double v[4] = {0, u[1]/u[0], u[2]/u[0], u[3]/u[0]};
  double F[4] = VEC4_CROSS(B, v);
  E[1] = F[1];
  E[2] = F[2];
  E[3] = F[3];
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



void bc_inflow(struct gusto_sim *sim)
{

  for (int n=0; n<sim->num_rows; ++n) {

    struct mesh_cell *C0 = sim->rows[n].cells;
    struct mesh_cell *C1 = sim->rows[n].cells->prev;

    gusto_default_aux(&C0->aux);
    sim->initial_data(&sim->user, &C0->aux, C0->x);
    gusto_complete_aux(&C0->aux);
    gusto_to_conserved(&C0->aux, C0->U, C0->dA);

    gusto_default_aux(&C1->aux);
    sim->initial_data(&sim->user, &C1->aux, C1->x);
    gusto_complete_aux(&C1->aux);
    gusto_to_conserved(&C1->aux, C1->U, C1->dA);

  }
}



void bc_inflow_outflow(struct gusto_sim *sim)
{

  for (int n=0; n<sim->num_rows; ++n) {

    struct mesh_cell *C0 = sim->rows[n].cells;
    struct mesh_cell *C1 = sim->rows[n].cells->prev;

    gusto_default_aux(&C0->aux);
    sim->initial_data(&sim->user, &C0->aux, C0->x);
    gusto_complete_aux(&C0->aux);
    gusto_to_conserved(&C0->aux, C0->U, C0->dA);

    double tmpR = C1->aux.R;
    gusto_default_aux(&C1->aux);
    C1->aux = C1->prev->aux;
    C1->aux.R = tmpR;
    gusto_complete_aux(&C1->aux);
    gusto_to_conserved(&C1->aux, C1->U, C1->dA);

  }
}



OpBoundaryCon gusto_lookup_boundary_con(const char *user_key)
{
  const char *keys[] = {
    "none",
    "inflow",
    "inflow_outflow",
    NULL
  } ;
  OpBoundaryCon vals[] = {
    bc_none,
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
