#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gusto.h"
#include "srmhd_c2p.h"
#include "quartic.h"


#include <sys/time.h>

void *gusto_start_clock()
{
  struct timeval *T = (struct timeval *) malloc(sizeof(struct timeval));
  gettimeofday(T, NULL);
  return T;
}

double gusto_stop_clock(void *clock_s)
{
  struct timeval *T0 = (struct timeval *) clock_s;
  struct timeval *T1 = (struct timeval *) malloc(sizeof(struct timeval));
  gettimeofday(T1, NULL);
  double ds = T1->tv_sec  - T0->tv_sec;
  double du = T1->tv_usec - T0->tv_usec;
  double dt = ds + 1e-6*du;
  free(T0);
  free(T1);
  return dt;
}



void gusto_vars_complete_aux(struct aux_variables *A);


void initial_data_function(struct aux_variables *A, double *X)
{
  A->velocity_four_vector[1] = 0.0;
  A->velocity_four_vector[2] = 0.0;
  A->velocity_four_vector[3] = 0.0;
  A->magnetic_four_vector[1] = 0.0;
  A->magnetic_four_vector[2] = 0.0;
  A->magnetic_four_vector[3] = 0.0;
  double x = X[1] - 0.5;
  double y = X[3] - 0.5;
  if (sqrt(x*x + y*y) < 0.125) {
    A->comoving_mass_density = 1.0;
    A->gas_pressure = 1.0;
  }
  else {
    A->comoving_mass_density = 0.1;
    A->gas_pressure = 1.0;//0.125;
  }
  gusto_vars_complete_aux(A);
}



void gusto_vars_complete_aux(struct aux_variables *A)
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
  b[0] = u[0] * ub3 / (1 + u[0] * uu3);

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
}



void gusto_vars_to_conserved(struct aux_variables *A, double U[8], double dA[4])
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
  double bb = b1*b1 + b2*b2 + b3*b3 - b0*b0;
  double dg = A->comoving_mass_density;
  double pg = A->gas_pressure;
  double ug = pg / (gamma_law_index - 1.0);
  double pb = 0.5 * bb;
  double ub = 0.5 * bb;
  double H0 = dg + (ug + ub) + (pg + pb);
  U[DDD] = dA[0] * (dg * u0);
  U[TAU] = dA[0] * (H0 * u0 * u0 - b0 * b0 - (pg + pb) - dg * u0);
  U[S11] = dA[0] * (H0 * u0 * u1 - b0 * b1);
  U[S22] = dA[0] * (H0 * u0 * u2 - b0 * b2);
  U[S33] = dA[0] * (H0 * u0 * u3 - b0 * b3);
  U[B11] = dA[1] * (b1 * u0 - u1 * b0);
  U[B22] = dA[2] * (b2 * u0 - u2 * b0);
  U[B33] = dA[3] * (b3 * u0 - u3 * b0);
}



int gusto_vars_from_conserved(struct aux_variables *A, double U[8], double dA[4])
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

  srmhd_c2p_set_pressure_floor(c2p, -1.0);
  srmhd_c2p_set_gamma(c2p, gamma_law_index);
  srmhd_c2p_new_state(c2p, Uin);
  srmhd_c2p_estimate_from_cons(c2p);

  int error = srmhd_c2p_solve_noble1dw(c2p, Pin);
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
  gusto_vars_complete_aux(A);

  srmhd_c2p_del(c2p);

  if (error != 0) {
    printf("c2p failed: %s\n", srmhd_c2p_get_error(c2p, error));
    printf("%f %f %f %f %f %f %f %f\n", Uin[0], Uin[1], Uin[2], Uin[3], Uin[4], Uin[5], Uin[6], Uin[7]);
  }
  
  return error;
}



void gusto_generate_verts(struct gusto_sim *sim)
{
  sim->num_rows = sim->user.N[0];
  sim->row_size = (int *) malloc(sim->num_rows * sizeof(int));
  sim->verts = (struct mesh_vert **)
    malloc(sim->num_rows * sizeof(struct mesh_vert*));

  for (int n=0; n<sim->num_rows; ++n) {
    sim->row_size[n] = sim->user.N[1];

    sim->verts[n] = (struct mesh_vert *)
      malloc(sim->row_size[n] * sizeof(struct mesh_vert));

    for (int i=0; i<sim->row_size[n]; ++i) {

      double R0 = 0.0;
      double R1 = 1.0;
      double z0 = 0.0;
      double z1 = 1.0;
      double dR = (R1 - R0) / (sim->num_rows - 1);
      double dz = (z1 - z0) / (sim->row_size[n] - 1);

      struct mesh_vert p;
      p.x[0] = 0.0;
      p.x[1] = R0 + n * dR + 0.05 * sin(4 * M_PI * (z0 + i * dz));
      p.x[2] = 0.0;
      p.x[3] = z0 + i * dz;

      p.u[0] = 0.0;
      p.u[1] = 0.0;
      p.u[2] = 0.0;
      p.u[3] = 0.0;
      p.cell = NULL;
      
      sim->verts[n][i] = p;
    }
  }
}



#define VEC4_SUB(x,y) {0,x[1]-y[1],x[2]-y[2],x[3]-y[3]}
#define VEC4_ADD(x,y) {0,x[1]+y[1],x[2]+y[2],x[3]+y[3]}
#define VEC4_DOT(x,y) (x[1]*y[1]+x[2]*y[2]+x[3]*y[3])
#define VEC4_MOD(x) sqrt(VEC4_DOT(x,x))
#define VEC4_CROSS(x,y) {0,			\
			 x[2]*y[3]-x[3]*y[2],	\
			 x[3]*y[1]-x[1]*y[3],	\
			 x[1]*y[2]-x[2]*y[1]}
#define VEC4_NORMALIZE(x) do {					\
    double norm = sqrt(x[1]*x[1] + x[2]*x[2] + x[3]*x[3]);	\
    x[1] /= norm;						\
    x[2] /= norm;						\
    x[3] /= norm;						\
  } while (0)


void gusto_generate_cells(struct gusto_sim *sim)
{
  ARRAY_INIT(struct mesh_cell, sim->cells, sim->num_cells, sim->num_cells_max);

  sim->smallest_cell_length = 1.0;
  
  for (int n=0; n<sim->num_rows-1; ++n) {
    for (int i=0; i<sim->row_size[n]-1; ++i) {

      struct mesh_cell C;

      C.verts[0] = &sim->verts[n+0][i+0];
      C.verts[1] = &sim->verts[n+0][i+1];
      C.verts[2] = &sim->verts[n+1][i+0];
      C.verts[3] = &sim->verts[n+1][i+1];

      /*
       * These vectors define the 2-forms on the cell. There are four, one for
       * each corner.
      */
      double dR0[4] = VEC4_SUB(C.verts[2]->x, C.verts[0]->x);
      double dR1[4] = VEC4_SUB(C.verts[3]->x, C.verts[1]->x);
      double dz0[4] = VEC4_SUB(C.verts[1]->x, C.verts[0]->x);
      double dz1[4] = VEC4_SUB(C.verts[3]->x, C.verts[2]->x);
      double dphi[4] = {0, 0, 1, 0};
      double dAf0[4] = VEC4_CROSS(dz0, dR0);
      double dAf1[4] = VEC4_CROSS(dz0, dR1);
      double dAf2[4] = VEC4_CROSS(dz1, dR0);
      double dAf3[4] = VEC4_CROSS(dz1, dR1);
      double dAR0[4] = VEC4_CROSS(dphi, dz0);
      double dAR1[4] = VEC4_CROSS(dphi, dz1);
      double dAz0[4] = VEC4_CROSS(dR0, dphi);
      double dAz1[4] = VEC4_CROSS(dR1, dphi);

      C.dA[1] = 0.50 * (VEC4_MOD(dAR0) + VEC4_MOD(dAR1));
      C.dA[2] = 0.25 * (VEC4_MOD(dAf0) + VEC4_MOD(dAf1) +
			VEC4_MOD(dAf2) + VEC4_MOD(dAf3));
      C.dA[3] = 0.50 * (VEC4_MOD(dAz0) + VEC4_MOD(dAz1));
      C.dA[0] = C.dA[2]; /* Volume and phi cross section are the same */

      double dAR = gusto_min3(VEC4_MOD(dAR0), VEC4_MOD(dAR1), 1.0);
      double dAz = gusto_min3(VEC4_MOD(dAz0), VEC4_MOD(dAz1), 1.0);

      if (dAR < sim->smallest_cell_length) sim->smallest_cell_length = dAR;
      if (dAz < sim->smallest_cell_length) sim->smallest_cell_length = dAz;

      ARRAY_APPEND(struct mesh_cell,
		   sim->cells, sim->num_cells, sim->num_cells_max, C);
    }
  }

  int j = 0;
  for (int n=0; n<sim->num_rows-1; ++n) {
    for (int i=0; i<sim->row_size[n]-1; ++i) {
      sim->verts[n][i].cell = &sim->cells[j];
      j += 1;
    }
  }
 
  printf("num_cells: %d\n", sim->num_cells);
  printf("num_cells_max: %d\n", sim->num_cells_max);
}



void gusto_generate_faces(struct gusto_sim *sim)
{
  ARRAY_INIT(struct mesh_face, sim->faces, sim->num_faces, sim->num_faces_max);

  /* R faces */
  for (int n=0; n<sim->num_rows; ++n) {
    for (int i=0; i<sim->row_size[n]-1; ++i) {

      struct mesh_face F;
      F.verts[0] = &sim->verts[n][i+0];
      F.verts[1] = &sim->verts[n][i+1];
      F.cells[0] = n == 0 ? NULL : sim->verts[n-1][i].cell;
      F.cells[1] =                 sim->verts[n+0][i].cell; /* could be NULL */
      
      double df[4] = {0, 0, 1, 0};
      double dz[4] = VEC4_SUB(F.verts[1]->x, F.verts[0]->x);
      double dA[4] = VEC4_CROSS(df, dz);

      F.nhat[0] = VEC4_MOD(dA);
      F.nhat[1] = dA[1] / F.nhat[0];
      F.nhat[2] = dA[2] / F.nhat[0];
      F.nhat[3] = dA[3] / F.nhat[0];
	 
      ARRAY_APPEND(struct mesh_face,
		   sim->faces, sim->num_faces, sim->num_faces_max, F);
    }
  }

  /* z faces */
  for (int n=0; n<sim->num_rows-1; ++n) {
    for (int i=0; i<sim->row_size[n]; ++i) {

      struct mesh_face F;
      F.verts[0] = &sim->verts[n+0][i];
      F.verts[1] = &sim->verts[n+1][i];
      F.cells[0] = i == 0 ? NULL : sim->verts[n][i-1].cell;
      F.cells[1] =                 sim->verts[n][i+0].cell; /* could be NULL */

      double df[4] = {0, 0, 1, 0};
      double dR[4] = VEC4_SUB(F.verts[1]->x, F.verts[0]->x);
      double dA[4] = VEC4_CROSS(dR, df);

      F.nhat[0] = VEC4_MOD(dA);
      F.nhat[1] = dA[1] / F.nhat[0];
      F.nhat[2] = dA[2] / F.nhat[0];
      F.nhat[3] = dA[3] / F.nhat[0];

      ARRAY_APPEND(struct mesh_face,
		   sim->faces, sim->num_faces, sim->num_faces_max, F);
    }
  }
  
  printf("num_faces: %d\n", sim->num_faces);
  printf("num_faces_max: %d\n", sim->num_faces_max);
}



void gusto_free(struct gusto_sim *sim)
{
  for (int n=0; n<sim->num_rows; ++n) {
    free(sim->verts[n]);
  }
  free(sim->verts);
  free(sim->row_size);
  ARRAY_FREE(struct mesh_cell, sim->cells, sim->num_cells, sim->num_cells_max);
  ARRAY_FREE(struct mesh_face, sim->faces, sim->num_faces, sim->num_faces_max);
}



void gusto_initial_data(struct gusto_sim *sim)
{
  /* set initial data on cell volumes */
  for (int j=0; j<sim->num_cells; ++j) {
    double *X0 = sim->cells[j].verts[0]->x;
    double *X1 = sim->cells[j].verts[1]->x;
    double *X2 = sim->cells[j].verts[2]->x;
    double *X3 = sim->cells[j].verts[3]->x;   
    double R = 0.25 * (X0[1] + X1[1] + X2[1] + X3[1]);
    double z = 0.25 * (X0[3] + X1[3] + X2[3] + X3[3]);
    double X[4] = {0, R, 0, z};
    struct aux_variables *A = &sim->cells[j].aux[0];
    initial_data_function(A, X);
    gusto_vars_to_conserved(A, sim->cells[j].U, sim->cells[j].dA);
  }

  /* set initial data on vertices */
  for (int n=0; n<sim->num_rows; ++n) {
    for (int i=0; i<sim->row_size[n]; ++i) {
      double *X = sim->verts[n][i].x;
      struct aux_variables *A = &sim->verts[n][i].aux[0];
      initial_data_function(A, X);
    }
  }
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
  const double d0 = A->comoving_mass_density;
  const double pg = A->gas_pressure;
  const double pb = A->magnetic_pressure;
  const double D0 = d0 * u0;
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
  return 0;
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
  gusto_vars_to_conserved(AL, UL, dA_unit);
  gusto_vars_to_conserved(AR, UR, dA_unit);
  gusto_fluxes(AL, nhat, FL);
  gusto_fluxes(AR, nhat, FR);

  double F_hll[8];
  double U_hll[8];
  double Sm = gusto_min3(lamL[0], lamR[0], 0.0);
  double Sp = gusto_max3(lamL[7], lamR[7], 0.0);

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



void gusto_compute_fluxes(struct gusto_sim *sim)
{
  for (int j=0; j<sim->num_faces; ++j) {

    struct mesh_face *face = &sim->faces[j];
    struct mesh_cell *CL = face->cells[0];
    struct mesh_cell *CR = face->cells[1];

    double vpar = 0.0 * face->nhat[1] + 0.0 * face->nhat[3];

    if (CL && CR) {
      gusto_riemann(CL->aux, CR->aux, face->nhat, face->Fhat, vpar);
    }
    else if (CL) {
      gusto_riemann(CL->aux, CL->aux, face->nhat, face->Fhat, vpar);
    }
    else if (CR) {
      gusto_riemann(CR->aux, CR->aux, face->nhat, face->Fhat, vpar);
    }
  }
}



void gusto_transmit_fluxes(struct gusto_sim *sim, double dt)
{
  for (int j=0; j<sim->num_faces; ++j) {
    struct mesh_face *face = &sim->faces[j];
    struct mesh_cell *CL = face->cells[0];
    struct mesh_cell *CR = face->cells[1];
    for (int q=0; q<8; ++q) {
      if (CL) CL->U[q] -= face->Fhat[q] * face->nhat[0] * dt;
      if (CR) CR->U[q] += face->Fhat[q] * face->nhat[0] * dt;
    }    
  }
}



void gusto_recover_variables(struct gusto_sim *sim)
{
  for (int j=0; j<sim->num_cells; ++j) {
    struct mesh_cell *C = &sim->cells[j];
    gusto_vars_from_conserved(C->aux, C->U, C->dA);
  }
}



void gusto_compute_variables_at_vertices(struct gusto_sim *sim)
{
  for (int n=0; n<sim->num_rows; ++n) {
    for (int i=0; i<sim->row_size[n]; ++i) {
      for (int d=1; d<4; ++d) {
	sim->verts[n][i].aux[0].velocity_four_vector[d] = 0.0;
	sim->verts[n][i].aux[0].magnetic_four_vector[d] = 0.0;
      }
      sim->verts[n][i].aux[0].comoving_mass_density = 0.0;
      sim->verts[n][i].aux[0].gas_pressure = 0.0;
    }
  }
  
  for (int j=0; j<sim->num_cells; ++j) {
    struct mesh_cell *C = &sim->cells[j];
    for (int v=0; v<4; ++v) {
      struct aux_variables *A = &C->verts[v]->aux[0];
      for (int d=1; d<4; ++d) {
	A->velocity_four_vector[d] += C->aux[0].velocity_four_vector[d];
	A->magnetic_four_vector[d] += C->aux[0].magnetic_four_vector[d];
      }
      A->comoving_mass_density += C->aux[0].comoving_mass_density;
      A->gas_pressure += C->aux[0].gas_pressure;
    }
  }

  for (int n=0; n<sim->num_rows; ++n) {
    for (int i=0; i<sim->row_size[n]; ++i) {
      for (int d=1; d<4; ++d) {
	sim->verts[n][i].aux[0].velocity_four_vector[d] /= 4;
	sim->verts[n][i].aux[0].magnetic_four_vector[d] /= 4;
      }
      sim->verts[n][i].aux[0].comoving_mass_density /= 4;
      sim->verts[n][i].aux[0].gas_pressure /= 4;
      gusto_vars_complete_aux(&sim->verts[n][i].aux[0]);
    }
  }
}



void gusto_write_to_ascii(struct gusto_sim *sim)
{
  FILE *outf = fopen("gusto.dat", "w");
  for (int j=0; j<sim->num_cells; ++j) {
    struct mesh_cell *C = &sim->cells[j];
    double *u = C->aux[0].velocity_four_vector;
    double *b = C->aux[0].magnetic_four_vector;
    double z = 0.25 * (C->verts[0]->x[3] + C->verts[1]->x[3] +
		       C->verts[2]->x[3] + C->verts[3]->x[3]);
    fprintf(outf, "%f "
	    "%f %f %f %f "
	    "%f %f %f %f "
	    "%f %f\n", z,
	    u[0], u[1], u[2], u[3],
	    b[0], b[1], b[2], b[3],
	    C->aux[0].comoving_mass_density,
	    C->aux[0].gas_pressure);
  }
  fclose(outf);
}



/*
 * Main function
 * =====================================================================
 */
int main(int argc, char **argv)
{
  struct gusto_sim sim;
  gusto_user_set_defaults(&sim.user);
  gusto_status_set_defaults(&sim.status);

  
  for (int n=1; n<argc; ++n) {
    gusto_user_set_from_arg(&sim.user, argv[n]);
  }


  gusto_user_report(&sim.user);
  gusto_generate_verts(&sim);
  gusto_generate_cells(&sim);
  gusto_generate_faces(&sim);
  gusto_initial_data(&sim);


  double dt = 0.25 * sim.smallest_cell_length;
  sim.status.time_step = dt;

  while (sim.status.time_simulation < sim.user.tmax) {

    void *start_cycle = gusto_start_clock();

    gusto_compute_fluxes(&sim);
    gusto_transmit_fluxes(&sim, dt);
    gusto_recover_variables(&sim);

    double seconds = gusto_stop_clock(start_cycle);

    sim.status.kzps = 1e-3 * sim.num_cells / seconds;

    if (sim.status.iteration % 1 == 0) {
      printf("[ffe] n=%06d t=%6.4e dt=%6.4e %3.2f kzps\n",
	     sim.status.iteration,
	     sim.status.time_simulation,
	     sim.status.time_step,
	     sim.status.kzps);
    }
    
    sim.status.time_simulation += dt;
    sim.status.iteration += 1;
  }

  gusto_compute_variables_at_vertices(&sim);
  gusto_write_to_ascii(&sim);
  gusto_write_checkpoint(&sim, "chkpt.0000.h5");
  gusto_read_write_status(&sim.status, "chkpt.0000.h5", 'a');
  gusto_read_write_user(&sim.user, "chkpt.0000.h5", 'a');

  gusto_free(&sim);
  return 0;
}
