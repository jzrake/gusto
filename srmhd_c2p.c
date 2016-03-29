/*------------------------------------------------------------------------------
 * FILE: srmhd_c2p.c
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * REFERENCES:
 *
 * Anton, L., Zanotti, O., Miralles, J. A., Martı, J. M., Ibanez, J. M., Font,
 * J. A., & Pons, J. A. 2006, The Astrophysical Journal, 637, 296
 *
 * Noble, S. C., Gammie, C. F., McKinney, J. C., & Zanna, L. D.  2006, The
 * Astrophysical Journal, 641, 626
 *
 *
 * DESCRIPTION:
 *
 * This module solves for the primitive state (rho, pre, v, B) given the
 * conserved state (ddd, tau, S, B). In other words it solves for the 5 unknowns
 * rho, pre, vx, vy, vz. An adiabatic equation of state is assumed throughout,
 * with the index being provided by srmhd_c2p_set_gamma(). The solve functions
 * promise not to modify the pointer to result primitives unless the execution
 * is successful.
 *
 * ------------------------------------------------------------------------------
 */


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "srmhd_c2p.h"


enum { ddd, tau, Sx, Sy, Sz, Bx, By, Bz }; // Conserved
enum { rho, pre, vx, vy, vz };             // Primitive

struct srmhd_c2p {
  int MaxIterations;
  double Tolerance;
  double bigZ;
  double bigW;
  double smlZ;
  double smlW;

  int AppliedPressureFloor;
  int Iterations;
  double AdiabaticGamma;
  double gamf;
  double D,Tau;
  double S2,B2,BS,BS2;
  double Cons[8];
  double Z_start, W_start;
  double PressureFloor;
} ;


srmhd_c2p *srmhd_c2p_new()
{
  srmhd_c2p *c2p = (srmhd_c2p*) malloc(sizeof(srmhd_c2p));
  c2p->MaxIterations = 250;
  c2p->Tolerance = 1e-12;
  c2p->bigZ = 1e20;
  c2p->bigW = 1e12;
  c2p->smlZ = 0.0;
  c2p->smlW = 1.0;
  c2p->Iterations = 0;
  c2p->AppliedPressureFloor = 0;
  c2p->AdiabaticGamma = 1.4;
  c2p->PressureFloor = -1.0;
  return c2p;
}
void srmhd_c2p_del(srmhd_c2p *c2p)
{
  free(c2p);
}

int srmhd_c2p_put_pressure_floor(srmhd_c2p *c2p)
{
  return c2p->AppliedPressureFloor;
}

// Set method for the adiabatic index (defaults to 1.4)
// -----------------------------------------------------------------------------
void srmhd_c2p_set_gamma(srmhd_c2p *c2p, double adiabatic_gamma)
{
  c2p->AdiabaticGamma = adiabatic_gamma;
  c2p->gamf = (c2p->AdiabaticGamma - 1.0) / c2p->AdiabaticGamma;
}

// Get method for the adiabatic index (defaults to 1.4)
// -----------------------------------------------------------------------------
void srmhd_c2p_get_gamma(srmhd_c2p *c2p, double *adiabatic_gamma)
{
  *adiabatic_gamma = c2p->AdiabaticGamma;
}

// Get method for the iterations on the last execution
// -----------------------------------------------------------------------------
int srmhd_c2p_get_iterations(srmhd_c2p *c2p)
{
  return c2p->Iterations;
}

int srmhd_c2p_set_pressure_floor(srmhd_c2p *c2p, double pf)
{
  c2p->PressureFloor = pf;
  return pf < 0.0;
}

// Provide a new conserved state in memory. Before the solver is executed, the
// user must provide a guess for the initial primitive state, by calling either
// estimate_from_cons() or set_starting_prim(P).
// -----------------------------------------------------------------------------
void srmhd_c2p_new_state(srmhd_c2p *c2p, const double *U)
{
  c2p->D    = U[ddd];
  c2p->Tau  = U[tau];
  c2p->S2   = U[Sx]*U[Sx] + U[Sy]*U[Sy] + U[Sz]*U[Sz];
  c2p->B2   = U[Bx]*U[Bx] + U[By]*U[By] + U[Bz]*U[Bz];
  c2p->BS   = U[Bx]*U[Sx] + U[By]*U[Sy] + U[Bz]*U[Sz];
  c2p->BS2  = c2p->BS * c2p->BS;
  memcpy(c2p->Cons, U, 8*sizeof(double));
}

// This estimate becomes exact for no magnetic field in the NR limit.
// -----------------------------------------------------------------------------
void srmhd_c2p_estimate_from_cons(srmhd_c2p *c2p)
{
  c2p->Z_start = sqrt(c2p->S2 + c2p->D*c2p->D);
  c2p->W_start = c2p->Z_start / c2p->D;
}
void srmhd_c2p_get_starting_prim(srmhd_c2p *c2p, double *P)
{
  srmhd_c2p_reconstruct_prim(c2p, c2p->Z_start, c2p->W_start, P);
}

// Explicitly provide a guess to be used at the initial iteration.
// -----------------------------------------------------------------------------
void srmhd_c2p_set_starting_prim(srmhd_c2p *c2p, const double *P)
{
  const double V2 = P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz];
  const double W2 = 1.0 / (1.0 - V2);
  const double e = P[pre] / (P[rho] * (c2p->AdiabaticGamma - 1.0));
  const double h = 1.0 + e + P[pre] / P[rho];

  c2p->Z_start = P[rho] * h * W2;
  c2p->W_start = sqrt(W2);
}

// Using Z=rho*h*W^2, and W, get the primitive variables.
// -----------------------------------------------------------------------------
int srmhd_c2p_reconstruct_prim(srmhd_c2p *c2p, double Z, double W, double *Pout)
{
  double P[8]; // Place result into temporary prim state for now.
  const double b0    = c2p->BS * W / Z;
  const double D     = c2p->D;
  const double B2    = c2p->B2;
  const double *Cons = c2p->Cons;

  P[rho] =  D/W;
  P[pre] = (D/W) * (Z/(D*W) - 1.0) * c2p->gamf;
  P[vx ] = (Cons[Sx] + b0*Cons[Bx]/W) / (Z+B2);
  P[vy ] = (Cons[Sy] + b0*Cons[By]/W) / (Z+B2);
  P[vz ] = (Cons[Sz] + b0*Cons[Bz]/W) / (Z+B2);
  P[Bx ] =  Cons[Bx];
  P[By ] =  Cons[By];
  P[Bz ] =  Cons[Bz];

  if (P[pre] <= 0.0 && c2p->PressureFloor > 0.0) {
    /* printf("setting pressure floor, p=%4.3e -> %2.1e\n", */
    /*     P[pre], PressureFloor); */
    P[pre] = c2p->PressureFloor;
    c2p->AppliedPressureFloor = 1;
  }

  int error = srmhd_c2p_check_prim(c2p, P);
  if (error) {
    return error;
  }

  // Don't actually modify the user's prim state until we're sure the state is
  // good.
  // ---------------------------------------------------------------------------
  memcpy(Pout, P, 8*sizeof(double));
  return SRMHD_C2P_SUCCESS;
}

// Routines to verify the health of conserved or primitive states
// -----------------------------------------------------------------------------
int srmhd_c2p_check_cons(srmhd_c2p *c2p, const double *U)
{
  int i;
  if (U[ddd] < 0.0) return SRMHD_C2P_CONS_NEGATIVE_DENSITY;
  if (U[tau] < 0.0) return SRMHD_C2P_CONS_NEGATIVE_ENERGY;
  for (i=0; i<8; ++i) {
    if (isnan(U[i])) {
      return SRMHD_C2P_CONS_CONTAINS_NAN;
    }
  }
  return SRMHD_C2P_SUCCESS;
}
int srmhd_c2p_check_prim(srmhd_c2p *c2p, const double *P)
{
  int i;
  const double v2 = P[vx]*P[vx] + P[vy]*P[vy] + P[vz]*P[vz];
  if (v2   >=  1.0) return SRMHD_C2P_PRIM_SUPERLUMINAL;
  if (P[pre] < 0.0) return SRMHD_C2P_PRIM_NEGATIVE_PRESSURE;
  if (P[rho] < 0.0) return SRMHD_C2P_PRIM_NEGATIVE_RESTMASS;
  for (i=0; i<8; ++i) {
    if (isnan(P[i])) {
      return SRMHD_C2P_PRIM_CONTAINS_NAN;
    }
  }
  return SRMHD_C2P_SUCCESS;
}

// Solution based on Anton & Zanotti (2006), equations 84 and 85.
// -----------------------------------------------------------------------------
int srmhd_c2p_solve_anton2dzw(srmhd_c2p *c2p, double *P)
{
  int bad_input = srmhd_c2p_check_cons(c2p, c2p->Cons);
  if (bad_input) {
    return bad_input;
  }
  // Starting values
  // ---------------------------------------------------------------------------
  c2p->Iterations = 0;
  const double D   = c2p->D;
  const double B2  = c2p->B2;
  const double S2  = c2p->S2;
  const double BS2 = c2p->BS2;
  const double Tau = c2p->Tau;

  double error = 1.0;
  double W = c2p->W_start;
  double Z = c2p->Z_start;

  while (error > c2p->Tolerance) {

    const double Z2 = Z*Z;
    const double Z3 = Z*Z2;
    const double W2 = W*W;
    const double W3 = W*W2;
    const double Pre = (D/W) * (Z/(D*W) - 1.0) * c2p->gamf;

    const double df0dZ = 2*(B2+Z)*(BS2*W2 + (W2-1)*Z3) / (W2*Z3);
    const double df0dW = 2*(B2+Z)*(B2+Z) / W3;
    const double df1dZ = 1.0 + BS2/Z3 - c2p->gamf/W2;
    const double df1dW = B2/W3 + (2*Z - D*W)/W3 * c2p->gamf;

    double f[2];
    double J[4], G[4];


    // Evaluation of the function, and its Jacobian
    // -------------------------------------------------------------------------
    f[0] = -S2  + (Z+B2)*(Z+B2)*(W2-1)/W2 - (2*Z+B2)*BS2/Z2;      // eqn (84)
    f[1] = -Tau +  Z+B2 - Pre - 0.5*B2/W2 -      0.5*BS2/Z2 - D;  // eqn (85)

    J[0] = df0dZ; J[1] = df0dW;
    J[2] = df1dZ; J[3] = df1dW;

    // G in the inverse Jacobian
    // -------------------------------------------------------------------------
    const double det = J[0]*J[3] - J[1]*J[2];
    G[0] =  J[3]/det; G[1] = -J[1]/det;
    G[2] = -J[2]/det; G[3] =  J[0]/det;      // G = J^{-1}

    const double dZ = -(G[0]*f[0] + G[1]*f[1]); // Matrix multiply, dx = -G . f
    const double dW = -(G[2]*f[0] + G[3]*f[1]);

    // Bracketing the root
    // -------------------------------------------------------------------------
    double Z_new = Z + dZ;
    double W_new = W + dW;

    Z_new = (Z_new > c2p->smlZ) ? Z_new : -Z_new;
    Z_new = (Z_new < c2p->bigZ) ? Z_new :  Z;

    W_new = (W_new > c2p->smlW) ? W_new : c2p->smlW;
    W_new = (W_new < c2p->bigW) ? W_new : c2p->bigW;

    Z = Z_new;
    W = W_new;

    // -------------------------------------------------------------------------
    error = fabs(dZ/Z) + fabs(dW/W);
    ++c2p->Iterations;

    if (c2p->Iterations == c2p->MaxIterations) {
      return SRMHD_C2P_MAXITER;
    }
  }

  return srmhd_c2p_reconstruct_prim(c2p, Z, W, P);
}


// Solution based on Noble et. al. (2006), using Z = rho h W^2 as the single
// unkown. Unfortunately, Noble uses 'W' for what I call Z. This function should
// really be called '1dz', but I use this name to reflect the name of the
// section in which it appears.
// -----------------------------------------------------------------------------
int srmhd_c2p_solve_noble1dw(srmhd_c2p *c2p, double *P)
{
  int bad_input = srmhd_c2p_check_cons(c2p, c2p->Cons);
  if (bad_input) {
    return bad_input;
  }
  // Starting values
  // ---------------------------------------------------------------------------
  c2p->Iterations = 0;
  const double D   = c2p->D;
  const double B2  = c2p->B2;
  const double S2  = c2p->S2;
  const double BS2 = c2p->BS2;
  const double Tau = c2p->Tau;

  double error = 1.0;
  double Z = c2p->Z_start;

  double f, g;

  while (error > c2p->Tolerance) {

    const double Z2  = Z*Z;
    const double Z3  = Z*Z2;
    const double a   = S2*Z2 + BS2*(B2 + 2*Z);
    const double b   = (B2 + Z)*(B2 + Z)*Z2;
    const double ap  = 2*(S2*Z + BS2);          // da/dZ
    const double bp  = 2*Z*(B2 + Z)*(B2 + 2*Z); // db/dZ
    const double V2  = a / b;
    const double W2  = 1.0 / (1.0 - V2);
    const double W   = sqrt(W2);
    const double W3  = W*W2;
    const double Pre = (D/W) * (Z/(D*W) - 1.0) * c2p->gamf;

    const double dv2dZ    = (ap*b - bp*a) / (b*b); // (a'b - b'a) / b^2
    const double delPdelZ = c2p->gamf/W2;
    const double delPdelW = c2p->gamf * (D/W2 - 2*Z/W3);
    const double dWdv2    = 0.5*W3;
    const double dPdZ     = delPdelW * dWdv2 * dv2dZ + delPdelZ;

    f = Tau + D - 0.5*B2*(1+V2) + 0.5*BS2/Z2 - Z + Pre; // equation (29)
    g = -0.5*B2*dv2dZ - BS2/Z3 - 1.0 + dPdZ;

    const double dZ = -f/g;

    // -------------------------------------------------------------------------
    double Z_new = Z + dZ;

    Z_new = (Z_new > c2p->smlZ) ? Z_new : -Z_new;
    Z_new = (Z_new < c2p->bigZ) ? Z_new :  Z;

    Z = Z_new;

    error = fabs(dZ/Z);
    ++c2p->Iterations;

    if (c2p->Iterations == c2p->MaxIterations) {
      return SRMHD_C2P_MAXITER;
    }
  }

  // Recover the W value from the converged Z value.
  // -------------------------------------------------------------------------
  const double Z2  = Z*Z;
  const double a   = S2*Z2 + BS2*(B2 + 2*Z);
  const double b   = (B2 + Z)*(B2 + Z)*Z2;
  const double V2  = a / b;
  const double W2  = 1.0 / (1.0 - V2);
  const double W   = sqrt(W2);

  return srmhd_c2p_reconstruct_prim(c2p, Z, W, P);
}

const char *srmhd_c2p_get_error(srmhd_c2p *c2p, int error)
{
  switch (error) {
  case SRMHD_C2P_SUCCESS:
    return "successful inversion from conserved to primitive";
  case SRMHD_C2P_CONS_CONTAINS_NAN:
    return "input conserved state contained nan's";
  case SRMHD_C2P_CONS_NEGATIVE_DENSITY:
    return "input conserved state has negative density";
  case SRMHD_C2P_CONS_NEGATIVE_ENERGY:
    return "input conserved state has negative energy";
  case SRMHD_C2P_PRIM_CONTAINS_NAN:
    return "derived primitive state contains nan's";
  case SRMHD_C2P_PRIM_NEGATIVE_PRESSURE:
    return "derived primitive state has negative pressure";
  case SRMHD_C2P_PRIM_NEGATIVE_RESTMASS:
    return "derived primitive state contains negative density";
  case SRMHD_C2P_PRIM_SUPERLUMINAL:
    return "derived primitive state contains superluminal velocity";
  case SRMHD_C2P_MAXITER:
    return "rootfinder hit maximum iterations";
  default:
    return "unkown error";
  }
}
