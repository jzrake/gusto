#include <math.h>
#include <stdio.h>
#include "gusto.h"


#define PARABOLOIDAL 1
#define MONOPOLE 0

static double poloidal_field(double R, double z, double *BR, double *Bz);
static double flux_function(double R, double z);


int gusto_validate_geometry(struct gusto_sim *sim)
{
  struct mesh_cell *C;

  for (int n=0; n<sim->num_rows; ++n) {
    for (C=sim->rows[n].cells; C; C=C->next) {

      double Y = C->x[1];
      double R = C->y[1];

      double dXB_tab = C->geom.dXB;
      double dXR_tab = C->geom.dXR;

      double dXB_tru = 0;
      double dXR_tru = 0;

      if (PARABOLOIDAL) {
	dXB_tru = (pow(1 + pow(Y/R, +2), -2.0) - 1) * 2 / (R * R);
	dXR_tru =  pow(1 + pow(Y/R, -2), -0.5);
      }
      if (MONOPOLE) {
	dXB_tru = -2 * pow(R, -3);
	dXR_tru = 1.0;
      }

      printf("dB/dX = [%e %e], dR/dX = [%e %e]\n",
	     dXB_tab, dXB_tru, dXR_tab, dXR_tru);

    }
  }

  return 0;
}



void gusto_geometry(struct aux_geometry *G, double y[4])
{
  double BR, Bz;
  double R = y[1];
  double z = y[3];
  double B = poloidal_field(R, z, &BR, &Bz);
  double h[4] = { -1, 1 / (B * R), R, 1 };
  G->cylindrical_radius = R;
  G->poloidal_field = B;
  G->volume_element = h[1] * h[2] * h[3];
  G->area_element[1] = h[2] * h[3];
  G->area_element[2] = h[3] * h[1];
  G->area_element[3] = h[1] * h[2];
  G->line_element[1] = h[1];
  G->line_element[2] = h[2];
  G->line_element[3] = h[3];
}



void gusto_geometry_source_terms(struct aux_variables *A,
				 struct aux_geometry *G, double Udot[8])
{
  double R = G->cylindrical_radius;
  double B = G->poloidal_field;
  double *u = A->velocity_four_vector;
  double *b = A->magnetic_four_vector;
  double H0 = A->enthalpy_density;
  double p0 = A->gas_pressure + A->magnetic_pressure;

  /* Udot[S33] = (2*p0 + u[2]*u[2]*H0 - b[2]*b[2]) / A->R; */

  Udot[S33] = -p0 * G->dXB / B + (u[2]*u[2]*H0 - b[2]*b[2]) * G->dXR / R;
}



double flux_function(double R, double z)
{
  if (PARABOLOIDAL) {
    return sqrt(R*R + z*z) - z;
  }
  if (MONOPOLE) {
    return 1 - z / sqrt(R*R + z*z);
  }
}



double poloidal_field(double R, double z, double *BR, double *Bz)
{
  double dR = R * 1e-6;
  double dz = z * 1e-6;

  double YRm2 = flux_function(R-2*dR, z);
  double YRm1 = flux_function(R-1*dR, z);
  double YRp1 = flux_function(R+1*dR, z);
  double YRp2 = flux_function(R+2*dR, z);
  double Yzm2 = flux_function(R, z-2*dz);
  double Yzm1 = flux_function(R, z-1*dz);
  double Yzp1 = flux_function(R, z+1*dz);
  double Yzp2 = flux_function(R, z+2*dz);

  double dRY = (-YRp2 + 8 * YRp1 - 8 * YRm1 + YRm2) / (12 * dR);
  double dzY = (-Yzp2 + 8 * Yzp1 - 8 * Yzm1 + Yzm2) / (12 * dz);

  *BR = -dzY / R;
  *Bz = +dRY / R;

  return sqrt((*BR) * (*BR) + (*Bz) * (*Bz));
}



double gusto_geometry_step_along_field(struct gusto_sim *sim,
				     double R, double z,
				     double *dR, double *dz, double *dB,
				     double dchi)
/*
 * Take an RK4 step along the field line. Return the poloidal field strenth at
 * the beginning of the step.
 *
 */
{
  double BR, Bz, Bp;
  double dR1, dz1;
  double dR2, dz2;
  double dR3, dz3;
  double dR4, dz4;

  Bp = poloidal_field(R, z, &BR, &Bz);
  dR1 = dchi * BR / Bp;
  dz1 = dchi * Bz / Bp;
  Bp = poloidal_field(R + 0.5 * dR1, z + 0.5 * dz1, &BR, &Bz);
  dR2 = dchi * BR / Bp;
  dz2 = dchi * Bz / Bp;
  Bp = poloidal_field(R + 0.5 * dR2, z + 0.5 * dz2, &BR, &Bz);
  dR3 = dchi * BR / Bp;
  dz3 = dchi * Bz / Bp;
  Bp = poloidal_field(R + 1.0 * dR3, z + 1.0 * dz3, &BR, &Bz);
  dR4 = dchi * BR / Bp;
  dz4 = dchi * Bz / Bp;

  *dR = (dR1 + 2 * dR2 + 2 * dR3 + dR4) / 6;
  *dz = (dz1 + 2 * dz2 + 2 * dz3 + dz4) / 6;

  /* Also compute the change of poloidal field magnitude for the given step dchi
     along the field line. */
  double Bp0 = poloidal_field(R, z, &BR, &Bz);
  double Bp1 = poloidal_field(R + *dR, z + *dz, &BR, &Bz);

  *dB = Bp1 - Bp0;

  return Bp0;
}
