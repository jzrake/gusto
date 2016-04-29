#include "gusto.h"



void gusto_geometry(struct aux_geometry *G, double x[4])
{
  double r = x[3];
  double B = 1.0 / (r * r);
  double h[4] = { -1, 1 / (B * r), r, 1 };
  G->cylindrical_radius = r;
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
  double *u = A->velocity_four_vector;
  double *b = A->magnetic_four_vector;
  double H0 = A->enthalpy_density;
  double p0 = A->gas_pressure + A->magnetic_pressure;

  Udot[S33] = (2*p0 + u[2]*u[2]*H0 - b[2]*b[2]) / A->R;

  /* Udot[S33] = -p0 * dBdx / B + (u[2]*u[2]*H0 - b[2]*b[2]) * dRdx / R; */
}
