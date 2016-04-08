#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "gusto.h"



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



double gusto_quad_area_centroid(double x0[4],
				double x1[4], double x2[4],
				double x3[4], double x4[4])
{
  /*
   *  2-------4
   *  |       |
   *  |   0   |
   *  |       |
   *  1-------3
   */

  x0[0] = 0;
  x0[1] = 0.25 * (x1[1] + x2[1] + x3[1] + x4[1]);
  x0[2] = 0.25 * (x1[2] + x2[2] + x3[2] + x4[2]);
  x0[3] = 0.25 * (x1[3] + x2[3] + x3[3] + x4[3]);

  double dx1[4] = VEC4_SUB(x1, x0);
  double dx2[4] = VEC4_SUB(x2, x0);
  double dx3[4] = VEC4_SUB(x3, x0);
  double dx4[4] = VEC4_SUB(x4, x0);
  double dA1[4] = VEC4_CROSS(dx1, dx2);
  double dA2[4] = VEC4_CROSS(dx2, dx4);
  double dA3[4] = VEC4_CROSS(dx4, dx3);
  double dA4[4] = VEC4_CROSS(dx3, dx1);

  return 0.5 * (VEC4_MOD(dA1) + VEC4_MOD(dA2) +
		VEC4_MOD(dA3) + VEC4_MOD(dA4));
}



void gusto_curlA1(struct mesh_cell *C, int i0, int i1, int i2,
		  double crlA[2])
/*
 * Evaluate the curl of A using the three vertices i0, i1, i2, where the vertex
 * i0 is treated as the origin and used twice.
 */
{
  double A0 = C->verts[i0]->aux[0].vector_potential;
  double A1 = C->verts[i1]->aux[0].vector_potential;
  double A2 = C->verts[i2]->aux[0].vector_potential;

  double R0 = C->verts[i0]->aux[0].R;
  double R1 = C->verts[i1]->aux[0].R;
  double R2 = C->verts[i2]->aux[0].R;

  double x0[2] = {C->verts[i0]->x[1], C->verts[i0]->x[3]};
  double x1[2] = {C->verts[i1]->x[1], C->verts[i1]->x[3]};
  double x2[2] = {C->verts[i2]->x[1], C->verts[i2]->x[3]};

  double n[2][2] = { {x1[0] - x0[0], x1[1] - x0[1]},
		     {x2[0] - x0[0], x2[1] - x0[1]} };
  double det = n[0][0] * n[1][1] - n[0][1] * n[1][0];
  double m[2][2] = { {+n[1][1]/det, -n[0][1]/det},
		     {-n[1][0]/det, +n[0][0]/det} }; /* inverse of n */

  double delY[2] = {A1*R1 - A0*R0, A2*R2 - A0*R0};
  double grdY[2] = {m[0][0] * delY[0] + m[0][1] * delY[1],
		    m[1][0] * delY[0] + m[1][1] * delY[1]}; /* (dR A, dz A) */

  crlA[0] = -grdY[1] / R0;
  crlA[1] = +grdY[0] / R0;
}



void gusto_curlA2(struct mesh_cell *C, double crlA[2])
/*
 * Evaluate the curl of A using the x-type stencil across the cell.
 */
{
  double A0 = C->verts[0]->aux[0].vector_potential;
  double A1 = C->verts[1]->aux[0].vector_potential;
  double A2 = C->verts[2]->aux[0].vector_potential;
  double A3 = C->verts[3]->aux[0].vector_potential;

  double R0 = C->verts[0]->aux[0].R;
  double R1 = C->verts[1]->aux[0].R;
  double R2 = C->verts[2]->aux[0].R;
  double R3 = C->verts[3]->aux[0].R;

  double x0[2] = {C->verts[0]->x[1], C->verts[0]->x[3]};
  double x1[2] = {C->verts[1]->x[1], C->verts[1]->x[3]};
  double x2[2] = {C->verts[2]->x[1], C->verts[2]->x[3]};
  double x3[2] = {C->verts[3]->x[1], C->verts[3]->x[3]};

  double n[2][2] = { {x1[0] - x2[0], x1[1] - x2[1]},
		     {x3[0] - x0[0], x3[1] - x0[1]} };

  double det = n[0][0] * n[1][1] - n[0][1] * n[1][0];
  double m[2][2] = { {+n[1][1]/det, -n[0][1]/det},
		     {-n[1][0]/det, +n[0][0]/det} }; /* inverse of n */

  double delY[2] = {A1*R1 - A2*R2, A3*R3 - A0*R0};
  double grdY[2] = {m[0][0] * delY[0] + m[0][1] * delY[1],
		    m[1][0] * delY[0] + m[1][1] * delY[1]}; /* (dR Y, dz Y) */

  crlA[0] = -grdY[1] / C->aux[0].R;
  crlA[1] = +grdY[0] / C->aux[0].R;
}
