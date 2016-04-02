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
