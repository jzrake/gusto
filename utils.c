#include <stdlib.h>
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
