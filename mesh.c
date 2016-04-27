#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "utlist.h"
#include "gusto.h"

#define STERADIAN (2 * M_PI)


void gusto_mesh_report(struct gusto_sim *sim)
{
  printf("[gusto] number rows  ... %d\n", sim->num_rows);
  printf("[gusto] number cells ... %d\n", gusto_mesh_count(sim, 'c', -1));
  printf("[gusto] number faces ... %d\n", gusto_mesh_count(sim, 'f', 0));
}



void gusto_mesh_clear(struct gusto_sim *sim, char which)
{
  struct mesh_cell *C0, *C1;
  struct mesh_face *F0, *F1;

  for (int n=0; n<sim->num_rows; ++n) {
    if (which == 'f' || which == 'a') {
      DL_FOREACH_SAFE(sim->rows[n].faces, F0, F1) {
	DL_DELETE(sim->rows[n].faces, F0);
	free(F0);
      }
    }
    if (which == 'c' || which == 'a') {
      DL_FOREACH_SAFE(sim->rows[n].cells, C0, C1) {
	DL_DELETE(sim->rows[n].cells, C0);
	free(C0);
      }
    }
  }

  if (which == 'v' || which == 'a') {
    free(sim->rows);
    sim->num_rows = 0;
  }
}



int gusto_mesh_count(struct gusto_sim *sim, char which, int n)
{
  struct mesh_cell *C;
  struct mesh_face *F;
  int count = 0;

  if (n == -1) {
    for (int m=0; m<sim->num_rows; ++m) {
      count += gusto_mesh_count(sim, which, m);
    }
  }

  else {
    switch (which) {
    case 'c':
      DL_COUNT(sim->rows[n].cells, C, count);
      break;
    case 'f':
      DL_COUNT(sim->rows[n].faces, F, count);
      break;
    }
  }

  return count;
}



OpInitialMesh gusto_lookup_initial_mesh(const char *user_key)
{
  const char *keys[] = {
    NULL
  } ;
  OpInitialMesh vals[] = {
    NULL } ;
  int n = 0;
  const char *key;
  OpInitialMesh val;
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
  printf("[gusto] ERROR: no such initial_mesh=%s\n", user_key);
  return NULL;
}



void gusto_mesh_generate(struct gusto_sim *sim)
{

}



void gusto_mesh_compute_geometry(struct gusto_sim *sim)
{

}



void gusto_mesh_advance_vertices(struct gusto_sim *sim, double dt)
{

}
