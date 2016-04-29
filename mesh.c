#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "utlist.h"
#include "gusto.h"



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
      if (n < sim->num_rows) {
	DL_COUNT(sim->rows[n].cells, C, count);
      }
      break;
    case 'f':
      if (n < sim->num_rows) {
	DL_COUNT(sim->rows[n].faces, F, count);
      }
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
  int row_size = sim->user.N;
  struct mesh_face *F;
  struct mesh_cell *C;

  gusto_mesh_clear(sim, 'a');
  sim->num_rows = 1;
  sim->rows = (struct mesh_row *) malloc(sim->num_rows*sizeof(struct mesh_row));

  for (int n=0; n<sim->num_rows; ++n) {

    sim->rows[n].faces = NULL;
    sim->rows[n].cells = NULL;

    double x0 = 1.0;
    double x1 = 2.0;
    double dx = (x1 - x0) / row_size;

    for (int i=0; i<row_size+1; ++i) {
      F = (struct mesh_face *) malloc(sizeof(struct mesh_face));
      F->x[1] = x0 + i * dx;
      gusto_default_aux(&F->aux);
      DL_APPEND(sim->rows[n].faces, F);
    }

    for (F=sim->rows[n].faces; F->next; F=F->next) {
      C = (struct mesh_cell *) malloc(sizeof(struct mesh_cell));
      C->faces[0] = F;
      C->faces[1] = F->next;
      gusto_default_aux(&C->aux);
      DL_APPEND(sim->rows[n].cells, C);
    }

    for (C=sim->rows[n].cells; C; C=C->next) {
      C->faces[0]->cells[1] = C;
      C->faces[1]->cells[0] = C;
    }
  }
}



void gusto_mesh_compute_geometry(struct gusto_sim *sim)
{
  struct mesh_face *F;
  struct mesh_cell *C;

  sim->smallest_cell_length = 1.0;

  for (int n=0; n<sim->num_rows; ++n) {

    for (F=sim->rows[n].faces; F; F=F->next) {
      gusto_geometry(&F->geom, F->x);

      F->dA[0] = F->geom.area_element[3];
      F->dA[1] = F->geom.line_element[1];
      F->dA[2] = F->geom.line_element[2];
      F->dA[3] = F->geom.line_element[3];
    }

    for (C=sim->rows[n].cells; C; C=C->next) {

      /* Set the cell centroid to the average of the face's coordinate. */
      C->x[1] = 0.5 * (C->faces[0]->x[1] + C->faces[1]->x[1]);

      /* Compute the scale factors for that position. */
      gusto_geometry(&C->geom, C->x);
      C->aux.R = C->geom.cylindrical_radius;

      double dx = C->faces[1]->x[1] - C->faces[0]->x[1];

      /* These are the cell's volume element (dA[0]), and area elements along
	 each axis. Only the */
      C->dA[0] = C->geom.volume_element  * dx;
      C->dA[1] = 1.0;
      C->dA[2] = C->geom.area_element[2] * dx;
      C->dA[3] = C->geom.area_element[3];

      if (dx < sim->smallest_cell_length) sim->smallest_cell_length = dx;
    }
  }
}



void gusto_mesh_advance_vertices(struct gusto_sim *sim, double dt)
{

}
