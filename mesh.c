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

    double R0 = sim->user.domain[0];
    double R1 = sim->user.domain[1];
    double dlogX = log(R1 / R0) / row_size;

    double X = 1.0;
    double R = R0;
    double z = 0.0;
    double dX, dR, dz, dB, Bp;

    for (int i=0; i<row_size+1; ++i) {

      dX = X * dlogX; /* change in distance along the field line */
      Bp = gusto_geometry_step_along_field(sim, R, z, &dR, &dz, &dB, dX);

      F = (struct mesh_face *) malloc(sizeof(struct mesh_face));
      F->y[1] = R;
      F->y[3] = z;
      F->x[1] = 1.0;
      F->x[3] = X;

      F->geom.dXB = dB / dX;
      F->geom.dXR = dR / dX;

      F->cells[0] = NULL;
      F->cells[1] = NULL;
      gusto_default_aux(&F->aux);
      DL_APPEND(sim->rows[n].faces, F);

      X += dX;
      R += dR;
      z += dz;
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

      if      (C->prev->next == NULL) C->cell_type = 'g';
      else if (C->next       == NULL) C->cell_type = 'g';
      else                            C->cell_type = 'i';

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
      gusto_geometry(sim, &F->geom, F->y);

      F->dA[0] = F->geom.area_element[3];
      F->dA[1] = F->geom.area_element[1];
      F->dA[2] = F->geom.area_element[2];
      F->dA[3] = F->geom.area_element[3];
    }

    for (C=sim->rows[n].cells; C; C=C->next) {

      double dX = C->faces[1]->x[3] - C->faces[0]->x[3];

      /* Set the cell centroid to the average of the face's coordinate. */
      C->y[1] = 0.5 * (C->faces[0]->y[1] + C->faces[1]->y[1]);
      C->y[3] = 0.5 * (C->faces[0]->y[3] + C->faces[1]->y[3]);
      C->x[1] = 0.5 * (C->faces[0]->x[1] + C->faces[1]->x[1]);
      C->x[3] = 0.5 * (C->faces[0]->x[3] + C->faces[1]->x[3]);

      /* Compute the scale factors for that position. */
      gusto_geometry(sim, &C->geom, C->y);

      C->geom.dXB = (C->faces[1]->geom.poloidal_field -
		     C->faces[0]->geom.poloidal_field) / dX;
      C->geom.dXR = (C->faces[1]->geom.cylindrical_radius -
		     C->faces[0]->geom.cylindrical_radius) / dX;

      /* These are the cell's volume element (dA[0]), and area elements along
	 each axis. */
      C->dA[0] = C->geom.volume_element  * dX;
      C->dA[1] = C->geom.area_element[1] * dX;
      C->dA[2] = C->geom.area_element[2] * dX;
      C->dA[3] = C->geom.area_element[3];

      if (dX < sim->smallest_cell_length) sim->smallest_cell_length = dX;
    }

  }
}



void gusto_mesh_advance_vertices(struct gusto_sim *sim, double dt)
{

}
