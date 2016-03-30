#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utlist.h"
#include "gusto.h"


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



int gusto_mesh_count(struct gusto_sim *sim, char which, int n)
{
  struct mesh_vert *V;
  struct mesh_cell *C;
  int count = 0;

  if (which != 'r' && n == -1) {
    for (int m=0; m<sim->num_rows; ++m) {
      count += gusto_mesh_count(sim, which, m);
    }
  }

  else {
    switch (which) {
    case 'v':
      DL_COUNT(sim->rows[n].verts, V, count);
      break;
    case 'c':
      DL_COUNT(sim->rows[n].cells, C, count);
      break;
    case 'r':
      count = sim->num_rows;
      break;
    }
  }

  return count;
}


void gusto_mesh_clear(struct gusto_sim *sim)
{
  struct mesh_vert *V0, *V1;
  struct mesh_cell *C0, *C1;

  for (int n=0; n<sim->num_rows; ++n) {
    DL_FOREACH_SAFE(sim->rows[n].verts, V0, V1) {
      DL_DELETE(sim->rows[n].verts, V0);
      free(V0);
    }
    DL_FOREACH_SAFE(sim->rows[n].cells, C0, C1) {
      DL_DELETE(sim->rows[n].cells, C0);
      free(C0);
    }
  }

  free(sim->rows);
  sim->num_rows = 0;
}



void gusto_mesh_generate_verts(struct gusto_sim *sim)
{
  sim->num_rows = sim->user.N[0];
  sim->rows = (struct mesh_row *) malloc(sim->num_rows*sizeof(struct mesh_row));

  int row_size = sim->user.N[1];
  int ngR = 1; /* number of ghost cells in the R direction */
  int ngz = 1; /* number of ghost cells in the z direction */
  double R0 = 0.0;
  double R1 = 1.0;
  double z0 = 0.0;
  double z1 = 1.0;

  for (int n=0; n<sim->num_rows; ++n) {

    /* Iteration variables for going over the linked list: */
    struct mesh_vert *VL, *VR;

    /*
     * Initiate the lists by setting them to NULL. Each row is lined with
     * vertices to its left and right. The list of vertices for a given row
     * alternates from left to right, so for a vertex VL on the left of the row,
     * its neighbor at the same (or similar) z coordinate is VR = VL->next. The
     * verts data member for a given cell is arranged like so:
     *
     *  1-------3
     *  |       |
     *  |       |
     *  0-------2
     *
     *  z
     *  ^
     *  |
     *  |
     *  x----> R
     *
     */
    sim->rows[n].verts = NULL;
    sim->rows[n].cells = NULL;

    /* ----------------------------------------------------------------------
     * Add the vertices, two for each i index (z coordinate).
     * ---------------------------------------------------------------------- */
    for (int i=0; i<row_size+1; ++i) {

      int num_interior_R = sim->num_rows - 2 * ngR;
      int num_interior_z = row_size      - 2 * ngz;
      int n_real = n - ngR;
      int i_real = i - ngz;
      double dR = (R1 - R0) / num_interior_R;
      double dz = (z1 - z0) / num_interior_z;

      VL = (struct mesh_vert *) malloc(sizeof(struct mesh_vert));
      VR = (struct mesh_vert *) malloc(sizeof(struct mesh_vert));

      VL->x[0] = 0.0;
      VR->x[0] = 0.0;

      VL->x[1] = R0 + (n_real + 0) * dR;
      VR->x[1] = R0 + (n_real + 1) * dR;
      VL->x[2] = 0.0;
      VR->x[2] = 0.0;
      VL->x[3] = z0 + i_real * dz;
      VR->x[3] = z0 + i_real * dz;

      for (int d=0; d<4; ++d) { /* vertex velocities */
	VL->v[d] = 0.0;
	VR->v[d] = 0.0;
      }

      DL_APPEND(sim->rows[n].verts, VL);
      DL_APPEND(sim->rows[n].verts, VR);
    }

    VL = sim->rows[n].verts;
    VR = sim->rows[n].verts->next;
  }

  printf("num_rows: %d\n", gusto_mesh_count(sim, 'r', -1));
  printf("num_cells[tot]: %d\n", gusto_mesh_count(sim, 'c', -1));
  printf("num_cells[row]: %d\n", gusto_mesh_count(sim, 'c', 0));
  printf("num_verts[tot]: %d\n", gusto_mesh_count(sim, 'v', -1));
  printf("num_verts[row]: %d\n", gusto_mesh_count(sim, 'v', 0));
}



void gusto_mesh_generate_cells(struct gusto_sim *sim)
{
  struct mesh_vert *VL, *VR;
  struct mesh_cell *C;
  /* ----------------------------------------------------------------------
   * Add the cells, one for each i index (z coordinate) except the last.
   * ---------------------------------------------------------------------- */
  for (int n=0; n<sim->num_rows; ++n) {
    int row_size = gusto_mesh_count(sim, 'v', n) / 2;
    for (int i=0; i<row_size; ++i) {
      C = (struct mesh_cell *) malloc(sizeof(struct mesh_cell));
      C->verts[0] = VL; VL = VL->next->next;
      C->verts[1] = VL;
      C->verts[2] = VR; VR = VR->next->next;
      C->verts[3] = VR;
      DL_APPEND(sim->rows[n].cells, C);
    }
  }
}



void gusto_mesh_generate_faces(struct gusto_sim *sim)
{

}



void gusto_mesh_compute_geometry(struct gusto_sim *sim)
{
  struct mesh_cell *C;

  DL_FOREACH(sim->cells, C) {

    /*
     * These vectors define the 2-forms on the cell. There are four, one for
     * each corner.
     */
    double dR0[4] = VEC4_SUB(C->verts[2]->x, C->verts[0]->x);
    double dR1[4] = VEC4_SUB(C->verts[3]->x, C->verts[1]->x);
    double dz0[4] = VEC4_SUB(C->verts[1]->x, C->verts[0]->x);
    double dz1[4] = VEC4_SUB(C->verts[3]->x, C->verts[2]->x);
    double dphi[4] = {0, 0, 1, 0};
    double dAf0[4] = VEC4_CROSS(dz0, dR0);
    double dAf1[4] = VEC4_CROSS(dz0, dR1);
    double dAf2[4] = VEC4_CROSS(dz1, dR0);
    double dAf3[4] = VEC4_CROSS(dz1, dR1);
    double dAR0[4] = VEC4_CROSS(dphi, dz0);
    double dAR1[4] = VEC4_CROSS(dphi, dz1);
    double dAz0[4] = VEC4_CROSS(dR0, dphi);
    double dAz1[4] = VEC4_CROSS(dR1, dphi);

    /* Cell's centroid position */
    C->x[1] = 0.25 * (C->verts[0]->x[0] + C->verts[1]->x[0] +
		      C->verts[2]->x[0] + C->verts[3]->x[0]);
    C->x[2] = 0.25 * (C->verts[0]->x[1] + C->verts[1]->x[1] +
		      C->verts[2]->x[1] + C->verts[3]->x[1]);
    C->x[3] = 0.25 * (C->verts[0]->x[2] + C->verts[1]->x[2] +
		      C->verts[2]->x[2] + C->verts[3]->x[2]);

    /* Cell's area and volume forms */
    C->dA[1] = 0.50 * (VEC4_MOD(dAR0) + VEC4_MOD(dAR1));
    C->dA[2] = 0.25 * (VEC4_MOD(dAf0) + VEC4_MOD(dAf1) +
		       VEC4_MOD(dAf2) + VEC4_MOD(dAf3));
    C->dA[3] = 0.50 * (VEC4_MOD(dAz0) + VEC4_MOD(dAz1));
    C->dA[0] = C->dA[2]; /* Volume and phi cross section are the same */

    /* Cell's longitudinal axis */
    C->zhat[0] = 0.0;
    C->zhat[1] = 0.5 * (dAz0[1] + dAz1[1]);
    C->zhat[2] = 0.5 * (dAz0[2] + dAz1[2]);
    C->zhat[3] = 0.5 * (dAz0[3] + dAz1[3]);

    VEC4_NORMALIZE(C->zhat);

    double dAR = gusto_min3(VEC4_MOD(dAR0), VEC4_MOD(dAR1), 1.0);
    double dAz = gusto_min3(VEC4_MOD(dAz0), VEC4_MOD(dAz1), 1.0);

    if (dAR < sim->smallest_cell_length) sim->smallest_cell_length = dAR;
    if (dAz < sim->smallest_cell_length) sim->smallest_cell_length = dAz;
  }
}
