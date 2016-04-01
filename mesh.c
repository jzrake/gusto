#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utlist.h"
#include "gusto.h"



void gusto_mesh_report(struct gusto_sim *sim)
{
  printf("[gusto] number rows  ... %d\n", sim->num_rows);
  printf("[gusto] number verts ... %d\n", gusto_mesh_count(sim, 'v', -1));
  printf("[gusto] number cells ... %d\n", gusto_mesh_count(sim, 'c', -1));
  printf("[gusto] number faces ... %d\n", gusto_mesh_count(sim, 'f', 0));
}



void gusto_mesh_clear(struct gusto_sim *sim)
{
  struct mesh_vert *V0, *V1;
  struct mesh_cell *C0, *C1;
  struct mesh_face *F0, *F1;

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

  DL_FOREACH_SAFE(sim->faces, F0, F1) {
    DL_DELETE(sim->faces, F0);
    free(F0);
  }

  free(sim->rows);
  sim->num_rows = 0;
}



int gusto_mesh_count(struct gusto_sim *sim, char which, int n)
{
  struct mesh_vert *V;
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
    case 'v':
      DL_COUNT(sim->rows[n].verts, V, count);
      break;
    case 'c':
      DL_COUNT(sim->rows[n].cells, C, count);
      break;
    case 'f':
      DL_COUNT(sim->faces, F, count);
      break;
    }
  }

  return count;
}



void gusto_mesh_generate_verts(struct gusto_sim *sim)
{
  sim->num_rows = sim->user.N[1];
  sim->rows = (struct mesh_row *) malloc(sim->num_rows*sizeof(struct mesh_row));

  int row_size = sim->user.N[0];
  int ngR = 0; /* number of ghost cells in the R direction */
  int ngz = 0; /* number of ghost cells in the z direction */
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
    for (int i=0; i<row_size+ngz+1; ++i) {

      int num_interior_R = sim->num_rows - 2 * ngR;
      int num_interior_z = row_size      - 2 * ngz;
      int n_real = n - ngR;
      int i_real = i - ngz;
      double dR = (R1 - R0) / num_interior_R;
      double dz = (z1 - z0) / num_interior_z;

      VL = (struct mesh_vert *) malloc(sizeof(struct mesh_vert));
      VR = (struct mesh_vert *) malloc(sizeof(struct mesh_vert));

      VL->row_index = n;
      VR->row_index = n;
      VL->col_index = 2*i + 0;
      VR->col_index = 2*i + 1;

      VL->x[0] = 0.0;
      VR->x[0] = 0.0;
      VL->x[1] = R0 + (n_real + 0) * dR;
      VR->x[1] = R0 + (n_real + 1) * dR;
      VL->x[2] = 0.0;
      VR->x[2] = 0.0;
      VL->x[3] = z0 + i_real * dz + 0.0 * dz * (n % 2 == 0);
      VR->x[3] = z0 + i_real * dz + 0.0 * dz * (n % 2 == 0);

      VL->x[1] += 0.0 * ((double) rand() / RAND_MAX - 0.5);
      VL->x[3] += 0.0 * ((double) rand() / RAND_MAX - 0.5);

      for (int d=0; d<4; ++d) { /* vertex velocities */
	VL->v[d] = 0.0;
	VR->v[d] = 0.0;
      }

      DL_APPEND(sim->rows[n].verts, VL);
      DL_APPEND(sim->rows[n].verts, VR);
    }
  }
}



void gusto_mesh_generate_cells(struct gusto_sim *sim)
{
  struct mesh_vert *VL, *VR;
  struct mesh_cell *C;

  /* ----------------------------------------------------------------------
   * Add the cells, one for each i index (z coordinate) except the last.
   * ---------------------------------------------------------------------- */
  for (int n=0; n<sim->num_rows; ++n) {

    VL = sim->rows[n].verts;
    VR = sim->rows[n].verts->next;

    while (VR->next) {
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
  struct mesh_vert *VL, *VR;
  struct mesh_cell *CL, *CR;
  struct mesh_face *F;

  sim->faces = NULL;

  for (int n=0; n<sim->num_rows; ++n) {

    VL = sim->rows[n].verts;
    VR = sim->rows[n].verts->next;
    CL = NULL;
    CR = sim->rows[n].cells;

    /* ----------------------------------------------------------------------
     * Add the faces, one for each z coordinate (pair of vertices).
     * ---------------------------------------------------------------------- */
    while (1) {

      F = (struct mesh_face *) malloc(sizeof(struct mesh_face));
      F->verts[0] = VL;
      F->verts[1] = VR;
      F->cells[0] = CL;
      F->cells[1] = CR;

      DL_APPEND(sim->faces, F);

      if (CR) { /* it was not the last row */
	VL = VL->next->next;
	VR = VR->next->next;
	CL = CR;
	CR = CR->next;
      }
      else {
	break;
      }
    }
  }


  for (int n=0; n<sim->num_rows-1; ++n) {

    struct mesh_vert *V0, *V1;
    struct mesh_vert *Vm, *Vp;
    struct mesh_cell *Cm, *Cp;
    struct mesh_face *F;

    Vm = sim->rows[n+0].verts->next;
    Vp = sim->rows[n+1].verts;
    Cm = NULL;
    Cp = NULL;

    if (Vm->x[3] < Vp->x[3]) {
      V0 = Vm;
      Vm = Vm->next->next;
      Cm = sim->rows[n+0].cells;
    }
    else {
      V0 = Vp;
      Vp = Vp->next->next;
      Cp = sim->rows[n+1].cells;
    }

    char which;
    int proceed = 1;

    /*
     * -------------------------------------------------------------------------
     * We start with a face F that is attached at its lower end to the vertex
     * V0. There are two possible vertices to which F might be attached at the
     * other end - Vm and Vp, which are on the left and right of the bounding
     * surface respectively. The correct one is whichever has the smaller z
     * coordinate. If that is Vp, then we attach the other end of F to V1 = Vp,
     * and advance Vp to the next vertex along its of the surface. We then add
     * that face to the list, reset V0 = V1, and begin a new face whose lower
     * end point is the new V0. If either Vp or Vm is NULL, then we are out of
     * vertices on that side. In that case, the construction is finished when
     * the next vertex on the remaining side is also NULL.
     * -------------------------------------------------------------------------
     */

    while (proceed) {

      if (Vm && Vp) {
	if (Vm->x[3] < Vp->x[3]) {
	  which = 'm';
	  V1 = Vm;
	  Vm = Vm->next ? Vm->next->next : NULL;
	}
	else {
	  which = 'p';
	  V1 = Vp;
	  Vp = Vp->next ? Vp->next->next : NULL;
	}
      }
      else if (Vm) {
	which = 'm';
	V1 = Vm;
	Vm = Vm->next ? Vm->next->next : NULL;
	proceed = Vm != NULL;
      }
      else if (Vp) {
	which = 'p';
	V1 = Vp;
	Vp = Vp->next ? Vp->next->next : NULL;
	proceed = Vp != NULL;
      }

      F = (struct mesh_face *) malloc(sizeof(struct mesh_face));
      DL_APPEND(sim->faces, F);
      F->verts[0] = V0;
      F->verts[1] = V1;
      F->cells[0] = Cm;
      F->cells[1] = Cp;
      V0 = V1;

      /* Now we advance the cell on the side of the vertex that advanced. The
	 first face always borders only one cell, thus on the first one or more
	 iterations, either Cm or Cp is NULL. Similarly, on the last one or more
	 iterations either Cm or Cp may be NULL.
	 */

      if (which == 'm') {
	if (Vm) { /* m side is still active */
	  Cm = Cm ? Cm->next : sim->rows[n+0].cells;
	}
	else {
	  Cm = NULL;
	}
      }

      if (which == 'p') {
	if (Vp) { /* p side is still active */
	  Cp = Cp ? Cp->next : sim->rows[n+1].cells;
	}
	else {
	  Cp = NULL;
	}
      }
    }
  }
}



void gusto_mesh_compute_geometry(struct gusto_sim *sim)
{
  struct mesh_cell *C;
  struct mesh_face *F;

  sim->smallest_cell_length = 1.0;

  /*
   * Set the cell geometry
   * ---------------------------------------------------------------------------
   */
  for (int n=0; n<sim->num_rows; ++n) {
    DL_FOREACH(sim->rows[n].cells, C) {

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
      C->x[0] = 0.0;
      C->x[1] = 0.25 * (C->verts[0]->x[1] + C->verts[1]->x[1] +
			C->verts[2]->x[1] + C->verts[3]->x[1]);
      C->x[2] = 0.25 * (C->verts[0]->x[2] + C->verts[1]->x[2] +
			C->verts[2]->x[2] + C->verts[3]->x[2]);
      C->x[3] = 0.25 * (C->verts[0]->x[3] + C->verts[1]->x[3] +
			C->verts[2]->x[3] + C->verts[3]->x[3]);

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

  /*
   * Set the face geometry
   * ---------------------------------------------------------------------------
   */
  DL_FOREACH(sim->faces, F) {
    double df[4] = {0, 0, 1, 0};
    double dl[4] = VEC4_SUB(F->verts[1]->x, F->verts[0]->x);
    double dA[4] = VEC4_CROSS(dl, df);
    F->nhat[0] = VEC4_MOD(dA);
    F->nhat[1] = dA[1] / F->nhat[0];
    F->nhat[2] = dA[2] / F->nhat[0];
    F->nhat[3] = dA[3] / F->nhat[0];
  }
}



void gusto_mesh_advance_vertices(struct gusto_sim *sim, double dt)
{
  struct mesh_vert *V;
  for (int n=0; n<sim->num_rows; ++n) {
    DL_FOREACH(sim->rows[n].verts, V) {
      V->x[1] += V->v[1] * dt;
      V->x[2] += V->v[2] * dt;
      V->x[3] += V->v[3] * dt;
    }
  }
}
