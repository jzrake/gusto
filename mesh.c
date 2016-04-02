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



void gusto_mesh_clear(struct gusto_sim *sim, char which)
{
  struct mesh_vert *V0, *V1;
  struct mesh_cell *C0, *C1;
  struct mesh_face *F0, *F1;

  for (int n=0; n<sim->num_rows; ++n) {
    if (which == 'v' || which == 'a') {
      DL_FOREACH_SAFE(sim->rows[n].verts, V0, V1) {
	DL_DELETE(sim->rows[n].verts, V0);
	free(V0);
      }
    }
    if (which == 'c' || which == 'a') {
      DL_FOREACH_SAFE(sim->rows[n].cells, C0, C1) {
	DL_DELETE(sim->rows[n].cells, C0);
	free(C0);
      }
    }
  }

  if (which == 'f' || which == 'a') {
    DL_FOREACH_SAFE(sim->faces, F0, F1) {
      DL_DELETE(sim->faces, F0);
      free(F0);
    }
  }

  if (which == 'v' || which == 'a') {
    free(sim->rows);
    sim->num_rows = 0;
  }
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



void vert_pos_staggered_cartesian(struct gusto_sim *sim, struct mesh_vert *V,
				  int n, int i, int NR, int Nz, char pm)
{
  double R0 = sim->user.domain[0];
  double R1 = sim->user.domain[1];
  double z0 = sim->user.domain[2];
  double z1 = sim->user.domain[3];
  double dR = (R1 - R0) / NR;
  double dz = (z1 - z0) / Nz;

  V->x[0] = 0.0;
  V->x[1] = R0 + (n + (pm == 'p')) * dR;
  V->x[2] = 0.0;
  V->x[3] = z0 + i * dz + 0.5 * dz * (n % 2 == 0);
}



void vert_pos_radial(struct gusto_sim *sim, struct mesh_vert *V,
		     int n, int i, int Nt, int Nr, char pm)
{
  double r0 = sim->user.domain[0];
  double r1 = sim->user.domain[1];
  double t0 = sim->user.domain[2];
  double t1 = sim->user.domain[3];
  double dr = log(r1 / r0) / Nr;
  double dt = (t1 - t0) / Nt;

  double r = r0 * exp(i * dr);
  double t = t0 + dt * (n + (pm == 'p'));

  /* Theta is measured from the equatorial plane, not the pole. */
  V->x[0] = 0.0;
  V->x[1] = r * cos(t);
  V->x[2] = 0.0;
  V->x[3] = r * sin(t);
}



void gusto_mesh_generate_verts(struct gusto_sim *sim)
{
  sim->num_rows = sim->user.N[1];
  sim->rows = (struct mesh_row *) malloc(sim->num_rows*sizeof(struct mesh_row));

  int row_size = sim->user.N[0];
  int ngR = sim->user.ng[0]; /* number of ghost cells in the R direction */
  int ngz = sim->user.ng[1]; /* number of ghost cells in the z direction */

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

      int NR = sim->num_rows - 2 * ngR; /* number of interior R cells */
      int Nz = row_size      - 2 * ngz; /* number of interior z cells */
      int n_real = n - ngR;
      int i_real = i - ngz;

      VL = (struct mesh_vert *) malloc(sizeof(struct mesh_vert));
      VR = (struct mesh_vert *) malloc(sizeof(struct mesh_vert));

      VL->row_index = n;
      VR->row_index = n;
      VL->col_index = 2*i + 0;
      VR->col_index = 2*i + 1;

      vert_pos_radial(sim, VL, n_real, i_real, NR, Nz, 'm');
      vert_pos_radial(sim, VR, n_real, i_real, NR, Nz, 'p');

      /* vert_pos_staggered_cartesian(sim, VL, n_real, i_real, NR, Nz, 'm'); */
      /* vert_pos_staggered_cartesian(sim, VR, n_real, i_real, NR, Nz, 'p'); */

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

  gusto_mesh_clear(sim, 'f');

  for (int n=0; n<sim->num_rows; ++n) {

    VL = sim->rows[n].verts;
    VR = sim->rows[n].verts->next;
    CL = NULL;
    CR = sim->rows[n].cells;

    /* ----------------------------------------------------------------------
     * Add longitudinal faces, one for each pair of vertices along a row.
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


  /* ----------------------------------------------------------------------
   * Add lateral faces.
   * ---------------------------------------------------------------------- */
  for (int n=0; n<sim->num_rows-1; ++n) {

    struct mesh_vert *V0, *V1;
    struct mesh_vert *Vm, *Vp;
    struct mesh_cell *Cm, *Cp;
    struct mesh_face *F;

    Vm = sim->rows[n+0].verts->next;
    Vp = sim->rows[n+1].verts;
    Cm = NULL;
    Cp = NULL;

    /* vector pointing along the local longitudinal axis */
    double dl[4] = VEC4_SUB(Vp->next->next->x, Vp->x);
    double lm = VEC4_DOT(Vm->x, dl);
    double lp = VEC4_DOT(Vp->x, dl);

    if (lm < lp) {
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
     * surface respectively. The correct one is whichever has the smaller 'l'
     * coordinate. If that is Vp, then we attach the other end of F to V1 = Vp,
     * and advance Vp to the next vertex along its of the surface. We then add
     * that face to the list, reset V0 = V1, and begin a new face whose lower
     * end point is the new V0. If either Vp or Vm is NULL, then we are out of
     * vertices on that side. In that case, the construction is finished when
     * the next vertex on the remaining side is also NULL.
     *
     * The 'l' coordinate is the measured along the local longitudinal
     * axis. That would be the vertex's z coordinate if the rows are along z.
     *
     -------------------------------------------------------------------------
     */

    while (proceed) {

      if (Vm && Vp) {

	double lm = VEC4_DOT(Vm->x, dl);
	double lp = VEC4_DOT(Vp->x, dl);

	if (lm < lp) {
	  //printf("m\n");
	  which = 'm';
	  V1 = Vm;
	  Vm = Vm->next ? Vm->next->next : NULL;
	}
	else {
	  //printf("p\n");
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
      F->verts[0] = V0;
      F->verts[1] = V1;
      F->cells[0] = Cm;
      F->cells[1] = Cp;

      VEC4_SUB2(V1->x, V0->x, dl);

      if (VEC4_MOD(dl) > 1e-12) { /* So we don't create faces with zero area */
	DL_APPEND(sim->faces, F);
      }
      else {
	free(F);
      }

      V0 = V1;

      /* Now we advance the cell on the side of the vertex that advanced. The
	 first and last faces always border only one cell, so on the first and
	 last iteration(s), either Cm or Cp is NULL.
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
      double dAR0[4] = VEC4_CROSS(dphi, dz0);
      double dAR1[4] = VEC4_CROSS(dphi, dz1);
      double dAz0[4] = VEC4_CROSS(dR0, dphi);
      double dAz1[4] = VEC4_CROSS(dR1, dphi);
      double area = gusto_quad_area_centroid(C->x, /* also sets C->x */
					     C->verts[0]->x, C->verts[1]->x,
					     C->verts[2]->x, C->verts[3]->x);

      C->dA[0] = area;
      C->dA[1] = 0.50 * (VEC4_MOD(dAR0) + VEC4_MOD(dAR1));
      C->dA[2] = area;
      C->dA[3] = 0.50 * (VEC4_MOD(dAz0) + VEC4_MOD(dAz1));

      if (sim->user.coordinates == 'p') {
	C->dA[0] *= C->x[1];
	C->dA[1] *= C->x[1];
	C->dA[3] *= C->x[1];
      }

      /* Cell's longitudinal axis */
      /* C->zhat[0] = 0.0; */
      /* C->zhat[1] = 0.5 * (dAz0[1] + dAz1[1]); */
      /* C->zhat[2] = 0.5 * (dAz0[2] + dAz1[2]); */
      /* C->zhat[3] = 0.5 * (dAz0[3] + dAz1[3]); */
      /* VEC4_NORMALIZE(C->zhat); */


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
    double dx[4] = VEC4_SUB(F->verts[1]->x, F->verts[0]->x);
    double x0[4] = VEC4_AVG(F->verts[1]->x, F->verts[0]->x);
    double dA[4] = VEC4_CROSS(df, dx);

    F->nhat[0] = VEC4_MOD(dA);
    F->nhat[1] = dA[1] / F->nhat[0];
    F->nhat[2] = dA[2] / F->nhat[0];
    F->nhat[3] = dA[3] / F->nhat[0];



    /* If in cylindrical coordinates then the face area is multiplied by
       R. There's a bug here: why does dividing by 2 make it work? */
    if (sim->user.coordinates == 'p') {
      F->nhat[0] *= 0.5 * x0[1];
    }



    /* This ensures that the face's normal vector points from cells[0] ->
       cells[1]. */
    double dz[4] = {0,0,0,0};
    if (F->cells[0] && F->cells[1]) {
      VEC4_SUB2(F->cells[1]->x, F->cells[0]->x, dz);
    }
    else if (F->cells[0]) {
      VEC4_SUB2(x0, F->cells[0]->x, dz);
    }
    else if (F->cells[1]) {
      VEC4_SUB2(F->cells[1]->x, x0, dz);
    }
    if (VEC4_DOT(dz, F->nhat) < 0.0) {
      F->nhat[1] *= -1;
      F->nhat[2] *= -1;
      F->nhat[3] *= -1;
    }


    /* This is a self-consistency check. It should be guarenteed by the previous
       lines. */
    if (F->cells[0] && F->cells[1]) {
      double dL[4] = VEC4_SUB(F->cells[1]->x, F->cells[0]->x);
      double aligned = VEC4_DOT(dL, F->nhat);
      if (aligned < 0) {
	printf("[gusto] ERROR: a face is on backwards\n");
      }
    }
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
