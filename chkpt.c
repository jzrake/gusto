#include <stdlib.h>
#include <hdf5.h>
#include "utlist.h"
#include "gusto.h"



#define WRITE_VARIABLE_VERT(member, alias) do {				\
    struct mesh_vert *V = sim->rows[n].verts;				\
    for (int i=0; i<data_size; ++i) {					\
      data[i] = V->member;						\
      V = V->next;							\
    }									\
    hid_t set = H5Dcreate(vg, alias, H5T_NATIVE_DOUBLE, spc,		\
			  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);	\
    H5Dwrite(set, H5T_NATIVE_DOUBLE, spc, spc, H5P_DEFAULT, data);	\
    H5Dclose(set);							\
  } while(0)								\



#define WRITE_VARIABLE_CELL(member, alias) do {				\
    struct mesh_cell *C = sim->rows[n].cells;				\
    for (int i=0; i<data_size; ++i) {					\
      data[i] = C->member;						\
      C = C->next;							\
    }									\
    hid_t set = H5Dcreate(cg, alias, H5T_NATIVE_DOUBLE, spc,		\
			  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);	\
    H5Dwrite(set, H5T_NATIVE_DOUBLE, spc, spc, H5P_DEFAULT, data);	\
    H5Dclose(set);							\
  } while(0)								\



void gusto_write_checkpoint(struct gusto_sim *sim, const char *fname)
{
  char default_name[64];
  char row_name[64];

  if (fname == NULL) {
    snprintf(default_name, 64, "chkpt.%04d.h5", sim->status.checkpoint_number);
    fname = default_name;
  }
  printf("[gusto] writing checkpoint %s\n", fname);

  hid_t h5f = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hid_t grp = H5Gcreate(h5f, "rows", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  for (int n=0; n<sim->num_rows; ++n) {

    snprintf(row_name, 64, "row_%06d", n);
    hid_t row = H5Gcreate(grp, row_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t vg = H5Gcreate(row, "verts", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t cg = H5Gcreate(row, "cells", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Write the vertex data
     * -------------------------------------------------------------------------
     */
    if (1) {
      hsize_t data_size = gusto_mesh_count(sim, 'v', n);
      hid_t spc = H5Screate_simple(1, &data_size, NULL);
      double *data = (double *) malloc(data_size * sizeof(double));

      WRITE_VARIABLE_VERT(x[1], "x1");
      WRITE_VARIABLE_VERT(x[2], "x2");
      WRITE_VARIABLE_VERT(x[3], "x3");
      free(data);
      H5Sclose(spc);
    }

    /*
     * Write the cell data
     * -------------------------------------------------------------------------
     */
    if (1) {
      hsize_t data_size = gusto_mesh_count(sim, 'c', n);
      hid_t spc = H5Screate_simple(1, &data_size, NULL);
      double *data = (double *) malloc(data_size * sizeof(double));

      WRITE_VARIABLE_CELL(x[1], "x1");
      WRITE_VARIABLE_CELL(x[2], "x2");
      WRITE_VARIABLE_CELL(x[3], "x3");
      WRITE_VARIABLE_CELL(aux[0].velocity_four_vector[1], "u1");
      WRITE_VARIABLE_CELL(aux[0].velocity_four_vector[2], "u2");
      WRITE_VARIABLE_CELL(aux[0].velocity_four_vector[3], "u3");
      WRITE_VARIABLE_CELL(aux[0].magnetic_four_vector[1], "b1");
      WRITE_VARIABLE_CELL(aux[0].magnetic_four_vector[2], "b2");
      WRITE_VARIABLE_CELL(aux[0].magnetic_four_vector[3], "b3");
      WRITE_VARIABLE_CELL(aux[0].comoving_mass_density, "dg");
      WRITE_VARIABLE_CELL(aux[0].gas_pressure, "pg");

      WRITE_VARIABLE_CELL(verts[0]->x[1], "v0.x1");
      WRITE_VARIABLE_CELL(verts[0]->x[2], "v0.x2");
      WRITE_VARIABLE_CELL(verts[0]->x[3], "v0.x3");
      WRITE_VARIABLE_CELL(verts[1]->x[2], "v1.x2");
      WRITE_VARIABLE_CELL(verts[1]->x[3], "v1.x3");
      WRITE_VARIABLE_CELL(verts[1]->x[1], "v1.x1");
      WRITE_VARIABLE_CELL(verts[2]->x[3], "v2.x3");
      WRITE_VARIABLE_CELL(verts[2]->x[1], "v2.x1");
      WRITE_VARIABLE_CELL(verts[2]->x[2], "v2.x2");
      WRITE_VARIABLE_CELL(verts[3]->x[1], "v3.x1");
      WRITE_VARIABLE_CELL(verts[3]->x[2], "v3.x2");
      WRITE_VARIABLE_CELL(verts[3]->x[3], "v3.x3");

      free(data);
      H5Sclose(spc);
    }

    H5Gclose(vg);
    H5Gclose(cg);
    H5Gclose(row);
  }

  /*
   * Write the face data
   * -------------------------------------------------------------------------
   */
  if (1) {

    /*
     * Face data has shape [Nfaces, 2, 2]
     *
     *       [face, vert0/vert1, row_index/col_index]
     * -----------------------------------------------------------------------
     */
    hsize_t data_size[3] = { gusto_mesh_count(sim, 'f', 0), 2, 2 };
    hid_t spc = H5Screate_simple(3, data_size, NULL);
    int *data = (int *) malloc(data_size[0] * 4 * sizeof(int));

    struct mesh_face *F;

    int data_index = 0;

    DL_FOREACH(sim->faces, F) {
      data[data_index++] = F->verts[0]->row_index;
      data[data_index++] = F->verts[0]->col_index;
      data[data_index++] = F->verts[1]->row_index;
      data[data_index++] = F->verts[1]->col_index;
    }

    hid_t set = H5Dcreate(h5f, "faces", H5T_NATIVE_INT, spc,
			  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(set, H5T_NATIVE_INT, spc, spc, H5P_DEFAULT, data);
    H5Dclose(set);
    H5Sclose(spc);
    free(data);
  }

  H5Gclose(grp);
  H5Fclose(h5f);

  gusto_read_write_status(&sim->status, fname, 'a');
  gusto_read_write_user(&sim->user, fname, 'a');
}
