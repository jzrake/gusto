#include <stdlib.h>
#include <math.h>
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



#define WRITE_VARIABLE_WAVESPEED(member, alias) do {			\
    struct mesh_cell *C = sim->rows[n].cells;				\
    double u, evals[8];							\
    double Xhat[4] = {0, 0, 0, 1};					\
    for (int i=0; i<data_size; ++i) {					\
      gusto_wavespeeds(&C->aux, Xhat, evals);				\
      u = evals[member] / sqrt(1 - pow(evals[member], 2));		\
      data[i] = u;							\
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
    hid_t cg = H5Gcreate(row, "cells", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Write the cell data
     * -------------------------------------------------------------------------
     */
    if (1) {
      hsize_t data_size = gusto_mesh_count(sim, 'c', n);
      hid_t spc = H5Screate_simple(1, &data_size, NULL);
      double *data = (double *) malloc(data_size * sizeof(double));

      WRITE_VARIABLE_CELL(y[1], "R");
      WRITE_VARIABLE_CELL(y[3], "z");
      WRITE_VARIABLE_CELL(x[3], "X");
      WRITE_VARIABLE_CELL(aux.velocity_four_vector[0], "u0");
      WRITE_VARIABLE_CELL(aux.velocity_four_vector[1], "u1");
      WRITE_VARIABLE_CELL(aux.velocity_four_vector[2], "u2");
      WRITE_VARIABLE_CELL(aux.velocity_four_vector[3], "u3");
      WRITE_VARIABLE_CELL(aux.magnetic_four_vector[0], "b0");
      WRITE_VARIABLE_CELL(aux.magnetic_four_vector[1], "b1");
      WRITE_VARIABLE_CELL(aux.magnetic_four_vector[2], "b2");
      WRITE_VARIABLE_CELL(aux.magnetic_four_vector[3], "b3");
      WRITE_VARIABLE_CELL(aux.comoving_mass_density, "dg");
      WRITE_VARIABLE_CELL(aux.gas_pressure, "pg");
      WRITE_VARIABLE_WAVESPEED(0, "-fast");
      WRITE_VARIABLE_WAVESPEED(1, "-alfv");
      WRITE_VARIABLE_WAVESPEED(2, "-slow");
      WRITE_VARIABLE_WAVESPEED(3, "fluid");
      WRITE_VARIABLE_WAVESPEED(5, "+slow");
      WRITE_VARIABLE_WAVESPEED(6, "+alfv");
      WRITE_VARIABLE_WAVESPEED(7, "+fast");

      free(data);
      H5Sclose(spc);
    }

    H5Gclose(cg);
    H5Gclose(row);
  }

  H5Gclose(grp);
  H5Fclose(h5f);

  gusto_read_write_status(&sim->status, fname, 'a');
  gusto_read_write_user(&sim->user, fname, 'a');
}
