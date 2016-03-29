#include <stdlib.h>
#include <hdf5.h>
#include "gusto.h"


#define WRITE_VARIABLE(member, alias) do {				\
    for (int i=0; i<sim->row_size[n]; ++i) {				\
      data[i] = sim->verts[n][i].member;				\
    }									\
    hid_t set = H5Dcreate(row, alias, H5T_NATIVE_DOUBLE, spc,		\
			  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);	\
    H5Dwrite(set, H5T_NATIVE_DOUBLE, spc, spc, H5P_DEFAULT, data);	\
    H5Dclose(set);							\
  } while(0)								\


void gusto_write_checkpoint(struct gusto_sim *sim, const char *fname)
{  
  char rowname[64];
  hid_t h5f = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hid_t grp = H5Gcreate(h5f, "rows", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  for (int n=0; n<sim->num_rows; ++n) {
    snprintf(rowname, 64, "row_%06d", n);
    double *data = (double *) malloc(sim->row_size[n] * sizeof(double));    
    hsize_t row_size = sim->row_size[n];
    hid_t spc = H5Screate_simple(1, &row_size, NULL);
    hid_t row = H5Gcreate(grp, rowname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    WRITE_VARIABLE(x[1], "x1");
    WRITE_VARIABLE(x[2], "x2");
    WRITE_VARIABLE(x[3], "x3");
    WRITE_VARIABLE(aux[0].velocity_four_vector[1], "u1");
    WRITE_VARIABLE(aux[0].velocity_four_vector[2], "u2");
    WRITE_VARIABLE(aux[0].velocity_four_vector[3], "u3");
    WRITE_VARIABLE(aux[0].magnetic_four_vector[1], "b1");
    WRITE_VARIABLE(aux[0].magnetic_four_vector[2], "b2");
    WRITE_VARIABLE(aux[0].magnetic_four_vector[3], "b3");
    WRITE_VARIABLE(aux[0].comoving_mass_density, "dg");
    WRITE_VARIABLE(aux[0].gas_pressure, "pg");

    H5Gclose(row);
    H5Sclose(spc);
    free(data);
  }
  
  H5Gclose(grp);
  H5Fclose(h5f);
}
