#ifndef GUSTO_HEADER
#define GUSTO_HEADER


#include "ser.h"
#define gamma_law_index (4./3)

/* indices into conserved variable arrays */
enum { DDD, TAU, S11, S22, S33, B11, B22, B33 };



/*
 * Data structures
 * =====================================================================
 */

struct aux_variables
{
  double velocity_four_vector[4];
  double magnetic_four_vector[4];
  double momentum_density[4]; /* T^{0,:} */
  double comoving_mass_density;
  double gas_pressure;
  double magnetic_pressure;
} ;

struct mesh_vert
{
  double x[4];
  double v[4];
  struct aux_variables aux[4];
  struct mesh_cell *cell;
  int num; /* temporary solution for averaging cells to verts */
  int row_index;
  int col_index;
  struct mesh_vert *next;
  struct mesh_vert *prev;
} ;

struct mesh_face
{
  double Fhat[8]; /* Godunov fluxes */
  double nhat[4]; /* 0: face area, 1,2,3: unit normal */
  struct mesh_vert *verts[2];
  struct mesh_cell *cells[2];
  struct mesh_face *next;
  struct mesh_face *prev;
} ;

struct mesh_cell
{
  double x[4];
  double dA[4];                /* 0: volume, dA1, dA2, dA3: cross-sections */
  double zhat[4];              /* 1,2,3: cell's longitudinal axis (unit) */
  double U[8];                 /* total mass, energy, momentum, magnetic flux */
  struct mesh_vert *verts[4];
  struct aux_variables aux[5]; /* aux vars at different cell locations */
  struct mesh_cell *next;
  struct mesh_cell *prev;
} ;

struct mesh_row
{
  struct mesh_vert *verts;
  struct mesh_cell *cells;
} ;

struct gusto_sim
{
  struct gusto_user user;
  struct gusto_status status;

  int num_cells_max;
  int num_cells;
  int num_faces_max;
  int num_faces;

  int num_rows;
  int *row_size;

  double smallest_cell_length; /* for CFL condition */

  struct mesh_vert **verts;
  struct mesh_face *faces;
  struct mesh_cell *cells;
  struct mesh_row *rows;
} ;




/*
 * Function prototypes
 * =====================================================================
 */
void gusto_write_checkpoint(struct gusto_sim *sim, const char *fname);
void *gusto_start_clock();
double gusto_stop_clock(void *clock_s);


int gusto_mesh_count(struct gusto_sim *sim, char which, int n);
void gusto_mesh_clear(struct gusto_sim *sim);
void gusto_mesh_generate_verts(struct gusto_sim *sim);
void gusto_mesh_generate_cells(struct gusto_sim *sim);
void gusto_mesh_generate_faces(struct gusto_sim *sim);
void gusto_mesh_compute_geometry(struct gusto_sim *sim);


#define gusto_max3(a,b,c)(((a)>(b))?(((a)>(c))?(a):(c)):(((b)>(c))?(b):(c)))
#define gusto_min3(a,b,c)(((a)<(b))?(((a)<(c))?(a):(c)):(((b)<(c))?(b):(c)))


#endif /* GUSTO_HEADER */
