#ifndef GUSTO_HEADER
#define GUSTO_HEADER
#include "ser.h"



/*
 * =====================================================================
 * Data structures
 * =====================================================================
 */

struct gusto_sim
{
  int num_rows;
  double smallest_cell_length; /* for CFL condition */
  struct mesh_face *faces;
  struct mesh_row *rows;
  struct gusto_user user;
  struct gusto_status status;
} ;

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
  int row_index;
  int col_index;
  struct aux_variables aux[4];
  struct mesh_cell *cell;
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



/*
 * =====================================================================
 * Function prototypes
 * =====================================================================
 */
/* General utilities */
void *gusto_start_clock();
double gusto_stop_clock(void *clock_s);


/* Mesh operations */
int gusto_mesh_count(struct gusto_sim *sim, char which, int n);
void gusto_mesh_report(struct gusto_sim *sim);
void gusto_mesh_clear(struct gusto_sim *sim, char which);
void gusto_mesh_generate_verts(struct gusto_sim *sim);
void gusto_mesh_generate_cells(struct gusto_sim *sim);
void gusto_mesh_generate_faces(struct gusto_sim *sim);
void gusto_mesh_compute_geometry(struct gusto_sim *sim);
void gusto_mesh_advance_vertices(struct gusto_sim *sim, double dt);


/* Initial data functions */
void initial_data_cylindrical_shock(struct aux_variables *A, double *X);
void initial_data_sound_wave(struct aux_variables *A, double *X);


/* Operations on one or more variable states */
void gusto_riemann(struct aux_variables *AL, struct aux_variables *AR,
		   double nhat[4], double Fhat[8], double s);
void gusto_vars_to_conserved(struct aux_variables *A, double U[8], double dA[4]);
void gusto_vars_complete_aux(struct aux_variables *A);
int gusto_vars_from_conserved(struct aux_variables *A, double U[8], double dA[4]);
int gusto_wavespeeds(struct aux_variables *A, double n[4], double evals[8]);


/* Operations on the whole simulation */
void gusto_recover_variables(struct gusto_sim *sim);
void gusto_initial_data(struct gusto_sim *sim);
void gusto_compute_vertex_velocities(struct gusto_sim *sim);
void gusto_enforce_boundary_condition(struct gusto_sim *sim);
void gusto_compute_variables_at_vertices(struct gusto_sim *sim);
void gusto_compute_fluxes(struct gusto_sim *sim);
void gusto_transmit_fluxes(struct gusto_sim *sim, double dt);


void gusto_free(struct gusto_sim *sim);
void gusto_write_checkpoint(struct gusto_sim *sim, const char *fname);


/*
 * =====================================================================
 * Macro definitions
 * =====================================================================
 */

#define DDD 0 /* indices into conserved variable arrays */
#define TAU 1
#define S11 2
#define S22 3
#define S33 4
#define B11 5
#define B22 6
#define B33 7
#define gamma_law_index (4./3)


#define gusto_max3(a,b,c)(((a)>(b))?(((a)>(c))?(a):(c)):(((b)>(c))?(b):(c)))
#define gusto_min3(a,b,c)(((a)<(b))?(((a)<(c))?(a):(c)):(((b)<(c))?(b):(c)))


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



#endif /* GUSTO_HEADER */
