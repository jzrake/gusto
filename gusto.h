#ifndef GUSTO_HEADER
#define GUSTO_HEADER
#include "ser.h"



/*
 * =====================================================================
 * Data structures
 * =====================================================================
 */
struct aux_variables;
struct gusto_sim;
struct mesh_face;
struct mesh_cell;


typedef void (*OpBoundaryCon)(struct gusto_sim *sim);
typedef void (*OpInitialMesh)(struct gusto_user *user);
typedef const char **(*OpInitialData)(struct gusto_user *user,
				      struct aux_variables *A,
				      double *X);

struct gusto_sim
{
  int num_rows;
  double smallest_cell_length; /* for CFL condition */
  struct mesh_face *faces;
  struct mesh_row *rows;
  struct gusto_user user;
  struct gusto_status status;
  OpBoundaryCon boundary_con;
  OpInitialMesh initial_mesh;
  OpInitialData initial_data;
} ;


struct aux_variables
{
  double velocity_four_vector[4];
  double magnetic_four_vector[4];
  double momentum_density[4];
  double comoving_mass_density;
  double gas_pressure;
  double magnetic_pressure;
  double enthalpy_density; /* includes magnetic part: h + b^2 */
  double vector_potential; /* vector potential */
  double R;                /* cylindrical radius */
} ;


struct aux_geometry
{
  double cylindrical_radius;
  double volume_element;  /* sqrt(full metric determinant) */
  double area_element[4]; /* sqrt(reduced metric determinant) */
  double line_element[4]; /* scale factors */
} ;


struct mesh_face
{
  double x[4];
  double x_rk[4];
  double Fhat[8]; /* Godunov fluxes */
  double dA[4];
  struct aux_variables aux;
  struct aux_geometry geom;
  struct mesh_cell *cells[2];
  struct mesh_face *next;
  struct mesh_face *prev;
} ;


struct mesh_cell
{
  double x[4];    /* cell centroid */
  double U[8];    /* conserved masses */
  double U_rk[8]; /* cached values for RK stepping */
  double dA[4];   /* conversion factors for conserved masses and densities */
  char cell_type;
  struct aux_variables aux;
  struct aux_geometry geom;
  struct mesh_face *faces[2];
  struct mesh_cell *next;
  struct mesh_cell *prev;
} ;


struct mesh_row
{
  struct mesh_face *faces;
  struct mesh_cell *cells;
  struct mesh_cell *cells_end;
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
void gusto_mesh_generate(struct gusto_sim *sim);
void gusto_mesh_compute_geometry(struct gusto_sim *sim);
void gusto_mesh_advance_vertices(struct gusto_sim *sim, double dt);
void gusto_cache_rk(struct gusto_sim *sim);
void gusto_average_rk(struct gusto_sim *sim, double b);


/* Operations on one or more variable states */
void gusto_riemann(struct aux_variables *AL, struct aux_variables *AR,
		   double nhat[4], double Fhat[8], double s);
void gusto_to_conserved(struct aux_variables *A, double U[8], double dA[4]);
void gusto_default_aux(struct aux_variables *A);
void gusto_complete_aux(struct aux_variables *A);
void gusto_geometric_source_terms(struct aux_variables *A, double Udot[8]);
void gusto_electric_field(struct aux_variables *A, double E[4]);
int gusto_from_conserved(struct aux_variables *A, double U[8], double dA[4]);
int gusto_wavespeeds(struct aux_variables *A, double n[4], double evals[8]);
int gusto_fluxes(struct aux_variables *A, double n[4], double F[8]);


void gusto_geometry(struct aux_geometry *G, double x[4]);


/* Operations on the whole simulation */
void gusto_recover_variables(struct gusto_sim *sim);
void gusto_initial_data(struct gusto_sim *sim);
void gusto_compute_face_velocities(struct gusto_sim *sim);
void gusto_compute_variables_at_vertices(struct gusto_sim *sim);
void gusto_compute_fluxes(struct gusto_sim *sim);
void gusto_compute_face_magnetic_flux(struct gusto_sim *sim);
void gusto_compute_cell_field_from_faces(struct gusto_sim *sim);
void gusto_compute_cell_magnetic_field(struct gusto_sim *sim);
void gusto_transmit_fluxes(struct gusto_sim *sim, double dt);
void gusto_smooth_electric_field(struct gusto_sim *sim);
void gusto_advance_vector_potential(struct gusto_sim *sim, double dt);
void gusto_add_source_terms(struct gusto_sim *sim, double dt);


void gusto_init(struct gusto_sim *sim);
void gusto_free(struct gusto_sim *sim);
void gusto_config_from_user(struct gusto_sim *sim);
void gusto_write_checkpoint(struct gusto_sim *sim, const char *fname);
int gusto_is_valid(struct gusto_sim *sim);



OpBoundaryCon gusto_lookup_boundary_con(const char *user_key);
OpInitialMesh gusto_lookup_initial_mesh(const char *user_key);
OpInitialData gusto_lookup_initial_data(const char *user_key);


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
#define VEC4_AVG(x,y) {0,0.5*(x[1]+y[1]),0.5*(x[2]+y[2]),0.5*(x[3]+y[3])}
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

#define VEC4_SUB2(x,y,z) z[0]=0;z[1]=x[1]-y[1];z[2]=x[2]-y[2];z[3]=x[3]-y[3];


#endif /* GUSTO_HEADER */
