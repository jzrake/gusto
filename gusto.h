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
struct mesh_vert;
struct mesh_face;
struct mesh_cell;


typedef void (*OpBoundaryCon)(struct gusto_sim *sim);
typedef void (*OpInitialMesh)(struct gusto_user *user, struct mesh_vert *V);
typedef const char **(*OpInitialData)(struct gusto_user *user,
				      struct aux_variables *A, double *X);


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
  double R;                /* factor for cylindrical coordinates */
} ;

struct mesh_vert
{
  double x[4];
  double v[4];
  double Efield;
  int row_index;
  int col_index;
  int num_counts;
  struct aux_variables aux[4];
  struct mesh_cell *cell;
  struct mesh_vert *next;
  struct mesh_vert *prev;
} ;

struct mesh_face
{
  double Fhat[8]; /* Godunov fluxes */
  double nhat[4]; /* 0: face area, 1,2,3: unit normal */
  double Bflux;   /* magnetic flux through face */
  double length;  /* distance between face's 2 vertices */
  struct mesh_vert *verts[2];
  struct mesh_cell *cells[2];
  struct mesh_face *next;
  struct mesh_face *prev;
} ;

struct mesh_cell
{
  char cell_type;
  double x[4];
  double dA[4];                /* 0: volume, dA1, dA2, dA3: cross-sections */
  double zhat[4];              /* 1,2,3: cell's longitudinal axis (unit) */
  double U[8];                 /* total mass, energy, momentum, magnetic flux */
  double weightA;
  double weightB;
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
double gusto_quad_area_centroid(double x0[4],
				double x1[4], double x2[4],
				double x3[4], double x4[4]);
void gusto_curlA1(struct mesh_cell *C, int i0, int i1, int i2,
		  double crlA[2]);
void gusto_curlA2(struct mesh_cell *C, double crlA[2]);

/* Mesh operations */
int gusto_mesh_count(struct gusto_sim *sim, char which, int n);
void gusto_mesh_report(struct gusto_sim *sim);
void gusto_mesh_clear(struct gusto_sim *sim, char which);
void gusto_mesh_generate_verts(struct gusto_sim *sim);
void gusto_mesh_generate_cells(struct gusto_sim *sim);
void gusto_mesh_generate_faces(struct gusto_sim *sim);
void gusto_mesh_compute_geometry(struct gusto_sim *sim);
void gusto_mesh_advance_vertices(struct gusto_sim *sim, double dt);



/* Operations on one or more variable states */
void gusto_riemann(struct aux_variables *AL, struct aux_variables *AR,
		   double nhat[4], double Fhat[8], double s);
void gusto_to_conserved(struct aux_variables *A, double U[8], double dA[4]);
void gusto_default_aux(struct aux_variables *A);
void gusto_complete_aux(struct aux_variables *A);
void gusto_cylindrical_source_terms(struct aux_variables *A, double Udot[8]);
void gusto_spherical_source_terms(struct aux_variables *A, double Udot[8]);
int gusto_from_conserved(struct aux_variables *A, double U[8], double dA[4]);
int gusto_wavespeeds(struct aux_variables *A, double n[4], double evals[8]);
int gusto_fluxes(struct aux_variables *A, double n[4], double F[8]);



/* Operations on the whole simulation */
void gusto_recover_variables(struct gusto_sim *sim);
void gusto_initial_data(struct gusto_sim *sim);
void gusto_compute_vertex_velocities(struct gusto_sim *sim);
void gusto_compute_variables_at_vertices(struct gusto_sim *sim);
void gusto_compute_fluxes(struct gusto_sim *sim);
void gusto_compute_cell_magnetic_field(struct gusto_sim *sim);
void gusto_transmit_fluxes(struct gusto_sim *sim, double dt);
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
