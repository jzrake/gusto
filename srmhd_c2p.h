#ifndef SRMHD_C2P_HEADER
#define SRMHD_C2P_HEADER

#ifdef __cplusplus
extern "C" {
#endif
  typedef struct srmhd_c2p srmhd_c2p;

  enum srmhd_c2p_error_code {
    SRMHD_C2P_SUCCESS,
    SRMHD_C2P_CONS_CONTAINS_NAN,
    SRMHD_C2P_CONS_NEGATIVE_DENSITY,
    SRMHD_C2P_CONS_NEGATIVE_ENERGY,
    SRMHD_C2P_PRIM_CONTAINS_NAN,
    SRMHD_C2P_PRIM_NEGATIVE_PRESSURE,
    SRMHD_C2P_PRIM_NEGATIVE_RESTMASS,
    SRMHD_C2P_PRIM_SUPERLUMINAL,
    SRMHD_C2P_MAXITER
  } ;

  srmhd_c2p *srmhd_c2p_new();
  void srmhd_c2p_del(srmhd_c2p *c2p);
  void srmhd_c2p_set_gamma(srmhd_c2p *c2p, double adiabatic_gamma);
  void srmhd_c2p_get_gamma(srmhd_c2p *c2p, double *adiabatic_gamma);
  void srmhd_c2p_new_state(srmhd_c2p *c2p, const double *U);
  void srmhd_c2p_estimate_from_cons(srmhd_c2p *c2p);
  void srmhd_c2p_set_starting_prim(srmhd_c2p *c2p, const double *P);
  void srmhd_c2p_get_starting_prim(srmhd_c2p *c2p, double *P);
  int srmhd_c2p_set_pressure_floor(srmhd_c2p *c2p, double pf);
  int srmhd_c2p_put_pressure_floor(srmhd_c2p *c2p);
  int srmhd_c2p_reconstruct_prim(srmhd_c2p *c2p, double Z, double W, double *P);
  int srmhd_c2p_solve_anton2dzw(srmhd_c2p *c2p, double *P);
  int srmhd_c2p_solve_noble1dw(srmhd_c2p *c2p, double *P);
  int srmhd_c2p_get_iterations(srmhd_c2p *c2p);
  int srmhd_c2p_check_cons(srmhd_c2p *c2p, const double *U);
  int srmhd_c2p_check_prim(srmhd_c2p *c2p, const double *P);
  const char *srmhd_c2p_get_error(srmhd_c2p *c2p, int error);

#ifdef __cplusplus
}
#endif

#endif /* SRMHD_C2P_HEADER */
