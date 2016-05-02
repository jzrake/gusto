#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gusto.h"


void gusto_init(struct gusto_sim *sim)
{
  sim->num_rows = 0;
  sim->smallest_cell_length = 0.0;
  sim->rows = NULL;
  gusto_user_set_defaults(&sim->user);
  gusto_status_set_defaults(&sim->status);
  sim->boundary_con = NULL;
  sim->initial_mesh = NULL;
  sim->initial_data = NULL;
}


void gusto_free(struct gusto_sim *sim)
{
  gusto_mesh_clear(sim, 'a');
}


void gusto_config_from_user(struct gusto_sim *sim)
{
  sim->boundary_con = gusto_lookup_boundary_con(sim->user.boundary_con);
  sim->initial_data = gusto_lookup_initial_data(sim->user.initial_data);
  /* sim->initial_mesh = gusto_lookup_initial_mesh(sim->user.initial_mesh); */
}


int gusto_validate_sim(struct gusto_sim *sim)
{
  int valid = 1;
  valid &= sim->boundary_con != NULL;
  valid &= sim->initial_data != NULL;
  /* valid &= sim->initial_mesh != NULL; */
  return valid;
}


void gusto_print_help(struct gusto_sim *sim)
{
  const char **help_lines = NULL;
  if (sim->initial_data) {
    help_lines = sim->initial_data(&sim->user, NULL, NULL);
  }
  if (help_lines) {
    int n = 0;
    printf("\n");
    while (help_lines[n]) {
      printf("\t%s\n", help_lines[n]);
      n += 1;
    }
    printf("\n");
  }
}



void gusto_advance_rk(struct gusto_sim *sim, double dt, double rkparam)
{
  /*
   * Evaluate derived data associated with the input state, namely Godunov
   * fluxes.
   * ---------------------------------------------------------------------------
   */
  gusto_compute_fluxes(sim);

  /*
   * Update conserved quantities, vector potential, and vertex locations.
   * ---------------------------------------------------------------------------
   */
  gusto_transmit_fluxes(sim, dt);
  gusto_add_source_terms(sim, dt);

  /*
   * Correct conserved quantities, vector potential, and vertex locations for
   * Runge-Kutta parameter.
   * ---------------------------------------------------------------------------
   */
  gusto_average_rk(sim, rkparam);

  /*
   * Recover primitive variables from the new conserved variables and vector
   * potential, enforce boundary condtitions.
   * ---------------------------------------------------------------------------
   */
  gusto_recover_variables(sim);
  sim->boundary_con(sim);
}



/*
 * Main function
 * =====================================================================
 */
int main(int argc, char **argv)
{
  int restarted_run = 0;

  struct gusto_sim sim;
  gusto_init(&sim);

  for (int n=1; n<argc; ++n) {
    gusto_user_set_from_arg(&sim.user, argv[n]);
  }


  if (restarted_run) {

    /* put restart code here */

  }

  else {

    sim.status.time_last_checkpoint = -sim.user.cpi;
    sim.status.checkpoint_number = -1;

  }


  gusto_user_report(&sim.user);
  gusto_config_from_user(&sim);
  gusto_print_help(&sim);

  if (!gusto_validate_sim(&sim)) {
    printf("[gusto] ERROR: invalid setup\n");
    return 0;
  }

  gusto_mesh_generate(&sim);
  gusto_mesh_compute_geometry(&sim);
  gusto_mesh_report(&sim);

  gusto_initial_data(&sim);
  /* gusto_validate_fluxes(&sim); */
  gusto_validate_geometry(&sim);


  while (sim.status.time_simulation < sim.user.tmax) {

    double dt = sim.user.cfl * sim.smallest_cell_length;
    sim.status.time_step = dt;

    /*
     * Write a checkpoint if it's time
     * =================================================================
     */
    if (sim.status.time_simulation - sim.status.time_last_checkpoint >=
	sim.user.cpi && sim.user.cpi > 0.0) {

      sim.status.time_last_checkpoint += sim.user.cpi;
      sim.status.checkpoint_number += 1;

      gusto_write_checkpoint(&sim, NULL);
    }

    void *start_cycle = gusto_start_clock();

    gusto_cache_rk(&sim);

    switch (sim.user.rk_order) {
    case 1:
      gusto_advance_rk(&sim, dt, 1.);
      break;
    case 2:
      gusto_advance_rk(&sim, dt, 1.);
      gusto_advance_rk(&sim, dt, 1./2);
      break;
    case 3:
      gusto_advance_rk(&sim, dt, 1.);
      gusto_advance_rk(&sim, dt, 1./4);
      gusto_advance_rk(&sim, dt, 2./3);
      break;
    }

    double seconds = gusto_stop_clock(start_cycle);

    sim.status.kzps = 1e-3 * gusto_mesh_count(&sim, 'c', -1) / seconds;

    if (sim.status.iteration % 1 == 0) {
      printf("[gusto] n=%06d t=%6.4e dt=%6.4e %3.2f kzps\n",
	     sim.status.iteration,
	     sim.status.time_simulation,
	     sim.status.time_step,
	     sim.status.kzps);
    }

    sim.status.time_simulation += dt;
    sim.status.iteration += 1;
  }

  gusto_write_checkpoint(&sim, "chkpt.final.h5");
  gusto_free(&sim);
  return 0;
}
