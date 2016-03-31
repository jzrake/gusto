#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gusto.h"


void gusto_free(struct gusto_sim *sim)
{
  gusto_mesh_clear(sim);
}



/*
 * Main function
 * =====================================================================
 */
int main1(int argc, char **argv)
{
  struct gusto_sim sim;

  gusto_user_set_defaults(&sim.user);
  gusto_status_set_defaults(&sim.status);


  for (int n=1; n<argc; ++n) {
    gusto_user_set_from_arg(&sim.user, argv[n]);
  }


  gusto_user_report(&sim.user);
  gusto_mesh_generate_verts(&sim);
  gusto_mesh_generate_cells(&sim);
  gusto_mesh_generate_faces(&sim);
  gusto_mesh_report(&sim);
  gusto_mesh_compute_geometry(&sim);
  gusto_write_checkpoint(&sim, NULL);

  gusto_free(&sim);

  return 0;
}




/*
 * Main function
 * =====================================================================
 */
int main(int argc, char **argv)
{
  struct gusto_sim sim;
  gusto_user_set_defaults(&sim.user);
  gusto_status_set_defaults(&sim.status);


  for (int n=1; n<argc; ++n) {
    gusto_user_set_from_arg(&sim.user, argv[n]);
  }


  gusto_user_report(&sim.user);
  gusto_mesh_generate_verts(&sim);
  gusto_mesh_generate_cells(&sim);
  gusto_mesh_generate_faces(&sim);
  gusto_mesh_compute_geometry(&sim);
  gusto_initial_data(&sim);


  double dt = 0.25 * sim.smallest_cell_length;
  sim.status.time_step = dt;

  while (sim.status.time_simulation < sim.user.tmax) {

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

    gusto_compute_variables_at_vertices(&sim);
    gusto_compute_vertex_velocities(&sim);
    gusto_compute_fluxes(&sim);
    gusto_transmit_fluxes(&sim, dt);
    gusto_mesh_advance_vertices(&sim, dt);
    gusto_mesh_compute_geometry(&sim);
    gusto_recover_variables(&sim);
    gusto_enforce_boundary_condition(&sim);

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

  gusto_free(&sim);
  return 0;
}
