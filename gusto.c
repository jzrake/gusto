#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gusto.h"



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
  gusto_mesh_compute_geometry(&sim);
  gusto_write_checkpoint(&sim, NULL);

  gusto_mesh_clear(&sim);

  return 0;
}
