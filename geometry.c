#include "gusto.h"


void gusto_geometry(struct aux_geometry *G, double x[4])
{
  double r = x[1];
  double B = 1.0 / (r * r);
  G->cylindrical_radius = r;
  G->volume_element = 1 / B;
  G->area_element[1] = r;
  G->area_element[2] = 1 / (B * r);
  G->area_element[3] = 1 / B;
  G->line_element[1] = 1 / (B * r);
  G->line_element[2] = r;
  G->line_element[3] = 1;
}
