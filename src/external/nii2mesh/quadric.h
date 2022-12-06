#ifndef QUADRIC_H
#define QUADRIC_H

#include "meshtypes.h"

void quadric_simplify_mesh(vec3d **vs, vec3i **ts, int* nvert, int *ntri, int target_count, double agressiveness, bool verbose, bool finishLossless);
void laplacian_smoothHC(vec3d *verts, vec3i *tris, int nvert, int ntri, double alpha, double beta, int iter, bool lockEdges);

#endif /* QUADRIC_H */