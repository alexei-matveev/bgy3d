/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

void bgy3d_vec_map1 (Vec ys, real (*f)(real x), Vec xs);
void bgy3d_vec_map2 (Vec zs, real (*f)(real x, real y), Vec xs, Vec ys);

