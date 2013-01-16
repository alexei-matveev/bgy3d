/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

void bgy3d_vec_map1 (Vec ys, real (*f)(real x), Vec xs);
void bgy3d_vec_map2 (Vec zs, real (*f)(real x, real y), Vec xs, Vec ys);

void bgy3d_vec_save (const char file[], const Vec vec);
Vec bgy3d_vec_load (const char file[]); /* Creates a new Vec */
void bgy3d_vec_read (const char file[], Vec vec);  /* Fills existing Vec */
void bgy3d_vec_save_ascii (const char file[], const Vec vec);


