/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

void bgy3d_vec_map1 (Vec ys, real (*f)(real x), Vec xs);
void bgy3d_vec_map2 (Vec zs, real (*f)(real x, real y), Vec xs, Vec ys);

real bgy3d_vec_mix (Vec dg, Vec dg_new, real a, Vec work);

void bgy3d_vec_save (const char file[], const Vec vec);
void bgy3d_vec_save1 (const char *format, int m, const Vec g[m]);
void bgy3d_vec_save2 (const char *format, int m, /* const */ Vec g2[m][m]);

Vec bgy3d_vec_load (const char file[]); /* Creates a new Vec */

void bgy3d_vec_read (const char file[], Vec vec);  /* Fills existing Vec */
void bgy3d_vec_read1 (const char *format, int m, const Vec g[m]);
void bgy3d_vec_read2 (const char *format, int m, /* const */ Vec g2[m][m]);

void bgy3d_vec_save_ascii (const char file[], const Vec vec);
void bgy3d_vec_save_ascii1 (const char *format, int m, const Vec vec[m]);
void bgy3d_vec_save_ascii2 (const char *format, int m, /* const */ Vec vec[m][m]);


void bgy3d_vec_read_radial2 (const State *BHD,
                             const char *format, int m, /* const */ Vec g2[m][m]);


