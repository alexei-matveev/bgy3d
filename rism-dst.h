/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2013 Alexei Matveev
*/

void rism_dst (size_t n, double out[n], const double in[n]);
void rism_dst_columns (int m, int n, double buf[m][n]);
void rism_dst_rows (int n, int m, double buf[n][m]);
