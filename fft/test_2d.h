/* PREC = 1 is single-precision complex (8-byte datum) */
/* PREC = 2 is double-precision complex (16-byte datum) */

#define PREC 2

#if PREC == 1

#define COMPLEX_DATA complex*8
#define REAL_DATA real*4
#define COMPLEX cmplx
#define REAL real
#define IMAG aimag
#define FLOAT real

#endif

#if PREC == 2

#define COMPLEX_DATA complex*16
#define REAL_DATA real*8
#define COMPLEX dcmplx
#define REAL dble
#define IMAG dimag
#define FLOAT dble

#endif
