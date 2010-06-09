/*==========================================================*/
/*  $Id: bgy3d_SolventParameters.h,v 1.2 2007-07-31 17:12:33 jager Exp $ */
/*==========================================================*/



/******************************************/
/*  2-site models */
/******************************************/

#ifdef N2
/*===================================================*/
/* Parametrisierung N2 nach Hirata JChemPhys(77) p.509:  */
/* Density : 0.01867 */
/* T= 72 K : beta = 6.9938 */
/* T=125 K : beta = 4.0284 */

#define sH  3.341
#define eH  0.08738
#define qH  0.2
#define sO  3.341
#define eO  0.08738
#define qO  -0.2

#define r_HH  -1
#define r_HO  1.1



#define EPSILON0INV 331.84164 //331.84164

/* smoothing parameter */
#define G   1.2 //1.2  



#endif

#ifdef HCl
/*===================================================*/
/* Parametrisierung HCl nach Hirata JChemPhys(77) p.509:  */
/* Density : 0.018 */
/* T= 210 K : beta = 2.39788 */
/* T= 420 K : beta = 1.1989 */
/* T= 315 K : beta = 1.5985 */

#define sH  2.735 //0.4
#define eH  0.03971
#define qH  0.2
#define sO  3.353
#define eO  0.51434
#define qO  -0.2

#define r_HH -1
#define r_HO  1.257 //1.3

#define EPSILON0INV 331.84164 //331.84164

/* smoothing parameter */
#define G   1.2 //1.2  



#endif




/******************************************/
/*  3-site models */
/******************************************/


#ifdef H2O

/*===================================================*/
/* Parametrisierung TIP3P :  */
/* sigma_H = 0.400   epsilon_H = 0.046  q_H = 0.417*/
/* sigma_O = 3.1506  epsilon_O = 0.1521 q_O =-0.834 */
/* r_OH = 0.9572 */
/* r_HH = 1.5139 */
/* theta_HOH = 104.52 */
/* mass H2O = 18.0154 u */
/* density = 1 kg/l = 0.6022142 u/A^3 => 0.033427745 / A^3 */
/* temperature : T= 298,15 K (25 C) => 0.5921 , => beta =1.6889 */

/* EPSILON0INV */
/* You have: e^2/4/pi/epsilon0/angstrom */
/* You want: kcal/avogadro/mol */
/* => 331.84164 */

#define sH 0.4 //2.5 //0.4 //0.400
#define eH 0.046 //0.046
#define qH 0.417
#define sO 3.1506 //2.8509 //3.1506
#define eO 0.1521 //0.1521
#define qO -0.834

#define r_HH  1.5139
#define r_HO  0.9572

#define EPSILON0INV 331.84164 //331.84164

/* smoothing parameter */
#define G   1.2 //1.2  

#endif

#ifdef CS2

/*===================================================*/
/* Parametrisierung CS2 nach Zhu:  H=S O=C*/
/* sigma_H = 3.520   epsilon_H = 0.39500  q_H = 0.154*/
/* sigma_O = 3.200  epsilon_O = 0.10128 q_O =-0.308 */
/* r_OH = 1.56 */
/* r_HH = 2*1.56 */
/* theta_HOH = 180.0 */
/* mass CS2 =  75.999432 u */
/* density = 1.263 g·cm^-3 =  0.76059653 u/A^3 => 0.010007924 / A^3 */
/* temperature : T= 298,15 K (25 C) => 0.5921 , => beta =1.6889 */
/* temperature : T= 360 K (? C) => 0.7149 , => beta =1.3988 */

/* EPSILON0INV */
/* You have: e^2/4/pi/epsilon0/angstrom */
/* You want: kcal/avogadro/mol */
/* => 331.84164 */

#define sH 3.520
#define eH 0.39500
#define qH 0.154
#define sO 3.200
#define eO 0.10128
#define qO -0.308

#define r_HH  3.12
#define r_HO  1.56

#define EPSILON0INV 331.84164 //331.84164

/* smoothing parameter */
#define G   0.8 //1.2  


#endif

#ifdef CS2_II

/*===================================================*/
/* Parametrisierung CS2 nach Tildesley:  H=S O=C*/
/* sigma_H = 3.520   epsilon_H = 0.36341  q_H = 0*/
/* sigma_O = 3.35  epsilon_O = 0.10168 q_O =0 */
/* r_OH = 1.57 */
/* r_HH = 2*1.57 */
/* theta_HOH = 180.0 */
/* mass CS2 =  75.999432 u */
/* density = 1.263 g·cm^-3 =  0.76059653 u/A^3 => 0.010007924 / A^3 */
/* temperature : T= 298,15 K (25 C) => 0.5921 , => beta =1.6889 */

/* EPSILON0INV */
/* You have: e^2/4/pi/epsilon0/angstrom */
/* You want: kcal/avogadro/mol */
/* => 331.84164 */

#define sH 3.520
#define eH 0.36341
#define qH 0.0
#define sO 3.35
#define eO 0.10168
#define qO -0.0

#define r_HH  3.14
#define r_HO  1.57

#define EPSILON0INV 331.84164 //331.84164

/* smoothing parameter */
#define G   0.8 //1.2 


#endif

#ifdef PROPANE

/*===================================================*/
/* Parametrisierung C3H8 nach Martin:  H=CH3 O=CH2*/
/* sigma_H = 3.75   epsilon_H = 0.194616  q_H = 0*/
/* sigma_O = 3.95  epsilon_O = 0.091350   q_O =0 */
/* r_OH = 1.57 */
/* r_HH = 2*1.57 */
/* theta_HOH = 114.0 */
/* mass C3H8 =  44.0099 u */
/* density = 0.5077 kg/L =  0.3011071 u/A^3 => 0.0068418038 /A^3 */
/* temperature : T= 200 K  => 0.3971755 , => beta =2.5178 */

/* EPSILON0INV */
/* You have: e^2/4/pi/epsilon0/angstrom */
/* You want: kcal/avogadro/mol */
/* => 331.84164 */

#define sH 3.75
#define eH 0.194616
#define qH 0.0
#define sO 3.95
#define eO 0.091350
#define qO -0.0

#define r_HH  2.5831
#define r_HO  1.54

#define EPSILON0INV 331.84164 //331.84164

/* smoothing parameter */
#define G   0.8 //1.2 


#endif


/******************************************/
/*  4-site models */
/******************************************/

#ifdef H2O2

/*===================================================*/
/* Parametrisierung TIP3P :  */
/* sigma_H = 0.400   epsilon_H = 0.046  q_H = 0.417*/
/* sigma_O = 3.1506  epsilon_O = 0.1521 q_O =-0.834 */
/* r_OH = 0.9572 */
/* r_HH = 1.5139 */
/* theta_HOH = 104.52 */
/* mass H2O = 18.0154 u */
/* density = 1 kg/l = 0.6022142 u/A^3 => 0.033427745 / A^3 */
/* temperature : T= 298,15 K (25 C) => 0.5921 , => beta =1.6889 */

/* EPSILON0INV */
/* You have: e^2/4/pi/epsilon0/angstrom */
/* You want: kcal/avogadro/mol */
/* => 331.84164 */

#define sH 0.4 //0.400
#define eH 0.046 //0.046
#define qH 0.350
#define sO 3.1506 //2.8509 //3.1506
#define eO 0.1521 //0.1521
#define qO -0.350

#define r_HH  1.5139
#define r_OO  1.453
#define r_HO  0.988
#define r_HO2 

#define EPSILON0INV 331.84164 //331.84164

/* smoothing parameter */
#define G   1.2 //1.2  

#endif
