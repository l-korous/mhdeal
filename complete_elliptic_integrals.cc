////////////////////////////////////////////////////////////////////////////////
// File: complete_elliptic_integrals.c                                        //
// Routine(s):                                                                //
//    Complete_Elliptic_Integral_First_Kind                                   //
//    Complete_Elliptic_Integral_Second_Kind                                  //
////////////////////////////////////////////////////////////////////////////////

#include"complete_elliptic_integrals.h"

////////////////////////////////////////////////////////////////////////////////
// double Complete_Elliptic_Integral_First_Kind(char arg, double x)           //
//                                                                            //
//  Description:                                                              //
//     The complete elliptic integral of the first kind is the integral from  //
//     0 to pi / 2 of the integrand                                           //
//                   dtheta / sqrt( 1 - k^2 sin^2(theta) ).                   //
//     The parameter k is called the modulus.  This integral is even in k.    //
//     The modulus, k, must satisfy |k| <= 1.  If k = 0 then the integral     //
//     can be readily evaluated.  If |k| = 1, then the integral is infinite.  //
//     Otherwise it must be approximated.                                     //
//                                                                            //
//     In practise the arguments the complete elliptic function of the first  //
//     kind are also given as F(pi/2 \ alpha) or F(pi/2 | m) where the angle  //
//     alpha, called the modular angle, satisfies k = sin(alpha) and the      //
//     argument m = k^2 is simply called the parameter.                       //
//     In terms of these arguments K = F(pi/2 \ alpha) = F(pi/2, sin(alpha))  //
//     and K = F(pi/2 | m) = F(pi/2, sqrt(m)), where                          //
//             K = Complete_Elliptic_Integral_First_Kind( k ).                //
//                                                                            //
//     Let K(k) be the complete elliptic integral of the second kind where    //
//     k is the modulus and  k' = sqrt(1-k^2) is the complementary modulus.   //
//                                                                            //
//     The common mean method, sometimes called the Gauss transform, is a     //
//     variant of the descending Landen transformation in which two sequences //
//     are formed: Setting a[0] = 1 and g[0] = k', the complementary modulus, //
//     a[i] is the arithmetic average and g[i] is the geometric mean of a[i-1]//
//     and g[i-1], i.e. a[i+1] = (a[i] + g[i])/2 and g[i+1] = sqrt(a[i]*g[i]).//
//     The sequences satisfy the inequalities g[0] < g[1] < ... < a[1] < a[0].//
//     Further, lim g[n] = lim a[n].                                          //
//     The value of the complete elliptic integral of the first kind is       //
//     (pi/2) lim (1/G[n]) as n -> inf.                                       //
//                                                                            //
//  Arguments:                                                                //
//     char    arg                                                            //
//                The type of argument of the second argument of F():         //
//                  If arg = 'k', then x = k, the modulus of F(pi/2,k).       //
//                  If arg = 'a', then x = alpha, the modular angle of        //
//                                F(pi/2 \ alpha), alpha in radians.          //
//                  If arg = 'm', then x = m, the parameter of F(pi/2 | m).   //
//                  The value of arg defaults to 'k'.                         //
//     double  x                                                              //
//                The second argument of the elliptic function F(pi/2,k),     //
//                F(pi/2 \ alpha) or F(pi/2 | m) corresponding to the value   //
//                of 'arg'.  Note that if arg = 'k', then | x | <= 1 and if   //
//                arg = 'm', then 0 <= x <= 1.                                //
//                                                                            //
//  Return Value:                                                             //
//     The value of the complete elliptic integral of the first kind for the  //
//     given modulus, modular angle, or parameter.  Note that if |k| = 1,     //
//     or m = 1 or a = (+/-) pi/2 then the integral is infinite and DBL_MAX   //
//     is returned.                                                           //
//                                                                            //
//  Example:                                                                  //
//     double K;                                                              //
//     double m, k, a;                                                        //
//                                                                            //
//     ( code to initialize a )                                               //
//                                                                            //
//     k = sin(a);                                                            //
//     m = k * k;                                                             //
//     K = Complete_Elliptic_Integral_First_Kind( 'a', a );                   //
//     printf("K(alpha) = %12.6f where angle(radians) = %12.6f\n",K, a);      //
//     K = Complete_Elliptic_Integral_First_Kind( 'k', k );                   //
//     printf("K(k) = %12.6f where k = %12.6f\n",K, k);                       //
//     K = Complete_Elliptic_Integral_First_Kind( 'm', m );                   //
//     printf("K(m) = %12.6f where m = %12.6f\n",K, m);                       //
////////////////////////////////////////////////////////////////////////////////

double Complete_Elliptic_Integral_First_Kind(char arg, double x)
{
   long double k;          // modulus
   long double m;          // parameter
   long double a;          // average
   long double g;          // geometric mean
   long double a_old;      // previous average
   long double g_old;      // previous geometric mean

   if ( x == 0.0 ) return M_PI_2;

   switch (arg) {
      case 'k': k = fabsl((long double) x);
                m = k * k;
                break;
      case 'm': m = (long double) x;
                k = sqrtl(fabsl(m));
                break;
      case 'a': k = sinl((long double)x);
                m = k * k;
                break;
      default:  k = fabsl((long double) x);
                m = k * k;
   }

   if ( m == 1.0 ) return DBL_MAX;

   a = 1.0L;
   g = sqrtl(1.0L - m);
   while (1) {
      g_old = g;
      a_old = a;
      a = 0.5L * (g_old + a_old);
      g = sqrtl(g_old * a_old);
      if ( fabsl(a_old - g_old) <= (a_old * LDBL_EPSILON) ) break;
   }
   return (double) (PI_2 / g);
}

////////////////////////////////////////////////////////////////////////////////
// double Complete_Elliptic_Integral_Second_Kind(cahr arg, double x)          //
//                                                                            //
//  Description:                                                              //
//     The complete elliptic integral of the second kind is the integral from //
//     0 to pi / 2 of the integrand                                           //
//                    sqrt( 1 - k^2 sin^2(theta) ) dtheta .                   //
//     The parameter k is called the modulus.  This integral is even in k.    //
//     The modulus, k, must satisfy |k| <= 1.  If k = 0 or |k| = 1 then the   //
//     integral can be readily evaluated.  Otherwise it must be approximated. //
//                                                                            //
//     In practise the arguments the elliptic function of the second kind are //
//     also given as E(pi/2 \ alpha) or E(pi/2 | m) where the angle alpha,    //
//     called the modular angle, satisfies k = sin(alpha) and the argument    //
//     m = k^2 is simply called the parameter.                                //
//     In terms of these arguments E = E(pi/2 \ alpha) = E(pi/2, sin(alpha))  //
//     and E = E(pi/2 | m) = E(pi/2, sqrt(m)), where                          //
//             E = Complete_Elliptic_Integral_Second_Kind( k ).               //
//                                                                            //
//     Let K(k) be the complete elliptic integral of the second kind where    //
//     k is the modulus and  k' = sqrt(1-k^2) is the complementary modulus.   //
//                                                                            //
//     The common mean method, sometimes called the Gauss transform, is a     //
//     variant of the descending Landen transformation in which two sequences //
//     are formed: Setting a[0] = 1 and g[0] = k', the complementary modulus, //
//     a[i] is the arithmetic average and g[i] is the geometric mean of a[i-1]//
//     and g[i-1], i.e. a[i+1] = (a[i] + g[i])/2 and g[i+1] = sqrt(a[i]*g[i]).//
//     The sequences satisfy the inequalities g[0] < g[1] < ... < a[1] < a[0].//
//     Further, lim g[n] = lim a[n] as n -> inf.                              //
//     The value of the complete elliptic integral of the second kind is      //
//           E(k) = lim (pi/8g[n]) (4 - 2k^2 - Sum(2^j(a[j]-g[j])^2)).        //
//     where the limit is as n -> inf and the sum extends from j = 0 to n.    //
//     The sum of 2^j (a[j]^2 - g[j]^2) from j = 1 to n equals                //
//           (1/2) Sum (2^i (a[i] - g[i])^2 for i = 0,...,n-1, so that        //
//           E(k) = lim (pi/4g[n]) (2 - k^2 - Sum(2^j(a[j]^2 -g[j]^2))).      //
//                                                                            //
//  Arguments:                                                                //
//     char    arg                                                            //
//                The type of argument of the second argument of E():         //
//                  If arg = 'k', then x = k, the modulus of E(pi/2,k).       //
//                  If arg = 'a', then x = alpha, the modular angle of        //
//                                E(pi/2 \ alpha), alpha in radians.          //
//                  If arg = 'm', then x = m, the parameter of E(pi/2 | m).   //
//                  The value of arg defaults to 'k'.                         //
//     double  x                                                              //
//                The second argument of the elliptic function E(pi/2,k),     //
//                E(pi/2 \ alpha) or E(pi/2 | m) corresponding to the value   //
//                of 'arg'.  Note that if arg = 'k', then | x | <= 1 and if   //
//                arg = 'm', then 0 <= x <= 1.                                //
//                                                                            //
//  Return Value:                                                             //
//     The value of the complete elliptic integral of the second kind for the //
//     given modulus, modular angle, or parameter.                            //
//                                                                            //
//  Example:                                                                  //
//     double E;                                                              //
//     double m, k, a;                                                        //
//                                                                            //
//     ( code to initialize a )                                               //
//                                                                            //
//     k = sin(a);                                                            //
//     m = k * k;                                                             //
//     E = Complete_Elliptic_Integral_Second_Kind( 'a', a );                  //
//     printf("E(alpha) = %12.6f where angle(radians) = %12.6f\n",E, a);      //
//     E = Complete_Elliptic_Integral_Second_Kind( 'k', k );                  //
//     printf("E(k) = %12.6f where k = %12.6f\n",E, k);                       //
//     E = Complete_Elliptic_Integral_Second_Kind( 'm', m );                  //
//     printf("E(m) = %12.6f where m = %12.6f\n",E, m);                       //
////////////////////////////////////////////////////////////////////////////////

double Complete_Elliptic_Integral_Second_Kind(char arg, double x)
{
   long double k;      // modulus
   long double m;      // the parameter of the elliptic function m = modulus^2
   long double a;      // arithmetic mean
   long double g;      // geometric mean
   long double a_old;  // previous arithmetic mean
   long double g_old;  // previous geometric mean
   long double two_n;  // power of 2
   long double Ek;

   if ( x == 0.0 ) return M_PI_2;

   switch (arg) {
      case 'k': k = fabsl((long double) x);
                m = k * k;
                break;
      case 'm': m = (long double) x;
                k = sqrtl(fabsl(m));
                break;
      case 'a': k = sinl((long double)x);
                m = k * k;
                break;
      default:  k = fabsl((long double) x);
                m = k * k;
   }

   if ( m == 1.0 ) return 1.0;

   a = 1.0L;
   g = sqrtl(1.0L - m);
   two_n = 1.0L;
   Ek = 2.0L - m;

   while (1) {
      g_old = g;
      a_old = a;
      a = 0.5L * (g_old + a_old);
      g = a_old * g_old;
      two_n += two_n;
      Ek -= two_n * (a * a - g);
      if ( fabsl(a_old - g_old) <= (a_old * LDBL_EPSILON) ) break;
      g = sqrtl(g);
   }

   return (double) ((PI_4 / a) * Ek);
}

void Complete_Elliptic_Integrals(char arg, double x, double& K, double& E)
{
   long double k;          // modulus
   long double m;          // parameter
   long double a;          // average
   long double g;          // geometric mean
   long double a_old;      // previous average
   long double g_old;      // previous geometric mean
   long double two_n;      // power of 2

   if ( x == 0.0 ){
     K = E = M_PI_2;
     return;
   }

   switch (arg) {
      case 'k': k = (long double) x;
                m = k * k;
                break;
      case 'm': m = (long double) x;
                k = sqrtl(fabsl(m));
                break;
      case 'a': k = sinl((long double)x);
                m = k * k;
                break;
      default:  k = (long double) x;
                m = k * k;
   }

   if ( m == 1.0 ){
    K = DBL_MAX;
    E = 1.0;
    return;
   }

   a = 1.0L;
   g = sqrtl(1.0L - m);
   two_n = 1.0L;
   E = 2.0L - m;
   while(true){
      g_old = g;
      a_old = a;
      a = 0.5L * (g_old + a_old);
      g = g_old * a_old;
      two_n += two_n;
      E -= two_n * (a * a - g);
      g = sqrtl(g);
      if ( fabsl(a_old - g_old) <= (a_old * LDBL_EPSILON) ) break;
   }

   K = (double) (PI_2 / g);
   E = (double) ((PI_4 / a) * E);
}

void Complete_Elliptic_Integrals_Modulus(double x, double& K, double& E)
{
   long double k;          // modulus
   long double m;          // parameter
   long double a;          // average
   long double g;          // geometric mean
   long double a_old;      // previous average
   long double g_old;      // previous geometric mean
   long double two_n;      // power of 2
   long double cE;         // sum

   if ( x == 0.0 ){
     K = E = M_PI_2;
     return;
   }

   k = (long double) x;
   m = k * k;

   if ( x == 1.0 ){
    K = DBL_MAX;
    E = 1.0;
    return;
   }

   a = 1.0L;
   g = sqrtl(1.0L - m);
   two_n = 1.0L;
   cE = 2.0L - m;
   while(true){
      g_old = g;
      a_old = a;
      a = 0.5L * (g_old + a_old);
      g = g_old * a_old;
      two_n += two_n;
      cE -= two_n * (a * a - g);
      g = sqrtl(g);
      if ( fabsl(a_old - g_old) <= (a_old * LDBL_EPSILON) ) break;
   }

   K = (double) (PI_2 / g);
   E = (double) ((PI_4 / a) * cE);
}

