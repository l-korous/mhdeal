////////////////////////////////////////////////////////////////////////////////
// File: complete_elliptic_integrals.h                                        //
// Routine(s):                                                                //
//    Complete_Elliptic_Integral_First_Kind                                   //
//    Complete_Elliptic_Integral_Second_Kind                                  //
//    Complete_Elliptic_Integral                                              //
//    Complete_Elliptic_Integral_Modulus                                      //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>       // required for fabs(), fabsl(), sqrtl(), and M_PI_2
#include <float.h>      // required for LDBL_EPSILON, DBL_MAX

static const long double PI_2 =  1.5707963267948966192313216916397514L; // pi/2
static const long double PI_4 = 0.7853981633974483096156608458198757L; // pi/4

////////////////////////////////////////////////////////////////////////////////
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
//                                                                            //
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
////////////////////////////////////////////////////////////////////////////////

double Complete_Elliptic_Integral_First_Kind(char, double);
double Complete_Elliptic_Integral_Second_Kind(char, double);
void Complete_Elliptic_Integrals(char, double, double&, double&);
void Complete_Elliptic_Integrals_Modulus(double, double&, double&);