/*------------------------------------------------------------------------------
 * FILE: quartic.c
 *
 * AUTHOR: Jonathan Zrake, KIPAC (adapted from Mathematica code)
 *
 * DESCRIPTION: Implementation of exact solution to the quartic polynomial
 *
 *------------------------------------------------------------------------------
 */

#include <math.h>
#include <complex.h>

static int do_root_polish = 0;
static void polish_newton(double as[4], double roots[4], int good[4]);

int solve_quartic_equation(double d4, double d3,
			   double d2, double d1, double d0,
			   double roots[4])
{
  typedef double complex Complex;

  Complex a3 = d3 / d4;
  Complex a2 = d2 / d4;
  Complex a1 = d1 / d4;
  Complex a0 = d0 / d4;
  Complex W = 1./3;
  Complex X = 12*a0 - 3*a1*a3 + a2*a2;
  Complex P = -72*a0*a2 - 9*a1*a2*a3 + 27*a1*a1 + 27*a0*a3*a3 + 2*a2*a2*a2;
  Complex Q = cpow(X,3);
  Complex S = csqrt(-4*Q + cpow(P,2));
  Complex T = -8*a1 + 4*a2*a3 - a3*a3*a3;
  Complex B = cabs(P + S) == 0.0 ? 0.0 : (cpow(2,W)*X)/(3.*cpow(P + S,W));
  Complex U = (-2*a2)/3. + (a3*a3)/4. + B;
  Complex C = csqrt(U + cpow(P + S,W)/(3.*cpow(2,W)))/2.;
  Complex D = cpow(P + S,W)/(3.*cpow(2,W));
  Complex E = T/(4.*csqrt(U + D));
  Complex F = csqrt((-4*a2)/3. + (a3*a3)/2. - B - D - E)/2.;
  Complex G = csqrt((-4*a2)/3. + (a3*a3)/2. - B - D + E)/2.;
  Complex r0 = -a3/4. - C - F;
  Complex r1 = -a3/4. - C + F;
  Complex r2 = -a3/4. + C - G;
  Complex r3 = -a3/4. + C + G;

  roots[0] = creal(r0);
  roots[1] = creal(r1);
  roots[2] = creal(r2);
  roots[3] = creal(r3);

  if (do_root_polish) {
    double as[4] = {a0, a1, a2, a3};
    int good[4];
    polish_newton(as, roots, good);
    return good[0] + good[1] + good[2] + good[3];
  }
  else {
    return 4;
  }
}

void polish_newton(double as[4], double roots[4], int good[4])
{
  int n, i;
  double f, g, h, x;
  for (n=0; n<4; ++n) {
    i = 0;
    x = roots[n];
    good[n] = 1;
    do {
      i += 1;
      f = x*x*x*x + as[3]*x*x*x + as[2]*x*x + as[1]*x + as[0];
      g = x*x*x*4 + as[3]*x*x*3 + as[2]*x*2 + as[1]*1;
      if (g != 0.0) {
	x -= f/g;
      }
      if (x != 0.0) {
	h = f / x;
      }
      else {
	h = f;
      }
      if (i == 100) { /* root is not a root */
	x = 0.0;
	good[n] = 0;
	break;
      }
    } while (fabs(h) > 1e-14);
    roots[n] = x;
  }
}


#ifdef __MAIN__
#include <stdio.h>
void test_roots(double a, double b, double c, double d)
{
  double d0 = a*b*c*d;
  double d1 = -a*b*c - a*b*d - a*c*d - b*c*d;
  double d2 = a*b + a*c + b*c + a*d + b*d + c*d;
  double d3 = -a - b - c - d;
  double d4 = 1.0;

  double r[4];
  solve_quartic_equation2(d4, d3, d2, d1, d0, r);

  printf("r0 = [%+12.10e %+12.10e] error=%6.4e\n", a, r[0], fabs((a - r[0])/a));
  printf("r1 = [%+12.10e %+12.10e] error=%6.4e\n", b, r[1], fabs((b - r[1])/b));
  printf("r2 = [%+12.10e %+12.10e] error=%6.4e\n", c, r[2], fabs((c - r[2])/c));
  printf("r3 = [%+12.10e %+12.10e] error=%6.4e\n", d, r[3], fabs((d - r[3])/d));
}
int main()
{
  int n=0;
  do {
    do_root_polish = n++;
    printf("\ndo_root_polish = %d\n", do_root_polish);
    test_roots(-2e+0, -1e+0, +1e+0, +2e+0);
    test_roots(-9e-1, -9e-1, +9e-1, +9e-1);
    test_roots(-9e-1, -9e-8, +9e-8, +9e-1);
    test_roots(-9e-1, -9e-10, +9e-10, +9e-1);
    test_roots(-9e-1, -9e-12, +9e-12, +9e-1);
    test_roots(-9e-1, -9e-14, +9e-14, +9e-1);

  } while (n != 2);
  return 0;
}
#endif /* __MAIN __ */
