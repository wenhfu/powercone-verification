# -*- coding: utf-8 -*-
"""
subsection 4.4.3: the sign of Pi_2.

Paper: subsection 4.4.3, with the auxiliary computations of Appendix B,
equations (4.32), (4.33), (4.34), (4.37).

The degree-9 factor Pi_2(u,s) is positive on (0,1)^2 by boundary values,
the discriminant in s, and three exact base-point checks, together with the
continuation argument of Lemma 4.3.

What is verified here
---------------------
* (4.32) : the two boundary values Pi_2(u,0) and Pi_2(u,1).
* (4.33) : the factorisation of Disc_s Pi_2, whose factors are A_6 and A_7
           (checked to be positive in Appendix B.1) together with simple
           positive factors on (0,1).
* (4.34) : Disc_s Pi_2(u,.) != 0 for 0<u<1/2 and for 1/2<u<1, which is the
           hypothesis of Lemma 4.3 on the two halves; only (2u-1)^3 vanishes,
           at u = 1/2.
* (4.37) : Pi_2(u,s) > 0 on (0,1)^2, using exact root
           counts of Pi_2(u0,.) on (0,1) for u0 = 1/4, 1/2, 3/4, and the
           Bernstein certificates for A_6 and A_7 above.

Run:  python checks/sec4_4_3.py
"""
from _common import *

reset_marks()
banner('subsection 4.4.3', '(4.32)-(4.34), (4.37)')

# (4.32): the two boundary values of Pi_2, both printed factored.
b0, b1 = Pi2.subs(s, 0), Pi2.subs(s, 1)
print("Pi_2(u,0) =", sp.factor(b0))
print("Pi_2(u,1) =", sp.factor(b1))
check("(4.32)   Pi_2(u,0)", is_zero(b0 - 4*(u+3)**9))
check("(4.32)   Pi_2(u,1)", is_zero(b1 - u*(1-u)**6*(7*u+9)**2))
print()

# Bernstein coefficients of the two auxiliary factors A_6 and A_7 on the two
# halves of [0,1], chosen to obtain positive coefficients on each subinterval.
# The factor 2u-1 in Disc_s Pi_2 separately determines the continuation intervals.
cA6a = bernstein_1d(A6, u, 0, sp.Rational(1, 2), 6)
cA6b = bernstein_1d(A6, u, sp.Rational(1, 2), 1, 6)
cA7a = bernstein_1d(A7, u, 0, sp.Rational(1, 2), 7)
cA7b = bernstein_1d(A7, u, sp.Rational(1, 2), 1, 7)

# (4.33): Disc_s Pi_2 = -2^12 u^8 (u-1)^24 (u+3)^72 (u+7)^6
#                             (2u-1)^3 (7u+9)^3 A_6(u)^3 A_7(u).
Disc2 = sp.discriminant(sp.Poly(Pi2, s))
Disc2_ref = (-2**12*u**8*(u-1)**24*(u+3)**72*(u+7)**6
             * (2*u-1)**3*(7*u+9)**3*A6**3*A7)
show_factors("factor(Disc_s Pi_2)", Disc2)
check("(4.33)   Disc_s Pi_2 factorisation", is_zero(sp.expand(Disc2 - Disc2_ref)))
print()

# (4.34): A_6, A_7 > 0 and the other bases have no roots in (0,1).
# Only (2u-1)^3 vanishes there, at u = 1/2; Lemma 4.3 applies on each side.
half = sp.Rational(1, 2)
bases = [u, u - 1, u + 3, u + 7, 7*u + 9]
crossing = [(b, rt) for b in bases
            for rt in sp.roots(sp.Poly(b, u).as_expr(), u)
            if 0 < rt < half or half < rt < 1]
print("roots in (0,1) of the bases u, u-1, u+3, u+7, 7u+9 :", crossing)
check("(4.34)   Disc_s Pi_2(u,.) != 0 for 0<u<1/2 and 1/2<u<1",
      not crossing and sp.roots(2*u - 1, u) == {half: 1}
      and all(co > 0 for co in cA6a + cA6b)
      and all(co > 0 for co in cA7a + cA7b))
print()

# (4.37): combine positive boundary values, Bernstein certificates for A_6, A_7,
# the discriminant identity, and exact root counts with a positive value at
# s = 0 for the base points u0 = 1/4, 1/2, 3/4.
pi2_ok = []
pi2_ok.append(silent(is_zero(b0 - 4*(u+3)**9)))
pi2_ok.append(silent(is_zero(b1 - u*(1-u)**6*(7*u+9)**2)))
pi2_ok.append(silent(all(co > 0 for co in cA6a + cA6b)))
pi2_ok.append(silent(all(co > 0 for co in cA7a + cA7b)))
pi2_ok.append(silent(is_zero(sp.expand(Disc2 - Disc2_ref))))
for u0 in [sp.Rational(1, 4), sp.Rational(3, 4), sp.Rational(1, 2)]:
    pi2_ok.append(silent(roots_in_unit(Pi2, u0) == 0
                         and Pi2.subs({u: u0, s: 0}) > 0))
check("(4.37)   Pi_2(u,s) > 0 on (0,1)^2   (exact Sturm + Bernstein certificates)",
      all(pi2_ok))
print()

report("subsection 4.4.3")
