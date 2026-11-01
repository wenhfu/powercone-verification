# -*- coding: utf-8 -*-
"""
subsection 4.4.3: the sign of Pi_2.

Checks the boundary values (4.32), discriminant factorisation (4.33),
and nonvanishing (4.34), together with exact root counts at
u = 1/4, 1/2, 3/4 used in the continuation argument for (4.37).
The positivity of A_6 and A_7 is checked by Bernstein coefficients
as in Appendix B.1.

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

# Appendix B.1: Bernstein coefficients of A_6 and A_7 on each half of [0,1].
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

# (4.34): all factors are nonzero for 0<u<1 except (2u-1)^3,
# which vanishes at u = 1/2. Lemma 4.3 applies on the two connected intervals.
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

# Base-point checks for the continuation argument in (4.37).
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
