# -*- coding: utf-8 -*-
"""
subsection 4.4.0: nonnegativity of the binary sextic L(u,s).

Paper: subsection 4.4.0, equations (4.13), (4.15), (4.17), (4.18), (4.20), (4.21).

Written as a polynomial in xi, d_{e,t}(xi) = Q(1,xi)^3 - P(1,xi)^2 has leading
coefficient

    [xi^6] d_{e,t} = 576 L(u,s),       u = e^2,  s = t^2,

with L the binary sextic of (4.17).  The proof that d_{e,t} is positive uses the
Bernstein representation of L on [0,1]^2.  The degenerate case t = 0 is treated
separately by the paper (Lemma 4.2, (4.13)) and is checked here first.

What is verified here
---------------------
* The two coefficients of (4.7) / (4.9) that are displayed again in this
  subsection, [xi^2] Q(1,xi) and [xi^3] P(1,xi).
* (4.13) : for t = 0 the identity
               Delta(p,r) = 108 r^2 A^3 (A p^2 - 4 r^2)^2,
           which is >= 0 because A = e^2 + 3 > 0 (Lemma 4.2).
* (4.15) : the closed form of the leading coefficient.
* (4.17) : [xi^6] d_{e,t} = 576 L(u,s) with u = e^2, s = t^2.
* (4.18) : the power-basis expansion of L(u,s).
* (4.20) : the 4 x 4 matrix of Bernstein coefficients of L, of degree 3 in each
           variable.
* (4.21) : L(u,s) > 0 on (0,1)^2.  Subtracting 81(1-u)^3 (the first row of the
           Bernstein matrix, by the partition of unity) leaves a polynomial
           whose Bernstein coefficients are the remaining rows, all >= 0.

Run:  python checks/sec4_4_0.py
"""
from _common import *

reset_marks()
banner('subsection 4.4.0', '(4.7), (4.9), (4.13), (4.15), (4.17), (4.18), (4.20), (4.21)')


def bern_coeffs(expr, du, ds):
    """Bernstein coefficients of bidegree (du,ds) of ``expr`` on [0,1]^2.

    For a polynomial written in the power basis, the Bernstein coefficient of
    bidegree (i,j) is

        l_ij = sum_{k<=i, l<=j} a_kl * C(i,k)/C(du,k) * C(j,l)/C(ds,l).
    """
    P = sp.Poly(sp.expand(expr), u, s)
    M = sp.zeros(du + 1, ds + 1)
    for i in range(du + 1):
        for j in range(ds + 1):
            M[i, j] = sp.expand(sum(
                P.coeff_monomial(u**k*s**l)*sp.Rational(comb(i, k), comb(du, k))
                * sp.Rational(comb(j, l), comb(ds, l))
                for k in range(i + 1) for l in range(j + 1)))
    return M

# Closed forms transcribed from the paper: the leading coefficient as it is
# written before L is introduced, the same quantity rewritten through L, and
# the power-basis expansion of L(u,s).
lead_pap = 64*((4*e**2*t**2+3*e**2+9)**3
               - 4*e**2*t**2*(9*(e**2+3)+t**2*(5*e**2-9))**2)
L_expl = sp.expand(sp.Rational(1, 9)*((4*u*s+3*u+9)**3
                                      - 4*u*s*(9*(u+3)+s*(5*u-9))**2))
L_pap = (81 + 81*u + 27*u**2 + 3*u**3
         - 216*u*s - 144*u**2*s - 24*u**3*s
         + 216*u*s**2 - 24*u**3*s**2
         - 36*u*s**3 + 40*u**2*s**3 - 4*u**3*s**3)

# The two coefficients of Q and P that the paper quotes again in this
# subsection; the printed value is what SymPy computes from f_e.
print("[xi^2] Q(1,xi) =", sp.factor(Q1.coeff(xi, 2)))
check("(4.7)   [xi^2] Q(1,xi)",
      is_zero(sp.expand(Q1.coeff(xi, 2) - 4*(4*e**2*t**2+3*e**2+9))))
print("[xi^3] P(1,xi) =", sp.factor(P1.coeff(xi, 3)))
check("(4.9)   [xi^3] P(1,xi)",
      is_zero(sp.expand(P1.coeff(xi, 3) + 16*e*t*(9*A+t**2*(5*e**2-9)))))
print()

# Lemma 4.2 / (4.13): the paper treats t = 0 separately.  For t = 0,
# Delta(p,r) = Q(p,r)^3 - P(p,r)^2 factors with a squared cofactor, hence is
# nonnegative as soon as A = e^2 + 3 > 0.
Delta0 = sp.expand(sp.expand(Q_der**3 - P_der**2).subs(t, 0))
print("Delta(p,r) at t = 0 =", sp.factor(Delta0))
check("(4.13)   Delta(p,r) for t = 0",
      is_zero(Delta0 - 108*r**2*A**3*(A*p**2 - 4*r**2)**2))
check("(4.13)   Delta(p,r) >= 0 for t = 0   (Delta = 108 r^2 A^3 (square))",
      is_zero(sp.cancel(Delta0/(108*r**2*A**3)) - (A*p**2 - 4*r**2)**2))
print()

# (4.15): the leading coefficient of d_{e,t} in the paper's form, that is,
# before the substitution u = e^2, s = t^2 is made.
check("(4.15)   [xi^6] d_{e,t}",
      is_zero(sp.expand(leading - lead_pap)))

# (4.17): L is *defined* as [xi^6] d_{e,t} / 576.  Recompute it from the
# monomials of the leading coefficient (only even powers occur), and compare it
# both with the compact closed form and with the power-basis expansion.
L_rec = sp.expand(sum(co*u**(i//2)*s**(j//2)
                      for (i, j), co in sp.Poly(leading, e, t).terms())/576)
check("(4.17)   [xi^6] d_{e,t} = 576 L(u,s)   with u=e^2, s=t^2",
      is_zero(sp.expand(leading - 576*L_expl.subs({u: e**2, s: t**2}))))
check("(4.18)   L(u,s) power-basis expansion", sp.expand(L_rec - L_pap) == 0)
print()

# (4.20): Bernstein coefficients of L of bidegree (3,3).
ell = bern_coeffs(L_rec, 3, 3)
ell_pap = sp.Matrix([[81, 81, 81, 81],
                     [108, 84, 84, 96],
                     [144, 80, 64, sp.Rational(256, 3)],
                     [192, 64, 0, 0]])
print("Bernstein coefficients of L(u,s), m=3 in u and in s:")
sp.pprint(ell)
check("(4.20)   Bernstein coefficient matrix (l_ij)", (ell - ell_pap).is_zero_matrix)
# (4.21): L(u,s) > 0 for 0<u,s<1.  The paper's argument uses only the
# first row of the matrix and the partition of unity,
#     81 B_0^(3)(u) sum_j B_j^(3)(s) = 81(1-u)^3,
# so it is enough to look at the Bernstein coefficients of L - 81(1-u)^3,
# which are the last three rows above.
ell_rest = bern_coeffs(sp.expand(L_rec - 81*(1-u)**3), 3, 3)
print("Bernstein coefficients of L(u,s) - 81(1-u)^3:")
sp.pprint(ell_rest)
check("(4.21)   L - 81(1-u)^3 is the sum of the Bernstein terms with i >= 1",
      ell_rest == sp.Matrix([[0 if i == 0 else ell[i, j] for j in range(4)]
                             for i in range(4)]))
check("(4.21)   L(u,s) > 0 for 0<u,s<1   (all those coefficients >= 0)",
      all(v >= 0 for v in ell_rest))
print()

report("subsection 4.4.0")
