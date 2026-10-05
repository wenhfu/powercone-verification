# -*- coding: utf-8 -*-
"""
subsection 4.4.0 -- nonnegativity of the binary sextic L(u,s).

Paper: subsection 4.4.0
(eq:leading-before-L, eq:leading, eq:Lpower, eq:Lmatrix, eq:Lpositive).

Writing d_{e,t}(xi) = Q(1,xi)^3 - P(1,xi)^2 as a polynomial in xi, its leading
coefficient factors as

    [xi^6] d_{e,t} = 576 L(u,s),       u = e^2,  s = t^2,

with L the binary sextic of eq:leading.  The proof of positivity of d_{e,t}
uses the Bernstein representation of L on [0,1]^2.

What is verified here
---------------------
* The two coefficients of eq:Q / eq:P that are displayed again in this
  subsection, [xi^2] Q(1,xi) and [xi^3] P(1,xi).
* eq:leading-before-L : the closed form of the leading coefficient.
* eq:leading         : [xi^6] d_{e,t} = 576 L(u,s) with u = e^2, s = t^2.
* eq:Lpower          : the power-basis expansion of L(u,s).
* eq:Lmatrix         : the 4 x 4 matrix of Bernstein coefficients of L for
                       degree 3 in each variable.
* eq:Lpositive       : all those Bernstein coefficients are >= 0, which
                       certifies L(u,s) >= 0 on [0,1]^2.

Run:  python sec4_4_0.py
"""
from _common import *

reset_marks()
banner('subsection 4.4.0', 'eq:Q, eq:P, eq:leading-before-L, eq:leading, eq:Lpower, eq:Lmatrix, eq:Lpositive')

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
check("eq:Q   [xi^2] Q(1,xi)",
      sp.simplify(sp.expand(Q1.coeff(xi, 2) - 4*(4*e**2*t**2+3*e**2+9))) == 0)
print("[xi^3] P(1,xi) =", sp.factor(P1.coeff(xi, 3)))
check("eq:P   [xi^3] P(1,xi)",
      sp.simplify(sp.expand(P1.coeff(xi, 3) + 16*e*t*(9*A+t**2*(5*e**2-9)))) == 0)
print()

# eq:leading-before-L: the leading coefficient of d_{e,t} in the form of the
# paper, i.e. before the substitution u = e^2, s = t^2 is made.
check("eq:leading-before-L   [xi^6] d_{e,t}",
      sp.simplify(sp.expand(leading - lead_pap)) == 0)

# eq:leading: L is *defined* as [xi^6] d_{e,t} / 576.  Recompute it from the
# monomials of the leading coefficient (only even powers occur) and compare
# both with the compact closed form and with the power-basis expansion.
L_rec = sp.expand(sum(co*u**(i//2)*s**(j//2)
                      for (i, j), co in sp.Poly(leading, e, t).terms())/576)
check("eq:leading   [xi^6] d_{e,t} = 576 L(u,s)   with u=e^2, s=t^2",
      sp.simplify(sp.expand(leading - 576*L_expl.subs({u: e**2, s: t**2}))) == 0)
check("eq:Lpower   L(u,s) power-basis expansion", sp.expand(L_rec - L_pap) == 0)
print()

# eq:Lmatrix: Bernstein coefficients of L of bidegree (3,3).  For a polynomial
# written in the power basis, the Bernstein coefficient of bidegree (i,j) is
#     l_ij = sum_{k<=i, l<=j} a_kl * C(i,k)/C(3,k) * C(j,l)/C(3,l).
PolyL = sp.Poly(L_rec, u, s)
ell = sp.zeros(4, 4)
for i in range(4):
    for j in range(4):
        tot = 0
        for kk in range(i+1):
            for ll in range(j+1):
                akl = PolyL.coeff_monomial(u**kk*s**ll)
                tot += (akl*sp.Rational(comb(i, kk), comb(3, kk))
                        * sp.Rational(comb(j, ll), comb(3, ll)))
        ell[i, j] = sp.nsimplify(sp.expand(tot))
ell_pap = sp.Matrix([[81, 81, 81, 81],
                     [108, 84, 84, 96],
                     [144, 80, 64, sp.Rational(256, 3)],
                     [192, 64, 0, 0]])
print("Bernstein coefficients of L(u,s), m=3 in u and in s:")
sp.pprint(ell)
check("eq:Lmatrix   Bernstein coefficient matrix (l_ij)", (ell - ell_pap).is_zero_matrix)
check("eq:Lpositive   all l_ij >= 0", all(v >= 0 for v in ell))
print()

report("subsection 4.4.0")
