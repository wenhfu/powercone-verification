# -*- coding: utf-8 -*-
"""
subsection 4.3 -- the quadratic form Q and the cubic form P.

Paper: subsection 4.3 (eq:Q and eq:P).

For the power-cone barrier f_e of eq:fdef the paper derives closed forms for
the quadratic and the cubic form

    Q(p,r) = G[h,h]        (self-concordance form),        (eq:Q)
    P(p,r) = T_f[h,h,h]    (third-order form),             (eq:P)

on the slice (X,Z) = (1,t), where h = p*d/dX + (c*t*p + r)*d/dZ.

What is verified here
---------------------
* The coefficients of Q(p,r) and P(p,r), recomputed from f_e by symbolic
  differentiation, coincide with the closed forms of eq:Q / eq:P.

Every coefficient is printed in factored form so that it can be compared with
the compact factors used in the paper.

Run:  python checks/sec4_3.py
"""
from _common import *

reset_marks()
banner('subsection 4.3', 'eq:Q, eq:P')


def show_factor_coeffs(name, expr, *vars):
    """Print the coefficients of ``expr``, one line per monomial, factored.

    Unlike ``show_coeffs`` this factors coefficients that are pure numbers as
    well, because subsection 4.3 compares its coefficients with the product
    form in which eq:Q / eq:P are written.
    """
    P = sp.Poly(sp.expand(expr), *vars)
    print("%s   (total degree %d, coefficients factored):"
          % (name, P.total_degree()))
    for mono, co in sorted(P.terms(), key=lambda tc: tc[0], reverse=True):
        tag = " ".join("%s^%d" % (v, m) for v, m in zip(vars, mono))
        print("   [%s] = %s" % (tag, sp.factor(co)))


# Closed forms of eq:Q and eq:P, transcribed as they appear in the paper.
Q_pap = ((1-t**2)*A*(A-4*e**2*t**2)*p**2
         + 8*e*t*A*(1-t**2)*p*r
         + 4*(4*e**2*t**2+3*e**2+9)*r**2)
P_pap = (A**3*(1-t**2)**3*p**3
         + 12*e*t*A*(1-t**2)*(A+2*t**2*(e**2-3))*p**2*r
         - 12*A*(1-t**2)*(3*A+t**2*(5*e**2-9))*p*r**2
         - 16*e*t*(9*A+t**2*(5*e**2-9))*r**3)

# Q_der and P_der come from the contractions
#     G[h,h]    = (f_xx - f_x^2) h_X^2 + 2 (f_xz - f_x f_z) h_X h_Z
#                 + (f_zz - f_z^2) h_Z^2,
#     T_f[h,h,h] = D^3 f[h,h,h] - 6 D^2 f[h,h] <grad f,h> + 4 <grad f,h>^3,
# evaluated on the slice (X,Z) = (1,t) and cleared of positive factors.
show_factor_coeffs("Q(p,r) recomputed from G[h,h] of f_e", Q_der, p, r)
check("eq:Q   Q(p,r)", is_zero(Q_der - Q_pap))
print()
show_factor_coeffs("P(p,r) recomputed from T_f[h,h,h] of f_e", P_der, p, r)
check("eq:P   P(p,r)", is_zero(P_der - P_pap))
print()

report("subsection 4.3")
