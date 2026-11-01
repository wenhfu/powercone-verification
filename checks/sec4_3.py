# -*- coding: utf-8 -*-
"""
subsection 4.3: the quadratic form Q and the cubic form P.

Derives Q and P from f_e by differentiation using (4.6) and (4.8),
then checks the coefficients against (4.7) and (4.9).

Run:  python checks/sec4_3.py
"""
from _common import *

reset_marks()
banner('subsection 4.3', '(4.7), (4.9)')


def show_factor_coeffs(name, expr, *vars):
    """Print polynomial coefficients in factored form, including numeric coefficients."""
    P = sp.Poly(sp.expand(expr), *vars)
    print("%s   (total degree %d, coefficients factored):"
          % (name, P.total_degree()))
    for mono, co in sorted(P.terms(), key=lambda tc: tc[0], reverse=True):
        tag = " ".join("%s^%d" % (v, m) for v, m in zip(vars, mono))
        print("   [%s] = %s" % (tag, sp.factor(co)))


# Q and P as given in (4.7) and (4.9).
Q_pap = ((1-t**2)*A*(A-4*e**2*t**2)*p**2
         + 8*e*t*A*(1-t**2)*p*r
         + 4*(4*e**2*t**2+3*e**2+9)*r**2)
P_pap = (A**3*(1-t**2)**3*p**3
         + 12*e*t*A*(1-t**2)*(A+2*t**2*(e**2-3))*p**2*r
         - 12*A*(1-t**2)*(3*A+t**2*(5*e**2-9))*p*r**2
         - 16*e*t*(9*A+t**2*(5*e**2-9))*r**3)

# Q_der and P_der are computed from f_e in _common.py using (4.6) and (4.8).
show_factor_coeffs("Q(p,r) recomputed from G[h,h] of f_e", Q_der, p, r)
check("(4.7)   Q(p,r)", is_zero(Q_der - Q_pap))
print()
show_factor_coeffs("P(p,r) recomputed from T_f[h,h,h] of f_e", P_der, p, r)
check("(4.9)   P(p,r)", is_zero(P_der - P_pap))
print()

report("subsection 4.3")
