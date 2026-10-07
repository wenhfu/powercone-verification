# -*- coding: utf-8 -*-
"""
Appendix C: the exact Sturm chain of the base-point polynomial d_0.

Paper: Appendix C, together with (4.39) of
subsection 4.4.4.

Starting from S_0 = d_0 and S_1 = d_0', the Sturm chain is built by taking
negative Euclidean remainders repeatedly.  The paper prints a rescaled version
of the chain, in which every member is multiplied by a positive constant; such a
rescaling changes no sign and therefore no Sturm count.

What is verified here
---------------------
* The seven members S_0..S_6, rescaled to primitive integer polynomials, against
  the chain printed in Appendix C.
* The signs of the chain at -oo and +oo, hence V(-oo) = V(+oo) = 3, which by
  Sturm's theorem means that d_0 has no real root, so that d_0 is strictly
  positive (d_0(0) = 66199 > 0).
* The same conclusion, independently, from an exact root count.

Run:  python checks/appendix_c.py
"""
from _common import *

reset_marks()
banner('Appendix C', 'Sturm chain of d_0 of (4.39), Appendix C')

# SymPy returns the chain; primitive_int() rescales each member by a positive
# rational so that its coefficients are coprime integers, exactly as the paper
# does when printing S_0..S_6.
st = [g.as_expr() for g in sp.sturm(sp.Poly(d0, xi))]
S_int = [primitive_int(g, xi) for g in st]
print("Sturm chain of d_0, each member rescaled by a positive constant to a")
print("primitive integer polynomial (readable against Appendix C):")
for i in range(7):
    show_coeffs("   S_%d(xi)" % i, S_int[i], xi)
print()

# The same chain as printed in Appendix C.
S_pap = [
    17914752*xi**6 - 5322240*xi**5 - 3190320*xi**4 + 3596208*xi**3
    + 3949008*xi**2 + 94668*xi + 66199,
    8957376*xi**5 - 2217600*xi**4 - 1063440*xi**3 + 899052*xi**2
    + 658168*xi + 7889,
    2850379920*xi**4 - 3760659672*xi**3 - 6046460336*xi**2
    - 320056730*xi - 148801387,
    -330462921868609056*xi**3 - 240471376030379788*xi**2
    - 23796321104735200*xi - 5496922181504621,
    138655679724599952107248724*xi**2 - 3642291457427814879633064*xi
    + 3576940775423415451698955,
    3292010764698075971634963590183272*xi
    - 140414792326549871129269285851577,
    -1,
]
for i in range(7):
    check("Appendix C   S_%d(xi)" % i, sp.expand(S_int[i] - S_pap[i]) == 0)
print()


def _signs_at(m):
    """Signs of the chain members at m = -oo or m = +oo.

    At +oo the sign of a polynomial is the sign of its leading coefficient; at
    -oo it is that sign multiplied by (-1)^degree.
    """
    out = []
    for g in S_int:
        Pg = sp.Poly(g, xi)
        lc, dg = Pg.LC(), Pg.degree()
        out.append(sp.sign(lc) if m == sp.oo else sp.sign(lc)*(-1)**dg)
    return out


def _var(v):
    """Number of sign variations in a sign sequence, ignoring zeros."""
    sv = [x for x in v if x != 0]
    return sum(1 for i in range(len(sv)-1) if sv[i] != sv[i+1])


# Sturm's theorem: the number of distinct real roots in (a,b] is V(a) - V(b).
vm, vp = _signs_at(-sp.oo), _signs_at(sp.oo)
print("signs of the Sturm chain at -oo :", pat(vm))
print("signs of the Sturm chain at +oo :", pat(vp))
print("V(-oo) = %d ,  V(+oo) = %d" % (_var(vm), _var(vp)))
check("Appendix C   sign pattern at -oo", [int(x) for x in vm] == [1, -1, 1, 1, 1, -1, -1])
check("Appendix C   sign pattern at +oo", [int(x) for x in vp] == [1, 1, 1, -1, 1, 1, -1])
check("Appendix C   V(-oo) = V(+oo) = 3  =>  d_0 has no real root",
      _var(vm) == 3 and _var(vp) == 3)
check("Appendix C   d_0 has no real root (count_roots)",
      sp.Poly(d0, xi).count_roots(-sp.oo, sp.oo) == 0)
print()

report("Appendix C")
