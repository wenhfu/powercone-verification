# -*- coding: utf-8 -*-
"""
Shared symbols and helpers for Sections 4.3-4.4 and Appendices A-C.

Defines f_e in (4.3), its derivatives, and the polynomials of Appendices A-B.
Provides exact polynomial checks, Bernstein coefficients, root counts,
and helpers for printing and reporting.
"""
import sys
from math import comb

import sympy as sp

# Use UTF-8 output in Windows consoles with other default encodings.
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


# ----------------------------------------------------------------------
# Check reporting
# ----------------------------------------------------------------------
# Field widths aligning the section and verdict in OVERALL summaries.
SUMMARY_WIDTH = 20
VERDICT_WIDTH = 4
ok_marks = []      # (name, passed) pairs collected by the running script
last_report = {}   # filled by report(); read by run_all.py
current_section = None    # part of the paper being checked (set by banner)
current_formulas = None   # equation numbers checked by the current script


def check(name, cond):
    """Record the check and print its name with ``verified`` or ``FAILED``."""
    ok_marks.append((name, bool(cond)))
    print("   %s   %s" % (("verified" if cond else "FAILED").ljust(8), name))
    return bool(cond)


def silent(cond):
    """Convert a condition to bool without recording or printing a check."""
    return bool(cond)


def banner(section, formulas=None):
    """Print the section and equation numbers, saving them for the summary."""
    global current_section, current_formulas
    current_section = section
    current_formulas = formulas
    print("=" * 72)
    print("Checking %s%s" % (section, (": " + formulas) if formulas else ""))
    print("=" * 72)
    print()


def reset_marks():
    """Discard the verdicts recorded so far (called at the top of each script)."""
    del ok_marks[:]


def report(section):
    """Print the summary, reset the verdicts, and exit with status 1 on failure."""
    global current_section, current_formulas
    bad = [n for n, ok in ok_marks if not ok]
    last_report.update(section=section, n=len(ok_marks), bad=bad)
    if current_formulas:
        tag = current_formulas
        if not tag.startswith("("):
            tag = "(" + tag + ")"
        tag = "   " + tag
    else:
        tag = ""
    print("=" * 72)
    print("OVERALL %s%s%s"
          % (("[%s]" % section).ljust(SUMMARY_WIDTH),
             ("PASS" if not bad else "FAIL").ljust(VERDICT_WIDTH), tag))
    for n in bad:
        print("   FAILED: %s" % n)
    reset_marks()
    current_section = None
    current_formulas = None
    if bad:
        raise SystemExit(1)


# ----------------------------------------------------------------------
# Printing helpers
# ----------------------------------------------------------------------
def show_coeffs(name, expr, *vars):
    """Print coefficients in descending degree or monomial order.

    Factor symbolic coefficients and leave numeric coefficients unchanged.
    """
    P = sp.Poly(sp.expand(expr), *vars)

    def fac(c):
        # Factor symbolic coefficients; leave numeric constants unchanged.
        return c if c.is_number else sp.factor(c)

    if len(vars) == 1:
        v = vars[0]
        print("%s   (degree %d in %s):" % (name, P.degree(v), v))
        for i in range(P.degree(v), -1, -1):
            print("   [%s^%d] = %s" % (v, i, fac(P.coeff_monomial(v**i))))
    else:
        print("%s   (total degree %d):" % (name, P.total_degree()))
        for mono, co in sorted(P.terms(), key=lambda tc: tc[0], reverse=True):
            tag = " ".join("%s^%d" % (v, m) for v, m in zip(vars, mono))
            print("   [%s] = %s" % (tag, fac(co)))


def num_factor_str(n):
    """Format a nonzero rational number as a product of prime powers."""
    n = sp.Rational(n)
    parts = []
    if n < 0:
        parts.append("(-1)")
        n = -n
    num, den = sp.fraction(n)
    for prime, exp in sorted(sp.factorint(num).items()):
        parts.append("(%d)^(%d)" % (prime, exp) if exp != 1 else "(%d)" % prime)
    for prime, exp in sorted(sp.factorint(den).items()):
        parts.append("(%d)^(-%d)" % (prime, exp))
    return "*".join(parts)


def show_factors(name, expr):
    """Print one factor per line, writing numeric constants as prime powers."""
    print("%s :" % name)
    for g in sp.Mul.make_args(sp.factor(expr)):
        if g.is_number:
            if g == 1:
                continue
            print("   %s" % num_factor_str(g))
            continue
        bs, ex = g.as_base_exp()
        print("   (%s)^(%s)" % (bs, ex) if ex != 1 else "   (%s)" % bs)


def pat(v):
    """Sign pattern of a list of nonzero reals, e.g. (+,-,+,+,+,-,-)."""
    return "(" + ",".join("+" if x > 0 else "-" for x in v) + ")"


# ----------------------------------------------------------------------
# Exact polynomial calculations
# ----------------------------------------------------------------------
def is_zero(expr):
    """Test whether the numerator vanishes after rational cancellation."""
    num, den = sp.fraction(sp.cancel(expr))
    if sp.expand(den) == 0:
        raise ZeroDivisionError("is_zero(): vanishing denominator")
    return sp.expand(num) == 0


def clear_denominator(expr, *vars):
    """Divide the numerator by the denominator and reject a nonzero remainder."""
    num, den = sp.fraction(sp.together(expr))
    q, rem = sp.div(sp.Poly(sp.expand(num), *vars),
                    sp.Poly(sp.expand(den), *vars))
    if not rem.is_zero:
        raise ArithmeticError("the denominator does not divide the numerator")
    return sp.expand(q.as_expr())


def _lcm(seq):
    seq = list(seq)
    return seq[0] if len(seq) == 1 else sp.ilcm(*seq)


def _gcd(seq):
    seq = list(seq)
    return seq[0] if len(seq) == 1 else sp.igcd(*seq)


def primitive_int(g, var):
    """Rescale a nonzero rational polynomial to coprime integer coefficients.

    Use a positive scale to preserve signs in the Sturm chain of Appendix C.
    """
    Pg = sp.Poly(sp.expand(g), var)
    nums, dens = [], []
    for k in Pg.all_coeffs():
        num, den = sp.fraction(sp.cancel(sp.together(k)))
        nums.append(abs(num))
        dens.append(den)
    return sp.expand(sp.Rational(_lcm(dens), _gcd(nums))*g)


def roots_in_unit(p_expr, u_val):
    """Count real roots of ``p_expr(u_val, s)`` in [0,1], including both endpoints."""
    return sp.Poly(sp.expand(p_expr.subs(u, u_val)), s).count_roots(0, 1)


def bernstein_1d(p_expr, var, lo, hi, m):
    """Return degree-m Bernstein coefficients on [lo,hi].

    Use var = lo + (hi-lo)*v. The degree m must be a nonnegative
    integer at least as large as the degree after substitution.
    """
    if not isinstance(m, (int, sp.Integer)) or m < 0:
        raise ValueError("the Bernstein degree must be a nonnegative integer")
    v = sp.Symbol("v")
    Qv = sp.Poly(sp.expand(p_expr.subs(var, lo + (hi - lo)*v)), v)
    if Qv.degree() > m:
        raise ValueError("the polynomial degree exceeds the Bernstein degree")
    return [sp.expand(sum(Qv.coeff_monomial(v**k)*sp.Rational(comb(i, k), comb(m, k))
                         for k in range(i + 1))) for i in range(m + 1)]


# ----------------------------------------------------------------------
# Symbols
# ----------------------------------------------------------------------
# Parameters e and t at the normalized point (X,Z) = (1,t).
# p, r are direction coordinates; xi is the univariate specialisation variable.
e, t, p, r, xi = sp.symbols("e t p r xi", positive=True)
# Normalized coordinates of f_e in (4.3).
X, Z = sp.symbols("X Z", positive=True)
# u = e^2 and s = t^2 are the variables in which Pi_1 and Pi_2 are polynomials.
u, s = sp.symbols("u s", positive=True)
# q = 1 - s in the base-point checks of Appendix B.2.
q = sp.Symbol("q")

A = e**2 + 3
c = (1 - e)/2     # h_Z = c*t*p + r in the direction coordinates


# ----------------------------------------------------------------------
# The barrier f_e of (4.3) and its derived quantities.
# ----------------------------------------------------------------------
f = -(3-e)/6*sp.log(X) - (9-e**2)/(6*A)*sp.log(1 - Z**2/X**(1-e))

D = lambda g, v: sp.diff(g, v)
fx, fz = D(f, X), D(f, Z)
fxx, fxz, fzz = D(fx, X), D(fx, Z), D(fz, Z)
fxxx, fxxz, fxzz, fzzz = D(fxx, X), D(fxx, Z), D(fxz, Z), D(fzz, Z)

# Tangent direction h = p * d/dX + (c*t*p + r) * d/dZ.
hX, hZ = p, c*t*p + r
# Directional derivatives and the forms G[h,h] and T_f[h,h,h].
Df = fx*hX + fz*hZ                                    # <grad f, h>
D2f = fxx*hX**2 + 2*fxz*hX*hZ + fzz*hZ**2             # D^2 f[h,h]
D3f = fxxx*hX**3 + 3*fxxz*hX**2*hZ + 3*fxzz*hX*hZ**2 + fzzz*hZ**3
Ghh = (fxx - fx**2)*hX**2 + 2*(fxz - fx*fz)*hX*hZ + (fzz - fz**2)*hZ**2
Tfff = D3f - 6*D2f*Df + 4*Df**3

# Recover Q and P using (4.6) and (4.8) at the normalized point (X,Z) = (1,t).
Q_der = clear_denominator(Ghh.subs({X: 1, Z: t})*36*A**2*(1-t**2)**2/(9-e**2),
                          e, t, p, r)
P_der = clear_denominator(Tfff.subs({X: 1, Z: t})
                          * (-54*A**3*(1-t**2)**3/(e*(9-e**2))), e, t, p, r)

# Set p = 1, r = xi in d_{e,t}(xi) = Q(1,xi)^3 - P(1,xi)^2.
Q1 = sp.expand(Q_der.subs({p: 1, r: xi}))
P1 = sp.expand(P_der.subs({p: 1, r: xi}))
d = sp.expand(Q1**3 - P1**2)
leading = sp.expand(d.coeff(xi, 6))                   # [xi^6] d_{e,t} = 576 L(u,s)

# Base-point polynomial d_0 of (4.39), cleared of denominators.
d0 = sp.expand(sp.Rational(729, 2)*d.subs({e: sp.Rational(1, 3), t: sp.Rational(1, 2)}))


# ----------------------------------------------------------------------
# Appendix A
#   Pi_1(u,s) = sum_{j=0..6} a_j(u) s^j
#   Pi_2(u,s) = sum_{j=0..9} b_j(u) s^j
# Coefficients a_j and b_j of Appendix A.
# The discriminant factorisation using Pi_1 and Pi_2 is checked in sec4_4_1.py.
# ----------------------------------------------------------------------
a = [-4*(u+3)**6,
     12*(u+3)**6,
     12*(u+3)**4*(3*u**2-10*u-9),
     -(u+3)**3*(23*u**3+294*u**2-465*u-108),
     12*u*(u-1)*(u+3)**2*(2*u**2-u+63),
     -12*u*(u-3)*(u-1)**2*(u+3)*(u+9),
     16*u**2*(u-1)**2*(u**2+2*u+9)]
b = [4*(u+3)**9,
     -(u+3)**8*(u+72),
     (u+3)**7*(u**2-45*u+576),
     (u+3)**6*(13*u**3-42*u**2+777*u-2700),
     -4*(u+3)**5*(5*u**4+27*u**3-177*u**2+1269*u-2052),
     4*(u+3)**4*(6*u**5+21*u**4+117*u**3-1237*u**2+4473*u-4212),
     4*(u+3)**3*(6*u**6-81*u**5+111*u**4-662*u**3+4392*u**2-9369*u+5859),
     4*(u-1)*(u+3)**2*(u**6-50*u**5+340*u**4-492*u**3+2049*u**2-6426*u+5346),
     -24*(u-1)**2*(u+3)*(u**5-22*u**4+82*u**3-108*u**2+405*u-486),
     4*(u-9)*(u-1)**3*(9*u**3-19*u**2+27*u-81)]

Pi1 = sp.expand(sum(a[j]*s**j for j in range(7)))
Pi2 = sp.expand(sum(b[j]*s**j for j in range(10)))


# ----------------------------------------------------------------------
# Auxiliary polynomials of Appendix B
# (used by the Appendix B checks and by Section 4.4.3)
# ----------------------------------------------------------------------
A6 = 1612*u**6 - 7822*u**5 + 35755*u**4 - 62928*u**3 + 62694*u**2 - 18954*u + 2187
A7 = (46301*u**7 - 430305*u**6 + 715218*u**5 + 2010910*u**4
      - 7240815*u**3 + 31031667*u**2 - 12772080*u + 1679616)
# Appendix B.2: R_7(s) and h_0(q), h_1(q), h_2(q), with q = 1-s.
R7 = (38692*s**7 - 539840*s**6 + 3911572*s**5 - 14174132*s**4
      + 27167315*s**3 - 28286181*s**2 + 15176721*s - 3294172)
h1 = 50379285*q**3 + 1563956*q**2 - 766260*q + 49923
h2 = 40628031*q**3 - 747684*q**2 - 14364*q + 361
h0 = 143927*q**3 - 1624*q**2 - 560*q + 25
