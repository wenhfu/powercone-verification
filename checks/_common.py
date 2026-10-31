# -*- coding: utf-8 -*-
"""
Shared infrastructure for the symbolic verification of the upper-bound
algebra of

    "The optimal barrier parameter of the three-dimensional power cone"

(subsection 4.3, subsection 4.4 and Appendices A-C).

The checks cover the algebraic claims of that scope: identities, factorisations
and sign certificates.  The steps of the proof that are not algebraic (the
real-root continuation of Lemmas 2.4 and 4.3, the Sylvester step of Lemma 4.1,
and the passage from an identity to a sign statement) remain the paper's own
arguments.  Everything is computed in exact (rational / polynomial) arithmetic,
and every equality verdict is decided by exact normalisation (see ``is_zero``),
never by a heuristic simplification.

What this module provides
-------------------------
* Verdict tracking  : ``check``, ``silent``, ``banner``, ``reset_marks``,
                      ``report``.
* Exact decisions   : ``is_zero``, the zero test behind every equality verdict,
                      and ``clear_denominator``, an exact, remainder-checked
                      division.
* Text output       : ``show_coeffs``, ``show_factors``, ``num_factor_str``,
                      ``pat`` (plain text; polynomials are never rendered with
                      multi-line ``pprint``).
* Exact / interval  : ``primitive_int``, ``roots_in_unit``, ``bernstein_1d``
                      (used for the sign certificates).
* The symbols, the barrier ``f_e`` of (4.3) and all quantities derived from it,
  together with the explicit factors ``Pi_1``, ``Pi_2`` of Appendix A and the
  auxiliary polynomials of Appendix B.

Output policy
-------------
Only values computed by SymPy are printed.  The paper's formulas enter only
through the verdict, which names the equation being checked and reduces the
expression to a short handle, so the output has the form

    <value computed by SymPy>
       verified   <equation>   <expression>

in which every computed value can be read off, line by line, against the
corresponding displayed formula of the paper.

Usage
-----
This module is imported by the checks; it is not meant to be run on its own.  It
is where the shared symbols, the barrier, the auxiliary polynomials, and the
printing, Sturm and Bernstein helpers live.

Each ``appendix_*.py`` / ``sec*_*.py`` script opens with ``banner(<section>,
<formulas>)`` and closes with ``report(<section>)``, so every run states which
part of the paper and which equation numbers it covers.  A script therefore
runs standalone

    python checks/sec4_4_0.py

and the whole suite can be driven by ``python run_all.py``.
"""
import sys
from math import comb

import sympy as sp

# Make UTF-8 output work when the scripts are launched from a Windows console
# whose default code page is not UTF-8.
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


# ----------------------------------------------------------------------
# Verdict bookkeeping
# ----------------------------------------------------------------------
# Field widths of the OVERALL line: the bracketed section, then the verdict.
# They keep the parts of a run_all.py report aligned.
SUMMARY_WIDTH = 20
VERDICT_WIDTH = 4
ok_marks = []      # (name, passed) pairs collected by the running script
last_report = {}   # filled by report(); read by run_all.py
current_section = None    # part of the paper being checked (set by banner)
current_formulas = None   # equation numbers checked by the current script


def check(name, cond):
    """Record one verdict and print "verified   <equation>   <expression>".

    The line opens with the verdict, ``verified`` (or ``FAILED`` when the check
    does not hold), and then gives the equation that was checked together with a
    short handle for the expression.  The verdict comes first, so that the lines
    of a script line up under one another.
    """
    ok_marks.append((name, bool(cond)))
    print("   %s   %s" % (("verified" if cond else "FAILED").ljust(8), name))
    return bool(cond)


def silent(cond):
    """Same as check(), but silent: an ingredient of a compound verdict.

    Several elementary facts are often combined into a single printed line, so
    that the list of verdicts stays readable.
    """
    return bool(cond)


def banner(section, formulas=None):
    """Announce the part of the paper covered by the current script.

    Prints a header naming the section and the equation numbers, and stores both
    so that the closing summary carries the same numbers.
    """
    global current_section, current_formulas
    current_section = section
    current_formulas = formulas
    print("=" * 72)
    print("Checking %s%s" % (section, (": " + formulas) if formulas else ""))
    print("only SymPy-computed values are printed; compare them with the paper")
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
# (everything printed here is computed by SymPy; nothing is copied from the
# paper)
# ----------------------------------------------------------------------
def show_coeffs(name, expr, *vars):
    """Print a polynomial one coefficient per line.

    ``expr`` is expanded and read as a polynomial in ``vars``.  With a single
    variable the coefficients run from the leading term down to the constant;
    with several variables they are listed in descending monomial order.  A
    coefficient that still depends on other symbols is printed in factored form,
    while a pure number is printed as it stands, which is how the paper writes
    the two kinds of coefficients.
    """
    P = sp.Poly(sp.expand(expr), *vars)

    def fac(c):
        # Factor genuine polynomials; leave integer/rational constants alone.
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
    """Write a rational number as a product of prime powers, as the paper does.

    For instance 134150929655722049475182592 becomes '(2)^(44)*(3)^(27)', and
    -4096 becomes '(-1)*(2)^(12)'.
    """
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
    """Print a factorisation, one factor per line.

    SymPy keeps the rational constant as an ordinary numeric factor; it is
    expanded into prime powers here, so that it can be read off against the
    paper's forms, such as 2^44 * 3^27.
    """
    print("%s :" % name)
    for g in sp.Mul.make_args(sp.factor(expr)):
        if g.is_number:
            if g == 1:
                continue                # a unit factor carries no information
            print("   %s" % num_factor_str(g))
            continue
        bs, ex = g.as_base_exp()
        print("   (%s)^(%s)" % (bs, ex) if ex != 1 else "   (%s)" % bs)


def pat(v):
    """Sign pattern of a list of nonzero reals, e.g. (+,-,+,+,+,-,-)."""
    return "(" + ",".join("+" if x > 0 else "-" for x in v) + ")"


# ----------------------------------------------------------------------
# Exact / interval helpers
# ----------------------------------------------------------------------
def is_zero(expr):
    """True if and only if ``expr`` is identically zero: an exact decision.

    The expression is first put over a common denominator (``sp.cancel``) and
    the numerator is then expanded, so the verdict follows from polynomial
    normalisation.  ``sp.simplify`` is never used as a criterion: it is a
    heuristic that looks for a simpler expression and may fail to prove an
    identity that does hold.
    """
    num, den = sp.fraction(sp.cancel(expr))
    if sp.expand(den) == 0:
        raise ZeroDivisionError("is_zero(): vanishing denominator")
    return sp.expand(num) == 0


def clear_denominator(expr, *vars):
    """Exact quotient of a rational expression by its own denominator.

    ``expr`` is combined over a single denominator, and the numerator is then
    divided by that denominator in the polynomial ring ``Q[vars]``.  The
    remainder is checked to vanish, so the division proves that the denominator
    divides the numerator and that the result is a genuine polynomial.  This is
    how the normalising factors of (4.6) and (4.8), including the minus sign
    in the cubic relation, are cleared in exact arithmetic rather than by a
    heuristic simplification.
    """
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
    """Rescale ``g`` by a *positive* rational so that its coefficients become
    coprime integers.

    Multiplying a Sturm chain by a positive constant changes no sign, and hence
    no Sturm count; this is the normalisation used in Appendix C to make the
    chain readable.
    """
    Pg = sp.Poly(sp.expand(g), var)
    nums, dens = [], []
    for k in Pg.all_coeffs():
        num, den = sp.fraction(sp.cancel(sp.together(k)))
        nums.append(abs(num))
        dens.append(den)
    return sp.expand(sp.Rational(_lcm(dens), _gcd(nums))*g)


def roots_in_unit(p_expr, u_val):
    """Exact number of roots of ``p_expr(u_val, s)`` on the interval [0,1].

    ``sp.Poly.count_roots(0, 1)`` takes the endpoints into account as well, so a
    count of 0 certifies that the polynomial has no root anywhere in [0,1], a
    stronger statement than the count on (0,1) that subsection 4.4.3 quotes.
    """
    return sp.Poly(sp.expand(p_expr.subs(u, u_val)), s).count_roots(0, 1)


def bernstein_1d(p_expr, var, lo, hi, m):
    """Bernstein coefficients of degree ``m`` of ``p_expr`` on ``[lo,hi]``.

    The polynomial is first pulled back to the unit interval by
    ``var = lo + (hi-lo)*v``.  If every returned coefficient is positive, then
    the polynomial is positive on the whole interval.
    """
    v = sp.Symbol("v")
    Qv = sp.Poly(sp.expand(p_expr.subs(var, lo + (hi - lo)*v)), v)
    return [sp.expand(sum(Qv.coeff_monomial(v**k)*sp.Rational(comb(i, k), comb(m, k))
                         for k in range(i + 1))) for i in range(m + 1)]


# ----------------------------------------------------------------------
# Symbols
# ----------------------------------------------------------------------
# e is the barrier parameter kept free in subsection 4.3; t parametrises the
# slice (X,Z) = (1,t); p and r are the components of the tangent direction h;
# xi is the variable used for the univariate specialisations.
e, t, p, r, xi = sp.symbols("e t p r xi", positive=True)
# Primal variables of the power-cone barrier f_e.
X, Z = sp.symbols("X Z", positive=True)
# u = e^2 and s = t^2 are the variables in which Pi_1 and Pi_2 are polynomials.
u, s = sp.symbols("u s", positive=True)
# q = 1 - s is the variable of the base-point certificates in Appendix B.2.
q = sp.Symbol("q")

A = e**2 + 3      # shorthand for e^2 + 3, used throughout the paper
c = (1 - e)/2     # offset of the affine slice Z = c*t*p + r


# ----------------------------------------------------------------------
# The barrier f_e of (4.3) and every quantity derived from it.
# They are inputs to the checks and are never printed.
# ----------------------------------------------------------------------
f = -(3-e)/6*sp.log(X) - (9-e**2)/(6*A)*sp.log(1 - Z**2/X**(1-e))

D = lambda g, v: sp.diff(g, v)        # compact partial-derivative helper
fx, fz = D(f, X), D(f, Z)
fxx, fxz, fzz = D(fx, X), D(fx, Z), D(fz, Z)
fxxx, fxxz, fxzz, fzzz = D(fxx, X), D(fxx, Z), D(fxz, Z), D(fzz, Z)

# Tangent direction h = p * d/dX + (c*t*p + r) * d/dZ.
hX, hZ = p, c*t*p + r
# The contractions that define the quadratic form Q and the cubic form P, and
# the two self-concordance expressions used in subsection 4.3.
Df = fx*hX + fz*hZ                                    # <grad f, h>
D2f = fxx*hX**2 + 2*fxz*hX*hZ + fzz*hZ**2             # D^2 f[h,h]
D3f = fxxx*hX**3 + 3*fxxz*hX**2*hZ + 3*fxzz*hX*hZ**2 + fzzz*hZ**3
Ghh = (fxx - fx**2)*hX**2 + 2*(fxz - fx*fz)*hX*hZ + (fzz - fz**2)*hZ**2
Tfff = D3f - 6*D2f*Df + 4*Df**3

# Closed forms of (4.7) and (4.9), obtained from G[h,h] and T_f[h,h,h] on the
# slice (X,Z) = (1,t) after clearing the normalising factors, including the
# minus sign in the cubic relation.  The clearing is an exact,
# remainder-checked polynomial division (never
# sp.simplify), so Q_der and P_der are the genuine polynomials of (4.7) / (4.9).
Q_der = clear_denominator(Ghh.subs({X: 1, Z: t})*36*A**2*(1-t**2)**2/(9-e**2),
                          e, t, p, r)
P_der = clear_denominator(Tfff.subs({X: 1, Z: t})
                          * (-54*A**3*(1-t**2)**3/(e*(9-e**2))), e, t, p, r)

# Univariate specialisation p = 1, r = xi, and the polynomial
# d_{e,t}(xi) = Q(1,xi)^3 - P(1,xi)^2 whose discriminant drives the proof.
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
# The factors a_j, b_j are the compact forms printed in Appendix A, and their
# product is matched against Disc_xi d_{e,t} under subsection 4.4.1.
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
# Auxiliary factors of Appendix B
# (used by the Appendix B checks and by subsection 4.4.3)
# ----------------------------------------------------------------------
A6 = 1612*u**6 - 7822*u**5 + 35755*u**4 - 62928*u**3 + 62694*u**2 - 18954*u + 2187
A7 = (46301*u**7 - 430305*u**6 + 715218*u**5 + 2010910*u**4
      - 7240815*u**3 + 31031667*u**2 - 12772080*u + 1679616)
# R_7 and the three cubics h_0, h_1, h_2 of the base-point checks (Appendix B.2);
# q = 1 - s is the variable in which h_0, h_1, h_2 are written.
R7 = (38692*s**7 - 539840*s**6 + 3911572*s**5 - 14174132*s**4
      + 27167315*s**3 - 28286181*s**2 + 15176721*s - 3294172)
h1 = 50379285*q**3 + 1563956*q**2 - 766260*q + 49923
h2 = 40628031*q**3 - 747684*q**2 - 14364*q + 361
h0 = 143927*q**3 - 1624*q**2 - 560*q + 25
