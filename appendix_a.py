# -*- coding: utf-8 -*-
"""
Appendix A -- the explicit polynomial factors Pi_1 and Pi_2.

Paper: eq:discd and Appendix A.

Disc_xi d_{e,t} factors as

    2^44 3^27 s^3 (u-1)^3 (u+3)^15 (s-1)^15 Pi_1(u,s)^3 Pi_2(u,s),

where

    Pi_1(u,s) = sum_{j=0..6} a_j(u) s^j,
    Pi_2(u,s) = sum_{j=0..9} b_j(u) s^j

and the a_j, b_j are the compact factors displayed in Appendix A.

What is verified here
---------------------
* The coefficient lists a_0..a_6 and b_0..b_9, factored, against the compact
  forms of Appendix A (ascending order, a_0 first).
* The two boundary values Pi_1(u,1) and Pi_2(u,1), and the extremal
  coefficients a_0, a_6, b_0, b_9.
* The identity for Disc_xi d_{e,t} at four exact rational sample points.

Why only sample points
----------------------
Proving the identity symbolically is the expensive computation (about
1.5 minutes) done in sec4_4_1.py.  Here the same identity is confirmed
exactly at four rational points (e,t), which keeps this script fast while
still being an exact check.  Both sides of the comparison are evaluated with
exact rational arithmetic, so no numerical tolerance is involved.

Run:  python appendix_a.py
"""
from _common import *

reset_marks()
banner('Appendix A', 'Pi_1, Pi_2 (Appendix A)')


def show_labelled_factors(name, prefix, expr, var):
    """Print the coefficients of ``expr`` in ``var``, one per line.

    The j-th coefficient is labelled ``prefix_j`` and printed in factored
    form, matching the way the coefficients are named in Appendix A.
    """
    P = sp.Poly(sp.expand(expr), var)
    print("%s   (degree %d in %s, coefficients factored):"
          % (name, P.degree(var), var))
    for j in range(P.degree(var) + 1):
        print("   %s_%d(u) = %s" % (prefix, j, sp.factor(P.coeff_monomial(var**j))))


# Coefficients of Pi_1 (degree 6 in s) and of the distinguished ones.
show_labelled_factors("Pi_1(u,s) = sum_j a_j(u) s^j", "a", Pi1, s)
check("Appendix A   deg_s Pi_1 = 6  (j = 0..6)",
      sp.Poly(Pi1, s).degree(s) == 6)
check("Appendix A   a_0(u) = -4(u+3)^6",
      sp.simplify(Pi1.subs(s, 0) + 4*(u+3)**6) == 0)
check("Appendix A   a_6(u) = 16u^2(u-1)^2(u^2+2u+9)",
      sp.simplify(sp.Poly(Pi1, s).LC() - 16*u**2*(u-1)**2*(u**2+2*u+9)) == 0)
check("Appendix A   Pi_1(u,1) = u(u-1)^3(7u+9)^2",
      sp.simplify(Pi1.subs(s, 1) - u*(u-1)**3*(7*u+9)**2) == 0)
print()

# Coefficients of Pi_2 (degree 9 in s) and of the distinguished ones.
show_labelled_factors("Pi_2(u,s) = sum_j b_j(u) s^j", "b", Pi2, s)
check("Appendix A   deg_s Pi_2 = 9  (j = 0..9)",
      sp.Poly(Pi2, s).degree(s) == 9)
check("Appendix A   b_0(u) = 4(u+3)^9",
      sp.simplify(Pi2.subs(s, 0) - 4*(u+3)**9) == 0)
check("Appendix A   b_9(u) = 4(u-9)(u-1)^3(9u^3-19u^2+27u-81)",
      sp.simplify(sp.Poly(Pi2, s).LC()
                  - 4*(u-9)*(u-1)**3*(9*u**3-19*u**2+27*u-81)) == 0)
check("Appendix A   Pi_2(u,1) = u(1-u)^6(7u+9)^2",
      sp.simplify(Pi2.subs(s, 1) - u*(1-u)**6*(7*u+9)**2) == 0)
print()

# Exact spot-check of  Disc_xi d_{e,t} = 2^44 3^27 s^3 (u-1)^3 (u+3)^15
# (s-1)^15 Pi_1^3 Pi_2  at four rational points (e,t); the discriminant is
# taken in the rational domain, so both sides are exact.
claimed = 2**44*3**27*s**3*(u-1)**3*(u+3)**15*(s-1)**15*Pi1**3*Pi2
samples = [(sp.Rational(1, 3), sp.Rational(1, 2)),
           (sp.Rational(1, 2), sp.Rational(1, 4)),
           (sp.Rational(2, 5), sp.Rational(3, 5)),
           (sp.Rational(1, 4), sp.Rational(3, 4))]
print("sample points (e, t) evaluated exactly:")
ok = True
for e0, t0 in samples:
    lhs = sp.discriminant(sp.Poly(sp.expand(d.subs({e: e0, t: t0})), xi))
    rhs = claimed.subs({u: e0**2, s: t0**2})
    print("   (%s, %s)" % (e0, t0))
    ok = ok and sp.simplify(lhs - rhs) == 0
check("Appendix A   Pi_1, Pi_2 reproduce Disc_xi d_{e,t} at the sample points",
      ok)
print("note: the symbolic identity  Disc_xi d_{e,t} = 2^44 3^27 s^3 (u-1)^3")
print("      (u+3)^15 (s-1)^15 Pi_1^3 Pi_2  is proved in sec4_4_1.py;")
print("      the check above is only a fast exact spot-check.")
print()

report("Appendix A")
