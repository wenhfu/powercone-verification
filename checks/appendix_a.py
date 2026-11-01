# -*- coding: utf-8 -*-
"""
Appendix A: exact polynomial factors Pi_1 and Pi_2.

Displays the coefficients a_0,...,a_6 and b_0,...,b_9 and checks
the degrees, constant and leading coefficients, and values at s = 1.
The full discriminant factorisation (4.22) is checked in sec4_4_1.py.

Run:  python checks/appendix_a.py
"""
from _common import *

reset_marks()
banner('Appendix A', 'Pi_1, Pi_2 of (4.22), Appendix A')


def show_named_factors(name, prefix, expr, var):
    """Print the coefficients of ``expr`` in ``var``, one per line.

    The j-th coefficient is labelled ``prefix_j`` and printed in factored form,
    matching the naming used in Appendix A.
    """
    P = sp.Poly(sp.expand(expr), var)
    print("%s   (degree %d in %s, coefficients factored):"
          % (name, P.degree(var), var))
    for j in range(P.degree(var) + 1):
        print("   %s_%d(u) = %s" % (prefix, j, sp.factor(P.coeff_monomial(var**j))))


# Appendix A: coefficients a_0,...,a_6 of Pi_1.
show_named_factors("Pi_1(u,s) = sum_j a_j(u) s^j", "a", Pi1, s)
check("Appendix A   deg_s Pi_1 = 6  (j = 0..6)",
      sp.Poly(Pi1, s).degree(s) == 6)
check("Appendix A   a_0(u) = -4(u+3)^6",
      is_zero(Pi1.subs(s, 0) + 4*(u+3)**6))
check("Appendix A   a_6(u) = 16u^2(u-1)^2(u^2+2u+9)",
      is_zero(sp.Poly(Pi1, s).LC() - 16*u**2*(u-1)**2*(u**2+2*u+9)))
check("Appendix A   Pi_1(u,1) = u(u-1)^3(7u+9)^2",
      is_zero(Pi1.subs(s, 1) - u*(u-1)**3*(7*u+9)**2))
print()

# Appendix A: coefficients b_0,...,b_9 of Pi_2.
show_named_factors("Pi_2(u,s) = sum_j b_j(u) s^j", "b", Pi2, s)
check("Appendix A   deg_s Pi_2 = 9  (j = 0..9)",
      sp.Poly(Pi2, s).degree(s) == 9)
check("Appendix A   b_0(u) = 4(u+3)^9",
      is_zero(Pi2.subs(s, 0) - 4*(u+3)**9))
check("Appendix A   b_9(u) = 4(u-9)(u-1)^3(9u^3-19u^2+27u-81)",
      is_zero(sp.Poly(Pi2, s).LC()
              - 4*(u-9)*(u-1)**3*(9*u**3-19*u**2+27*u-81)))
check("Appendix A   Pi_2(u,1) = u(1-u)^6(7u+9)^2",
      is_zero(Pi2.subs(s, 1) - u*(1-u)**6*(7*u+9)**2))
print()

print("Full discriminant factorisation is verified in sec4_4_1.py.")
print()

report("Appendix A")
