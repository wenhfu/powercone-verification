# -*- coding: utf-8 -*-
"""
subsection 4.4.1: factorisation of the discriminant of d_{e,t}.

Paper: subsection 4.4.1, equation (4.22), together with the explicit factors
Pi_1, Pi_2 of Appendix A.

The proof of the upper bound rests on the factorisation

    Disc_xi d_{e,t}
        = 2^44 3^27 s^3 (u-1)^3 (u+3)^15 (s-1)^15 Pi_1(u,s)^3 Pi_2(u,s),

where u = e^2, s = t^2 and Pi_1, Pi_2 are the polynomials of Appendix A.

What is verified here
---------------------
* The discriminant identity  Disc_xi d_{e,t} = -Res_xi(d, d') / [xi^6] d_{e,t}
  for a sextic.
* Evenness in e and t before replacing their squared powers by u and s,
  together with an exact check by reverse substitution.
* The exact factorisation of Disc_xi d_{e,t}, symbolically in (u,s), against
  the right-hand side above.

This is the most expensive computation of the suite, and the only place where
the full discriminant factorization is checked.  The Appendix A script prints
the encoded coefficient lists and checks degrees, extremal coefficients and
boundary values; it does not independently compare every coefficient.

Run:  python checks/sec4_4_1.py
"""
from _common import *

reset_marks()
banner('subsection 4.4.1', '(4.22)')

# For a polynomial of degree n with leading coefficient a_n one has
# Disc = (-1)^{n(n-1)/2} Res(p, p') / a_n; here n = 6, so the sign is -1 and
# the identity below is checked in the form Res + Disc * leading = 0.
res = sp.resultant(sp.Poly(d, xi), sp.Poly(sp.diff(d, xi), xi))
disc_et = sp.expand(sp.discriminant(sp.Poly(d, xi)))
check("(4.22)   Disc_xi d_{e,t} = -Res_xi(d, d')/[xi^6]d_{e,t}",
      is_zero(res + disc_et*leading))
print()

# Rewrite the discriminant in the variables u = e^2, s = t^2 (every monomial is
# an even power of e and of t), factor it, and compare with the paper's form.
Pdisc = sp.Poly(disc_et, e, t)
terms = Pdisc.terms()
if any(i % 2 or j % 2 for (i, j), co in terms):
    raise ValueError("The discriminant must be even in both e and t")
disc_us = sp.expand(sum(co*u**(i//2)*s**(j//2) for (i, j), co in terms))
check("(4.22)   substitution u=e^2, s=t^2",
      sp.expand(disc_us.subs({u: e**2, s: t**2}) - disc_et) == 0)
show_factors("factor(Disc_xi d_{e,t})", disc_us)
claimed = 2**44*3**27*s**3*(u-1)**3*(u+3)**15*(s-1)**15*Pi1**3*Pi2
check("(4.22)   Disc_xi d_{e,t} factorisation, with Pi_1, Pi_2 of Appendix A",
      sp.expand(disc_us - claimed) == 0)
print()

report("subsection 4.4.1")
