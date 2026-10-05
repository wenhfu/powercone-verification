# -*- coding: utf-8 -*-
"""
Appendix B.2 -- exact base-point checks for Pi_2.

Paper: Appendix B.2
(eq:pi2quarter, eq:h1, eq:disch1, eq:pi2threequarter, eq:h2, eq:disch2,
eq:pi2half, eq:R7, eq:R7q, eq:h0, eq:disch0).

At the three base points u = 1/4, 1/2, 3/4 the polynomial 2^18 Pi_2(u, 1-q) is
split into a cubic part h_i(q) plus a tail whose coefficients are all
nonnegative for q >= 0.  Each cubic has negative discriminant, positive
leading coefficient and a positive value at q = 0, hence is positive for
q >= 0; this yields Pi_2(u,s) > 0 for 0 <= s <= 1.

What is verified here
---------------------
* The coefficient lists of 2^18 Pi_2(1/4,1-q) and of its cubic part h_1(q),
  and the same for u = 3/4 with h_2(q).
* The discriminant values Disc h_1 and Disc h_2 and the positivity argument
  (negative discriminant + positive leading coefficient + positive value at 0).
* For u = 1/2 : the splitting Pi_2(1/2,s) = -((2s-7)^2/512) R_7(s), the
  coefficients of R_7(s), the coefficients of -R_7(1-q) and of its cubic part
  h_0(q), and the corresponding discriminant argument.

Run:  python checks/appendix_b2.py
"""
from _common import *

reset_marks()
banner('Appendix B.2', 'Pi_2 base-point certificates (Appendix B.2)')

# ---- u = 1/4 ---------------------------------------------------------
# The scaled polynomial 2^18 Pi_2(1/4,1-q) and the decomposition of
# eq:pi2quarter: 27 h_1(q) plus a tail with nonnegative coefficients.
lhs1 = sp.expand(2**18*Pi2.subs({u: sp.Rational(1, 4), s: 1 - q}))
rhs1 = sp.expand(27*h1 + 6450327560*q**4 + 12966173284*q**5 + 12808009492*q**6
                 + 6652881036*q**7 + 1866026592*q**8 + 291453120*q**9)
show_coeffs("2^18*Pi_2(1/4,1-q)", lhs1, q)
show_coeffs("h_1(q)", h1, q)
check("eq:pi2quarter   decomposition for u = 1/4", is_zero(lhs1 - rhs1))
# eq:disch1 and the positivity of h_1 on q >= 0.
print("Disc h_1 =", sp.discriminant(sp.Poly(h1, q)))
check("eq:disch1   Disc h_1", sp.discriminant(sp.Poly(h1, q)) == -133708464211246463833953147)
check("Appendix B.2   h_1: Disc<0, LC>0, h_1(0)=49923>0  =>  h_1(q)>0 for q>=0",
      sp.discriminant(sp.Poly(h1, q)) < 0 and sp.Poly(h1, q).LC() > 0
      and h1.subs(q, 0) > 0)
print()

# ---- u = 3/4 ---------------------------------------------------------
# Same construction as for u = 1/4, with the cubic h_2 of eq:h2.
lhs2 = sp.expand(2**18*Pi2.subs({u: sp.Rational(3, 4), s: 1 - q}))
rhs2 = sp.expand(27*h2 + 54014298552*q**4 + 69669314028*q**5 + 25199957628*q**6
                 + 3569396868*q**7 + 234936288*q**8 + 9142848*q**9)
show_coeffs("2^18*Pi_2(3/4,1-q)", lhs2, q)
show_coeffs("h_2(q)", h2, q)
check("eq:pi2threequarter   decomposition for u = 3/4", is_zero(lhs2 - rhs2))
print("Disc h_2 =", sp.discriminant(sp.Poly(h2, q)))
check("eq:disch2   Disc h_2", sp.discriminant(sp.Poly(h2, q)) == -1772210425734147611691)
check("Appendix B.2   h_2: Disc<0, LC>0, h_2(0)=361>0  =>  h_2(q)>0 for q>=0",
      sp.discriminant(sp.Poly(h2, q)) < 0 and sp.Poly(h2, q).LC() > 0
      and h2.subs(q, 0) > 0)
print()

# ---- u = 1/2 ---------------------------------------------------------
# Here 2^18 Pi_2(1/2,.) has a double root at s = 7/2, which is removed by the
# factor (2s-7)^2:  Pi_2(1/2,s) = -((2s-7)^2/512) R_7(s).  Writing q = 1 - s,
# -R_7(1-q) splits into the cubic h_0(q) plus a tail with positive
# coefficients.
pi2_half = sp.expand(Pi2.subs(u, sp.Rational(1, 2)))
R7_rec = sp.expand(sp.cancel(-512*pi2_half/(2*s - 7)**2))
show_coeffs("R_7(s)", R7, s)
show_coeffs("-R_7(1-q) recomputed from Pi_2(1/2,s)",
            sp.expand(-R7_rec.subs(s, 1 - q)), q)
check("eq:pi2half   Pi_2(1/2,s) = -((2s-7)^2/512) R_7(s)",
      is_zero(pi2_half + sp.Rational(1, 512)*(2*s - 7)**2*R7))
check("eq:R7q   -R_7(1-q) = h_0(q) + positive higher-order terms",
      is_zero(sp.expand(-R7_rec.subs(s, 1 - q))
              - sp.expand(h0 + 1359652*q**4 + 1485064*q**5
                          + 268996*q**6 + 38692*q**7)))
show_coeffs("h_0(q)", h0, q)
print("Disc h_0 =", sp.discriminant(sp.Poly(h0, q)))
check("eq:disch0   Disc h_0", sp.discriminant(sp.Poly(h0, q)) == -188304545926875)
check("Appendix B.2   h_0: Disc<0, LC>0, h_0(0)=25>0  =>  h_0(q)>0 for q>=0",
      sp.discriminant(sp.Poly(h0, q)) < 0 and sp.Poly(h0, q).LC() > 0
      and h0.subs(q, 0) > 0)
print()

report("Appendix B.2")
