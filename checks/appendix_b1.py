# -*- coding: utf-8 -*-
"""
Appendix B.1 -- positivity of the auxiliary factors A_6 and A_7.

Paper: Appendix B.1 (eq:A6, eq:A7 and the Bernstein coefficient lists below
them).

A_6 and A_7 are the two nontrivial factors of Disc_s Pi_2 (see sec4_4_3.py).
Both are shown to be positive on [0,1] by exhibiting their Bernstein
coefficients, which are all positive.

What is verified here
---------------------
* The coefficient lists of A_6(u) and A_7(u) against eq:A6 / eq:A7.
* The Bernstein coefficients of A_6(v/2), A_6((1+v)/2), A_7(v/2) and
  A_7((1+v)/2) against the four lists printed in Appendix B.1.
* Positivity of A_6 and A_7 on [0,1] from the positivity of all Bernstein
  coefficients on each half-interval.

The half-interval split at 1/2 is used because 2u-1 enters the discriminant of
Pi_2 as a factor; the two halves are handled separately.

Run:  python checks/appendix_b1.py
"""
from _common import *

reset_marks()
banner('Appendix B.1', 'A_6, A_7 positivity (Appendix B.1)')

# The two auxiliary polynomials themselves, coefficient by coefficient.
show_coeffs("A_6(u)", A6, u)
show_coeffs("A_7(u)", A7, u)
print()

# Bernstein coefficients of A_6 and A_7 of the proper degree on each half of
# [0,1]; the change of variable u = v/2 (resp. u = (1+v)/2) maps [0,1] onto
# the half-interval.
cA6a = bernstein_1d(A6, u, 0, sp.Rational(1, 2), 6)
cA6b = bernstein_1d(A6, u, sp.Rational(1, 2), 1, 6)
cA7a = bernstein_1d(A7, u, 0, sp.Rational(1, 2), 7)
cA7b = bernstein_1d(A7, u, sp.Rational(1, 2), 1, 7)

# The four lists as printed in Appendix B.1.
A6a_pap = [2187, sp.Rational(1215, 2), sp.Rational(729, 10),
           sp.Rational(1899, 10), sp.Rational(171403, 240),
           sp.Rational(48309, 32), sp.Rational(40527, 16)]
A6b_pap = [sp.Rational(40527, 16), sp.Rational(113799, 32),
           sp.Rational(1153753, 240), sp.Rational(31214, 5),
           7912, sp.Rational(29824, 3), 12544]
A7a_pap = [1679616, sp.Rational(5371272, 7), sp.Rational(897831, 4),
           sp.Rational(1408707, 56), sp.Rational(1177189, 8),
           sp.Rational(64156205, 112), sp.Rational(577617375, 448),
           sp.Rational(292871875, 128)]
A7b_pap = [sp.Rational(292871875, 128), sp.Rational(736242875, 224),
           sp.Rational(1023180785, 224), sp.Rational(42856724, 7),
           sp.Rational(55624276, 7), sp.Rational(70273152, 7),
           sp.Rational(86818176, 7), 15040512]

print("Bernstein coefficients of A_6(v/2)     :", [str(x) for x in cA6a])
print("Bernstein coefficients of A_6((1+v)/2) :", [str(x) for x in cA6b])
print("Bernstein coefficients of A_7(v/2)     :", [str(x) for x in cA7a])
print("Bernstein coefficients of A_7((1+v)/2) :", [str(x) for x in cA7b])
print()
check("Appendix B.1   A_6(v/2) Bernstein coefficients", list(cA6a) == A6a_pap)
check("Appendix B.1   A_6((1+v)/2) Bernstein coefficients", list(cA6b) == A6b_pap)
check("Appendix B.1   A_7(v/2) Bernstein coefficients", list(cA7a) == A7a_pap)
check("Appendix B.1   A_7((1+v)/2) Bernstein coefficients", list(cA7b) == A7b_pap)
# A polynomial all of whose Bernstein coefficients are > 0 is > 0 on the
# interval; hence A_6 > 0 and A_7 > 0 on [0,1].
check("Appendix B.1   A_6(u) > 0 on [0,1]  (all Bernstein coefficients > 0)",
      all(x > 0 for x in cA6a + cA6b))
check("Appendix B.1   A_7(u) > 0 on [0,1]  (all Bernstein coefficients > 0)",
      all(x > 0 for x in cA7a + cA7b))
print()

report("Appendix B.1")
