# -*- coding: utf-8 -*-
"""
subsection 4.4.4: real roots and strict positivity of d_{e,t}.

Paper: subsection 4.4.4, equation (4.39).

The base-point polynomial

    d_0(xi) = (729/2) * d_{1/3,1/2}(xi)

uses the positive scale 729/2, preserving signs and real roots.

What is verified here
---------------------
* The coefficients of d_0 against the explicit polynomial of (4.39).
* Its constant term d_0(0) = 66199.

Appendix C proves that d_0 has no real root.  The general family d_{e,t}
additionally requires the nonzero leading coefficient, nonvanishing
discriminant and root continuation argument of Section 4.4.4; the base-point
Sturm certificate alone does not cover the full parameter domain.

Run:  python checks/sec4_4_4.py
"""
from _common import *

reset_marks()
banner('subsection 4.4.4', '(4.39)')

# Compare the derived base-point polynomial with (4.39).
show_coeffs("d_0(xi) = (729/2)*d_{1/3,1/2}(xi)", d0, xi)
check("(4.39)   d_0",
      sp.expand(d0 - (17914752*xi**6 - 5322240*xi**5 - 3190320*xi**4
                      + 3596208*xi**3 + 3949008*xi**2 + 94668*xi + 66199)) == 0)
check("(4.39)   d_0(0) = 66199 > 0", d0.subs(xi, 0) == 66199)
print()

report("subsection 4.4.4")
