# -*- coding: utf-8 -*-
"""
subsection 4.4.4 -- real roots and strict positivity of d_{e,t}.

Paper: subsection 4.4.4 (eq:basepoly).

The base-point polynomial

    d_0(xi) = (729/2) * d_{1/3,1/2}(xi)

is the specialisation of d_{e,t} at (e,t) = (1/3,1/2) scaled by the positive
rational 729/2.  It is displayed explicitly as eq:basepoly, and its Sturm
chain is analysed in Appendix C, which shows that it has no real root.

What is verified here
---------------------
* The coefficients of d_0 against the explicit polynomial of eq:basepoly.
* Its constant term d_0(0) = 66199.

The absence of real roots (and hence the strict positivity of d_{e,t} on the
relevant domain) is then established in appendix_c.py from the Sturm chain.

Run:  python checks/sec4_4_4.py
"""
from _common import *

reset_marks()
banner('subsection 4.4.4', 'eq:basepoly')

# d_0 is defined as (729/2) * d_{1/3,1/2}; the multiplication by the positive
# constant does not change the sign or the number of real roots.
show_coeffs("d_0(xi) = (729/2)*d_{1/3,1/2}(xi)", d0, xi)
check("eq:basepoly   d_0",
      sp.expand(d0 - (17914752*xi**6 - 5322240*xi**5 - 3190320*xi**4
                      + 3596208*xi**3 + 3949008*xi**2 + 94668*xi + 66199)) == 0)
check("eq:basepoly   d_0(0) = 66199 > 0", d0.subs(xi, 0) == 66199)
print()

report("subsection 4.4.4")
