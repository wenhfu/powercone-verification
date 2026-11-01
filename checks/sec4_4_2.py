# -*- coding: utf-8 -*-
"""
subsection 4.4.2: the sign of Pi_1.

Checks the algebraic steps in Lemma 4.1 and (4.10), (4.24)-(4.25),
(4.29)-(4.31). H is the trace-free part of T with respect to M.

Uses u = e^2, s = t^2, w = e*t with w^2 = u*s.
With D = det(M), the calculations use Adj = D*M^{-1},
S_i = D*tau_i, G_ijk = D*H_ijk, and N = D^5*||H||_M^2.

Run:  python checks/sec4_4_2.py
"""
from _common import *

reset_marks()
banner('subsection 4.4.2', '(4.10), Lemma 4.1, (4.24)-(4.25), (4.29)-(4.31)')

# Auxiliary variable w = e*t.
w = sp.Symbol("w", positive=True)


def to_usw(g):
    """Rewrite a polynomial using u = e^2, s = t^2, w = e*t.

    Reject monomials whose e and t exponents have different parity,
    and verify the result by reverse substitution.
    """
    original = sp.expand(g)
    out = 0
    for (i, j), co in sp.Poly(original, e, t).terms():
        if i % 2 != j % 2:
            raise ValueError("monomial e^%d t^%d has mismatched parity" % (i, j))
        if i % 2 == 0:
            out += co*u**(i//2)*s**(j//2)
        else:
            out += co*w*u**((i-1)//2)*s**((j-1)//2)
    converted = sp.expand(out)
    if sp.expand(converted.subs({u: e**2, s: t**2, w: e*t}) - original) != 0:
        raise ArithmeticError("reverse substitution failed in to_usw")
    return converted


def reduce_w(g):
    """Apply the relation w^2 = u*s to a polynomial in (u,s,w)."""
    g = sp.expand(g)
    if not g.has(w):
        return g
    P = sp.Poly(g, w)
    out = 0
    for k in range(P.degree(w) + 1):
        co = P.coeff_monomial(w**k)
        out += co*(u*s)**(k//2) if k % 2 == 0 else co*w*(u*s)**((k-1)//2)
    return sp.expand(out)


# ----------------------------------------------------------------------
# Symmetric matrix M and tensor T of Q and P in the direction coordinates (p,r).
# ----------------------------------------------------------------------
Q2 = sp.Poly(sp.expand(Q_der), p, r)
P2 = sp.Poly(sp.expand(P_der), p, r)

# Q = M_11 p^2 + 2 M_12 p r + M_22 r^2.
M = sp.Matrix([[to_usw(Q2.coeff_monomial(p**2)),
                to_usw(Q2.coeff_monomial(p*r)/2)],
               [to_usw(Q2.coeff_monomial(p*r)/2),
                to_usw(Q2.coeff_monomial(r**2))]])

# P = T_111 p^3 + 3 T_112 p^2 r + 3 T_122 p r^2 + T_222 r^3.
T = {(1, 1, 1): to_usw(P2.coeff_monomial(p**3)),
     (1, 1, 2): to_usw(P2.coeff_monomial(p**2*r)/3),
     (1, 2, 2): to_usw(P2.coeff_monomial(p*r**2)/3),
     (2, 2, 2): to_usw(P2.coeff_monomial(r**3))}


def Tget(i, j, k):
    """Symmetric access to the components of T."""
    return T[tuple(sorted((i, j, k)))]


# Adjugate of M, that is, det(M) * M^{-1}, together with det(M) itself.
Adj = sp.Matrix([[M[1, 1], -M[0, 1]], [-M[0, 1], M[0, 0]]])
D = reduce_w(M.det())
cc = 12*(1-s)*(u+3)                     # det(M) = cc * R(u,s)
R = (u+3)**2 - 4*u*(u+3)*s + 4*u*(1-u)*s**2

# (4.24)-(4.25): S_i = D*tau_i and G_ijk = D*H_ijk.
S = [reduce_w(sum(Adj[j-1, k-1]*Tget(i, j, k) for j in (1, 2) for k in (1, 2)))
     for i in (1, 2)]


def G(i, j, k):
    return reduce_w(D*Tget(i, j, k)
                    - sp.Rational(1, 4)*(M[i-1, j-1]*S[k-1]
                                         + M[i-1, k-1]*S[j-1]
                                         + M[j-1, k-1]*S[i-1]))


# (4.10): det M = 12(1-s)(u+3)R(u,s).
print("R(u,s) =", sp.factor(R))
print("det M  =", sp.factor(D))
check("(4.10)   det M = 12(1-s)(u+3)R(u,s)",
      sp.expand(D - cc*R) == 0)
print()

# (4.24)-(4.25): M^{jk} H_ijk = 0 is equivalent to Adj_jk G_ijk = 0.
mtrace = [reduce_w(sum(Adj[j-1, k-1]*G(i, j, k) for j in (1, 2) for k in (1, 2)))
          for i in (1, 2)]
check("(4.24)-(4.25)   H is M-traceless: M^{jk}H_ijk = 0",
      mtrace == [0, 0])
print()

# (4.30): H_111 = G_111 / D = A^3(1-s)^3.
H111 = sp.factor(sp.cancel(G(1, 1, 1)/D))
print("H_111 =", H111)
check("(4.30)   H_111 = A^3(1-s)^3",
      sp.expand(G(1, 1, 1) - D*(u+3)**3*(1-s)**3) == 0)
print()

# (4.28)-(4.29): N = D^5*||H||_M^2.
N = reduce_w(sum(Adj[i-1, a-1]*Adj[j-1, b-1]*Adj[k-1, c-1]
                 * G(i, j, k)*G(a, b, c)
                 for i in (1, 2) for j in (1, 2) for k in (1, 2)
                 for a in (1, 2) for b in (1, 2) for c in (1, 2)))
show_factors("Pi_1(u,s) = -N/(cc^3 det(M)^2)   [= ||H||_M^2 * (-R^3)]",
             sp.cancel(-N/(cc**3*D**2)))
check("(4.29)   ||H||_M^2 = -Pi_1(u,s)/R(u,s)^3",
      sp.expand(N + Pi1*cc**3*D**2) == 0)
print()

# ----------------------------------------------------------------------
# Lemma 4.1: R > 0 and M positive definite, needed for the sign of Pi_1.
# ----------------------------------------------------------------------
dR = sp.diff(R, s)
print("M_22    =", sp.factor(M[1, 1]))
print("dR/ds   =", sp.expand(dR))
print("R(u,1)  =", sp.factor(R.subs(s, 1)))
check("Lemma 4.1   M_22 = 4(4us+3u+9) > 0",
      sp.expand(M[1, 1] - 4*(4*u*s + 3*u + 9)) == 0)
check("Lemma 4.1   dR/ds = -4u(3u+1) - 8u(1-u)(1-s)",
      sp.expand(dR + 4*u*(3*u + 1) + 8*u*(1 - u)*(1 - s)) == 0)
check("Lemma 4.1   R(u,1) = (1-u)(7u+9)",
      sp.expand(R.subs(s, 1) - (1 - u)*(7*u + 9)) == 0)
# For 0<u,s<1, R(u,s) >= R(u,1) > 0, det M = cc*R > 0,
# and M_22 > 0 give positive definiteness.
qpd = (sp.expand(M[1, 1] - 4*(4*u*s + 3*u + 9)) == 0
       and sp.expand(dR + 4*u*(3*u + 1) + 8*u*(1 - u)*(1 - s)) == 0
       and sp.expand(R.subs(s, 1) - (1 - u)*(7*u + 9)) == 0
       and sp.expand(D - cc*R) == 0)
check("Lemma 4.1   R(u,s) > 0 and det M > 0  =>  M is positive definite", qpd)
print()

# (4.31): positive definiteness, H_111 != 0, and R > 0 imply Pi_1 < 0
# via (4.29) for 0<u,s<1.
check("(4.31)   Pi_1(u,s) < 0 for 0<u,s<1   (||H||_M^2 > 0, R(u,s) > 0)",
      qpd and sp.expand(G(1, 1, 1) - D*(u+3)**3*(1 - s)**3) == 0)
print()

report("subsection 4.4.2")
