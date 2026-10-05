# -*- coding: utf-8 -*-
"""
subsection 4.4.2 -- the sign of Pi_1.

Paper: subsection 4.4.2
(eq:detM, eq:tensor-trace, eq:Htensor, eq:tensor-M-norm, eq:Hnormidentity,
eq:H111).

Write the quadratic form of eq:Q as Q(h) = M_ij h_i h_j and the cubic form of
eq:P as P(h) = T_ijk h_i h_j h_k, where h = (p,r), M is symmetric and positive
definite, T is symmetric, M^{-1} = (M^{ij}) and repeated indices are summed.
The M-trace of T and the M-traceless part of T are

    tau_i = M^{jk} T_ijk,
    H_ijk = T_ijk - (1/4)( M_ij tau_k + M_ik tau_j + M_jk tau_i ).

The M-inner product of eq:tensor-M-norm,

    <S,T>_M = M^{ia}M^{jb}M^{kc} S_ijk T_abc,

is positive definite, because in an M-orthonormal basis it is the ordinary
sum of products of tensor components.  Hence ||H||_M^2 > 0 whenever H != 0,
and the sign of Pi_1 follows from two exact identities:

    H_111     = A^3 (1-s)^3,                      (eq:H111)
    ||H||_M^2 = -Pi_1(u,s) / R(u,s)^3.            (eq:Hnormidentity)

What is verified here
---------------------
* eq:detM         : det M = 12 (1-s)(u+3) R(u,s), with
                    R(u,s) = (u+3)^2 - 4u(u+3)s + 4u(1-u)s^2.  This is where
                    R comes from; by Lemma lem:Qpd it is positive on (0,1)^2.
* eq:tensor-trace : H is M-traceless, M^{jk} H_ijk = 0, as constructed.
* eq:H111         : H_111 = A^3 (1-s)^3.
* eq:Hnormidentity: ||H||_M^2 = -Pi_1(u,s)/R(u,s)^3.

How M^{-1} is handled
---------------------
Instead of inverting the matrix M symbolically, everything is multiplied
through by det(M).  With Adj = det(M) * M^{-1} the adjugate and

    S_i = Adj_jk T_ijk = det(M) tau_i,
    G_ijk = det(M) H_ijk,

the norm becomes

    ||H||_M^2 = N / det(M)^5,   N = Adj_ia Adj_jb Adj_kc G_ijk G_abc,

and, since det(M) = cc * R with cc = 12(1-s)(u+3),

    ||H||_M^2 = -Pi_1/R^3   <=>   N = -Pi_1 cc^3 det(M)^2.

Both sides of the last identity are polynomials, so the comparison is exact
and does not involve any rational-function simplification.

The alternative variables
-------------------------
Every quantity below is even in the barrier parameter e and in the slice
parameter t separately, so it is written in

    u = e^2,  s = t^2,  w = e*t,   subject to  w^2 = u*s.

Each of Q, P contributes only even-even or odd-odd monomials, which is why
this rewriting is faithful; working in (u,s,w) instead of (e,t) halves the
degrees and keeps the expansion small.  The relation w^2 = u*s is applied by
reduce_w().

Run:  python sec4_4_2.py
"""
from _common import *

reset_marks()
banner('subsection 4.4.2', 'eq:detM, eq:tensor-trace, eq:Htensor, eq:Hnormidentity, eq:H111')

# w stands for the product e*t; see the module docstring.
w = sp.Symbol("w", positive=True)


def to_usw(g):
    """Rewrite a polynomial of (e,t) through u = e^2, s = t^2, w = e*t.

    Every monomial of ``g`` must have exponents of the same parity in e and in
    t, so that it is either even-even (mapped to a polynomial in u,s) or
    odd-odd (mapped to w times a polynomial in u,s).
    """
    out = 0
    for (i, j), co in sp.Poly(sp.expand(g), e, t).terms():
        assert i % 2 == j % 2, "monomial e^%d t^%d is not parity-homogeneous" % (i, j)
        if i % 2 == 0:
            out += co*u**(i//2)*s**(j//2)
        else:
            out += co*w*u**((i-1)//2)*s**((j-1)//2)
    return sp.expand(out)


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
# Symmetric matrix M of Q and symmetric third-order tensor T of P, in the
# direction coordinates (h_1,h_2) = (p,r).
# ----------------------------------------------------------------------
Q2 = sp.Poly(sp.expand(Q_der), p, r)
P2 = sp.Poly(sp.expand(P_der), p, r)

# Q = M_11 p^2 + 2 M_12 p r + M_22 r^2.
M = sp.Matrix([[to_usw(Q2.coeff_monomial(p**2)),
                to_usw(Q2.coeff_monomial(p*r)/2)],
               [to_usw(Q2.coeff_monomial(p*r)/2),
                to_usw(Q2.coeff_monomial(r**2))]])

# P = T_111 p^3 + 3 T_112 p^2 r + 3 T_122 p r^2 + T_222 r^3; the factors 3
# count the orderings of the indices of a symmetric tensor.
T = {(1, 1, 1): to_usw(P2.coeff_monomial(p**3)),
     (1, 1, 2): to_usw(P2.coeff_monomial(p**2*r)/3),
     (1, 2, 2): to_usw(P2.coeff_monomial(p*r**2)/3),
     (2, 2, 2): to_usw(P2.coeff_monomial(r**3))}


def Tget(i, j, k):
    """Symmetric access to T."""
    return T[tuple(sorted((i, j, k)))]


# Adjugate of M, i.e. det(M) * M^{-1}, and det(M) itself.
Adj = sp.Matrix([[M[1, 1], -M[0, 1]], [-M[0, 1], M[0, 0]]])
D = reduce_w(M.det())
cc = 12*(1-s)*(u+3)                     # det(M) = cc * R(u,s)
R = (u+3)**2 - 4*u*(u+3)*s + 4*u*(1-u)*s**2

# M-trace of T (eq:tensor-trace) and G = det(M) * H (eq:Htensor), both
# multiplied through by det(M) to avoid inverting M.
S = [reduce_w(sum(Adj[j-1, k-1]*Tget(i, j, k) for j in (1, 2) for k in (1, 2)))
     for i in (1, 2)]


def G(i, j, k):
    return reduce_w(D*Tget(i, j, k)
                    - sp.Rational(1, 4)*(M[i-1, j-1]*S[k-1]
                                         + M[i-1, k-1]*S[j-1]
                                         + M[j-1, k-1]*S[i-1]))


# eq:detM: the determinant of M, which produces the factor R(u,s) of the paper.
print("R(u,s) =", sp.factor(R))
print("det M  =", sp.factor(D))
check("eq:detM   det M = 12(1-s)(u+3)R(u,s)",
      sp.expand(D - cc*R) == 0)
print()

# eq:tensor-trace: M^{jk} H_ijk = 0 is equivalent to Adj_jk G_ijk = 0.
mtrace = [reduce_w(sum(Adj[j-1, k-1]*G(i, j, k) for j in (1, 2) for k in (1, 2)))
          for i in (1, 2)]
check("eq:tensor-trace   H is M-traceless: M^{jk}H_ijk = 0",
      mtrace == [0, 0])
print()

# eq:H111: H_111 = G_111 / det(M).  Here tau_1 = 0 identically, i.e. the
# M-trace of T vanishes in the p-direction, so the correction term of
# eq:Htensor contributes nothing and H_111 = T_111 = A^3(1-s)^3.
H111 = sp.factor(sp.cancel(G(1, 1, 1)/D))
print("H_111 =", H111)
check("eq:H111   H_111 = A^3(1-s)^3",
      sp.expand(G(1, 1, 1) - D*(u+3)**3*(1-s)**3) == 0)
print()

# eq:tensor-M-norm and eq:Hnormidentity: the M-norm of H.
N = reduce_w(sum(Adj[i-1, a-1]*Adj[j-1, b-1]*Adj[k-1, c-1]
                 * G(i, j, k)*G(a, b, c)
                 for i in (1, 2) for j in (1, 2) for k in (1, 2)
                 for a in (1, 2) for b in (1, 2) for c in (1, 2)))
show_factors("Pi_1(u,s) = -N/(cc^3 det(M)^2)   [= ||H||_M^2 * (-R^3)]",
             sp.cancel(-N/(cc**3*D**2)))
check("eq:Hnormidentity   ||H||_M^2 = -Pi_1(u,s)/R(u,s)^3",
      sp.expand(N + Pi1*cc**3*D**2) == 0)
print()

report("subsection 4.4.2")
