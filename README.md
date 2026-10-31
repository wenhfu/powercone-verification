# Symbolic verification of the upper-bound argument

The Python code verifies the algebraic identities and symbolic computations appearing in ***The optimal barrier parameter of the three-dimensional power cone***.

It covers subsections 4.3 and 4.4 and Appendices A-C of the upper-bound proof; the non-algebraic
steps and the lower-bound argument are proved by hand in the paper and are out of scope here.

## Requirements

* Python 3; tested with Python **3.12.14**.
* SymPy; tested with SymPy **1.14.0**.

These are the tested versions, not a claim that every Python 3 or SymPy version has been tested.

## How to run

From the repository root (the directory containing this README and `run_all.py`), install the
tested dependency and run all ten modules:

```sh
python -m pip install sympy==1.14.0
python run_all.py
```

From the same directory, run a single part of the paper with detailed output:

```sh
python checks/sec4_3.py
```

## Output

For each module, running it directly and through `run_all.py` executes the same checks.  The full
runner executes all ten modules and normally prints only their summaries and recorded failure
summaries; detailed computed values are suppressed in this mode, including those associated with
recorded failures.

Symbolic expressions, including the barrier's logarithms and symbolic powers, are manipulated and
differentiated exactly.  The final algebraic verdicts use exact rational and polynomial arithmetic:
polynomial normalisation, exact root counting, Bernstein coefficient signs, or remainder-checked
division, without floating-point sampling.

The examples and PASS table below show **expected results for the accompanying manuscript snapshot
dated 2026-10-31**.  Equation numbers refer to that snapshot.  For a published release, identify the
manuscript version together with the code release/tag or commit containing the tested files.

### Detailed output: a single script

Running one script prints a banner naming the part of the paper and the equation numbers it checks,
then the value SymPy computes for each ingredient, and one verdict per statement:

```
========================================================================
Checking subsection 4.4.4: (4.39)
symbolic values and encoded coefficients; compare them with the paper
========================================================================

d_0(xi) = (729/2)*d_{1/3,1/2}(xi)   (degree 6 in xi):
   [xi^6] = 17914752
   [xi^5] = -5322240
   [xi^4] = -3190320
   [xi^3] = 3596208
   [xi^2] = 3949008
   [xi^1] = 94668
   [xi^0] = 66199
   verified   (4.39)   d_0
   verified   (4.39)   d_0(0) = 66199 > 0

========================================================================
OVERALL [subsection 4.4.4]  PASS   (4.39)
```

Each group of computed values or encoded coefficient lists is accompanied by a verdict identifying
the relevant equation or appendix.  A verdict begins with `verified`, or with `FAILED` when its
condition does not hold.  Printing an encoded list is not an independent comparison of every entry
with the manuscript.  Each module's documentation distinguishes displayed inputs from the identities
and certificates it checks.  The non-algebraic steps combining those ingredients remain the paper's
own arguments.

### Summary output: run_all.py

For completed scripts, `python run_all.py` prints each `OVERALL` line and any recorded failure
summaries.  If every script completes its report, the runner then prints the grand total:

```
########################################################################
# running all verification scripts
########################################################################
OVERALL [Appendix A]        PASS   (Pi_1, Pi_2 of (4.22), Appendix A)
OVERALL [Appendix B.1]      PASS   (A_6, A_7 of (4.33), Appendix B.1)
OVERALL [Appendix B.2]      PASS   (B.1)-(B.11)
OVERALL [Appendix C]        PASS   (Sturm chain of d_0 of (4.39), Appendix C)
OVERALL [subsection 4.3]    PASS   (4.7), (4.9)
OVERALL [Section 4.4]       PASS   (4.7), (4.9), (4.13), (4.15), (4.17), (4.18), (4.20), (4.21)
OVERALL [subsection 4.4.1]  PASS   (4.22)
OVERALL [subsection 4.4.2]  PASS   (4.10), Lemma 4.1, (4.24)-(4.25), (4.29)-(4.31)
OVERALL [subsection 4.4.3]  PASS   (4.32)-(4.34), (4.37)
OVERALL [subsection 4.4.4]  PASS   (4.39)

########################################################################
# GRAND TOTAL: PASS
########################################################################
```

Both entry points exit with status `0` when all their checks pass and status `1` when a recorded
check fails.  The full runner continues after recorded check failures so that its final summary
includes every completed part.

Unexpected exceptions or premature script exits abort the run with a nonzero status.  Before
propagating the error, the runner writes the module name and all of that module's captured output to
standard error, followed by the exception traceback.  An aborted run has no grand total.  Run a
single script to inspect the detailed intermediate values of one part.

Before replacing powers of `e` and `t` by powers of `u = e^2` and `s = t^2`, the leading-coefficient
and discriminant checks explicitly verify that all exponents are even.  Both substitutions are also
checked by exact reverse substitution.

The tensor calculation separately rewrites coefficients whose exponents have the same parity using
`u = e^2`, `s = t^2`, and `w = et`.  Each conversion explicitly checks this parity and verifies its
reverse substitution, including under Python `-O`.  The Bernstein helpers reject requested degrees
that are negative, nonintegral, or too small for the input polynomial.

## Scripts

The table records expected results for the manuscript snapshot identified above.  The following
scope distinctions are relevant when running modules separately:

* `appendix_a.py` prints the 17 encoded coefficient polynomials of `Pi_1` and `Pi_2` and performs
  eight checks of degrees, extremal coefficients, and boundary values.  It does not independently
  compare all 17 entries with a second list.  `sec4_4_1.py` computes the sextic discriminant and
  verifies its complete factorization involving those polynomials.
* `appendix_b1.py` displays the encoded power-basis coefficients of `A_6` and `A_7`, recomputes
  their Bernstein coefficients on two subintervals, compares the four reference lists, and checks
  their positivity.  The subdivision obtains positive Bernstein certificates; the factor `2u-1`
  in the discriminant of `Pi_2` separately determines the continuation intervals.
* `sec4_4_4.py` verifies the base-point polynomial `d_0`; Appendix C establishes that it has no
  real root.  The general conclusion for `d_{e,t}` also uses the nonvanishing leading coefficient
  and discriminant and the paper's real-root continuation argument.

| Script             | Part of the paper | Equations                                                    | Expected result |
| ------------------ | ----------------- | ------------------------------------------------------------ | ------ |
| `appendix_a.py`  | Appendix A        | $\Pi_1$, $\Pi_2$ of (4.22)                               | PASS   |
| `appendix_b1.py` | Appendix B.1      | $A_6$, $A_7$ of (4.33)                                   | PASS   |
| `appendix_b2.py` | Appendix B.2      | (B.1)-(B.11)                                                 | PASS   |
| `appendix_c.py`  | Appendix C        | Sturm chain of $d_0$ of (4.39)                              | PASS   |
| `sec4_3.py`      | subsection 4.3    | (4.7), (4.9)                                                 | PASS   |
| `sec4_4_0.py`    | Preliminary reductions in Section 4.4 | (4.7), (4.9), (4.13), (4.15), (4.17), (4.18), (4.20), (4.21) | PASS   |
| `sec4_4_1.py`    | subsection 4.4.1  | (4.22)                                                       | PASS   |
| `sec4_4_2.py`    | subsection 4.4.2  | (4.10), Lemma 4.1, (4.24)-(4.25), (4.29)-(4.31)                     | PASS   |
| `sec4_4_3.py`    | subsection 4.4.3  | (4.32)-(4.34), (4.37)                                        | PASS   |
| `sec4_4_4.py`    | subsection 4.4.4  | (4.39)                                                       | PASS   |

## License

Released under the MIT License.  See `LICENSE`.

Copyright (c) 2026 Rui-Jin Zhang, Wenhao Fu, and Yu-Hong Dai.
