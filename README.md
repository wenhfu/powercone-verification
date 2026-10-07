# Symbolic verification of the upper-bound argument

The Python code verifies the algebraic identities and symbolic computations appearing in *The optimal barrier parameter of the three-dimensional power cone*.

It covers subsections 4.3 and 4.4 and Appendices A-C of the upper-bound proof; the non-algebraic
steps and the lower-bound argument are proved by hand in the paper and are out of scope here.

## Requirements

* Python 3
* SymPy

## How to run

Run a single part of the paper:

```sh
python checks/sec4_3.py
```

Run everything:

```sh
python run_all.py
```

## Output

The two entry points check the same statements and differ only in how much they print.  Everything is
computed in exact arithmetic (rational numbers and polynomials only), and every verdict is an exact
decision: polynomial normalisation, Sturm root counting, Bernstein signs, or a remainder-checked
division.

### Detailed output: a single script

Running one script prints a banner naming the part of the paper and the equation numbers it checks,
then the value SymPy computes for each ingredient, and one verdict per statement:

```
========================================================================
Checking subsection 4.4.4: (4.39)
only SymPy-computed values are printed; compare them with the paper
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
========================================================================
```

Each printed value carries its equation number and a short handle, so that it can be read off against
the paper; a check begins with `verified`, or with `FAILED` when it does not hold.  A banner lists the
numbered statements that receive a verdict, not every equation used inside the script: when a verdict
combines several proof steps, all of its symbolic ingredients are checked here, while the elementary
step that joins them remains the paper's own.

### Summary output: run_all.py

`python run_all.py` runs the same scripts with their detailed output captured, and prints only the `OVERALL` line of each part, followed by the grand total:

```
########################################################################
# running all verification scripts
########################################################################
OVERALL [Appendix A]        PASS   (Pi_1, Pi_2 of (4.22), Appendix A)
OVERALL [Appendix B.1]      PASS   (A_6, A_7 of (4.33), Appendix B.1)
OVERALL [Appendix B.2]      PASS   (B.1)-(B.11)
OVERALL [Appendix C]        PASS   (Sturm chain of d_0 of (4.39), Appendix C)
OVERALL [subsection 4.3]    PASS   (4.7), (4.9)
OVERALL [subsection 4.4.0]  PASS   (4.7), (4.9), (4.13), (4.15), (4.17), (4.18), (4.20), (4.21)
OVERALL [subsection 4.4.1]  PASS   (4.22)
OVERALL [subsection 4.4.2]  PASS   (4.10), Lemma 4.1, (4.24), (4.29)-(4.31)
OVERALL [subsection 4.4.3]  PASS   (4.32)-(4.34), (4.37)
OVERALL [subsection 4.4.4]  PASS   (4.39)

########################################################################
# GRAND TOTAL: PASS
########################################################################
```

Only the detailed lines are discarded; the `OVERALL` lines, the `FAIL` lines and the grand total are
always printed, so no failure can be hidden.  Run a single script for the intermediate values of one
part.

## Scripts

| Script             | Part of the paper | Equations                                                    | Result |
| ------------------ | ----------------- | ------------------------------------------------------------ | ------ |
| `appendix_a.py`  | Appendix A        | $\Pi_1$, $\Pi_2$ of (4.22)                               | PASS   |
| `appendix_b1.py` | Appendix B.1      | $A_6$, $A_7$ of (4.33)                                   | PASS   |
| `appendix_b2.py` | Appendix B.2      | (B.1)-(B.11)                                                 | PASS   |
| `appendix_c.py`  | Appendix C        | Sturm chain of$d_0$ of (4.39)                              | PASS   |
| `sec4_3.py`      | subsection 4.3    | (4.7), (4.9)                                                 | PASS   |
| `sec4_4_0.py`    | subsection 4.4.0  | (4.7), (4.9), (4.13), (4.15), (4.17), (4.18), (4.20), (4.21) | PASS   |
| `sec4_4_1.py`    | subsection 4.4.1  | (4.22)                                                       | PASS   |
| `sec4_4_2.py`    | subsection 4.4.2  | (4.10), Lemma 4.1, (4.24), (4.29)-(4.31)                     | PASS   |
| `sec4_4_3.py`    | subsection 4.4.3  | (4.32)-(4.34), (4.37)                                        | PASS   |
| `sec4_4_4.py`    | subsection 4.4.4  | (4.39)                                                       | PASS   |

## License

Released under the MIT License.  See `LICENSE`.

Copyright (c) 2026 Rui-Jin Zhang, Wenhao Fu, and Yu-Hong Dai.
