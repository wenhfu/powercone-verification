# Symbolic verification of the upper-bound proof

This repository accompanies the manuscript:

> Rui-Jin Zhang, Wenhao Fu, and Yu-Hong Dai.  
> *The optimal barrier parameter of the three-dimensional power cone*.

The scripts check the algebraic identities and certificates used in the
upper-bound proof in Sections 4.3–4.4 and Appendices A–C.

## Installation and Usage

The code requires Python 3 (tested with Python **3.12.14**) and SymPy (tested with **1.14.0**).

From the repository root (the directory containing this README and `run_all.py`), install SymPy and run all ten verification modules:

```sh
python -m pip install sympy==1.14.0
python run_all.py
```

To run an individual module with detailed output, use:

```sh
python checks/sec4_4_4.py
```

## Output

All checks use exact symbolic arithmetic.

### Detailed output: a single script

Each script prints computed values and a `verified` or `FAILED` verdict for each check:

```
========================================================================
Checking subsection 4.4.4: (4.39)
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

See each script's documentation for the scope of its checks.

### Summary output: run_all.py

A completed run of `run_all.py` prints a summary:

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

## Scripts

The table lists the scripts and expected results; each script documents the checks it performs.

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
