# -*- coding: utf-8 -*-
"""
Run the whole verification suite, in the order in which the material appears
in the paper.

A single script prints every value it computes together with the verdict for
each labelled statement of the paper.  This driver is quiet: it runs each
script with its output captured and prints only the per-part summary, so the
whole suite can be checked at a glance.  Nothing is hidden by this, because
only the printed details are discarded, never the verdict.

For the detailed report of one part, run that script directly, e.g.

    python sec4_4_4.py

Run:  python run_all.py
"""
import contextlib
import io
import os
import runpy
import sys

# Paper order.
FILES = ["appendix_a.py", "appendix_b1.py", "appendix_b2.py", "appendix_c.py",
         "sec4_3.py", "sec4_4_0.py", "sec4_4_1.py", "sec4_4_2.py",
         "sec4_4_3.py", "sec4_4_4.py"]

# Run each script as if it had been started from this directory, so that the
# shared module _common is importable.
here = os.path.dirname(os.path.abspath(__file__))
if here not in sys.path:
    sys.path.insert(0, here)

print("#" * 72)
print("# running all verification scripts")
print("#" * 72)

total = failed = 0
_common = None
for name in FILES:
    # Capture the detailed report and keep only the summary lines.
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        runpy.run_path(os.path.join(here, name), run_name="__main__")
    if _common is None:
        _common = sys.modules["_common"]
    # Accumulate the per-section tally recorded by _common.report().
    total += _common.last_report.get("n", 0)
    failed += len(_common.last_report.get("bad", []))
    for line in buf.getvalue().splitlines():
        if line.startswith("OVERALL ") or line.startswith("   FAILED:"):
            print(line)

print()
print("#" * 72)
print("# GRAND TOTAL: %s   (%d items checked, %d failed)"
      % ("ALL OK" if failed == 0 else "FAILED", total, failed))
print("#" * 72)
