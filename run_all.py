# -*- coding: utf-8 -*-
"""
Run all scripts in ``checks``, printing their summaries and recorded failures.
Run a script directly for its detailed values and verdicts, e.g.

    python checks/sec4_4_4.py

Recorded check failures do not stop the remaining scripts.  The process exits
with status 0 if all checks pass and status 1 if any fails.  Each completed
script must report at least one check and leave no unreported checks.
Unexpected exceptions or exits, missing reports and invalid reports abort the
run with a nonzero status, writing the module name and captured output to
stderr.  An aborted run has no grand total.

Run:  python run_all.py
"""
import contextlib
import io
import os
import runpy
import sys

# Check modules, in execution order.
CHECKS = "checks"
FILES = ["appendix_a.py", "appendix_b1.py", "appendix_b2.py", "appendix_c.py",
         "sec4_3.py", "sec4_4_0.py", "sec4_4_1.py", "sec4_4_2.py",
         "sec4_4_3.py", "sec4_4_4.py"]

# Make _common importable when running the suite or an individual script.
here = os.path.dirname(os.path.abspath(__file__))
checks_dir = os.path.join(here, CHECKS)
for d in (here, checks_dir):
    if d not in sys.path:
        sys.path.insert(0, d)

print("#" * 72)
print("# running all verification scripts")
print("#" * 72)

failed = 0
_common = sys.modules.get("_common")
for name in FILES:
    if _common is not None:
        _common.last_report.clear()
    # Capture the detailed report and keep only the summary lines.
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            try:
                runpy.run_path(os.path.join(checks_dir, name), run_name="__main__")
            except SystemExit as exc:
                # report() exits with status 1 after recording a failed check.
                # Any other exit is unexpected, including a premature exit(0).
                common = sys.modules.get("_common")
                if (exc.code != 1 or common is None
                        or not common.last_report.get("bad")):
                    raise RuntimeError(
                        "%s exited unexpectedly with status %r" % (name, exc.code)
                    ) from exc
        _common = sys.modules.get("_common")
        if _common is None or not _common.last_report:
            raise RuntimeError("%s finished without calling report()" % name)
        if _common.last_report.get("n", 0) <= 0:
            raise RuntimeError("%s reported no checks" % name)
        if _common.ok_marks:
            raise RuntimeError("%s left checks out of its final report" % name)
    except BaseException:
        # The redirect context has restored stdout before diagnostics are sent.
        print("ERROR [%s]: run aborted; captured output follows." % name,
              file=sys.stderr, flush=True)
        print(buf.getvalue(), end="", file=sys.stderr, flush=True)
        raise
    failed += len(_common.last_report.get("bad", []))
    for line in buf.getvalue().splitlines():
        if line.startswith("OVERALL ") or line.startswith("   FAILED:"):
            print(line)

print()
print("#" * 72)
print("# GRAND TOTAL: %s" % ("PASS" if failed == 0 else "FAIL"))
print("#" * 72)
sys.exit(1 if failed else 0)
