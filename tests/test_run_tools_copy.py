"""The orchestrator copies a fixed list of tools/*.py into every run directory and points SINEDERELLA_TOOLS at that copy.
A module that one of them imports but that is not in the list makes the satellite stage die with ModuleNotFoundError inside
the run (found 2026-10-05: satellite_kindB_verify.py was imported since b9a6750 but not copied). This test reads the list from
the SINEderella script, copies exactly those files to a temporary directory and imports the stage from there; it also checks
that every sibling import of the copied tools is itself in the list."""
import ast
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.join(os.path.dirname(__file__), "..")
TOOLS = os.path.join(ROOT, "tools")


def copied_tools():
    """the file names of every `for _tool in ...; do` list of the SINEderella script (one per mode)"""
    txt = open(os.path.join(ROOT, "SINEderella"), encoding="utf-8").read()
    lists = []
    for m in re.finditer(r"for _tool in((?:[^;]|\\\n)*?); do", txt):
        body = re.sub(r"\\\n", " ", m.group(1))
        body = re.sub(r"#.*", "", body)
        lists.append(sorted(set(body.split())))
    return lists


def local_imports(path):
    """names imported at module level that are files in tools/"""
    tree = ast.parse(open(path, encoding="utf-8").read())
    out = set()
    for node in ast.walk(tree):
        names = []
        if isinstance(node, ast.Import):
            names = [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom) and node.module:
            names = [node.module]
        for n in names:
            if os.path.exists(os.path.join(TOOLS, n.split(".")[0] + ".py")):
                out.add(n.split(".")[0] + ".py")
    return out


class RunToolsCopy(unittest.TestCase):
    def test_lists_present_and_identical(self):
        lists = copied_tools()
        self.assertGreaterEqual(len(lists), 3, "full, add, exclude (and resume) modes each copy the tools")
        for l in lists[1:]:
            self.assertEqual(l, lists[0], "the copy lists of the modes differ")

    def test_every_sibling_import_is_copied(self):
        lst = set(copied_tools()[0])
        for f in sorted(lst):
            p = os.path.join(TOOLS, f)
            self.assertTrue(os.path.exists(p), "%s is in the copy list but not in tools/" % f)
            missing = local_imports(p) - lst
            self.assertFalse(missing, "%s imports %s, which the orchestrator does not copy" % (f, sorted(missing)))

    def test_stage_runs_from_the_copied_set(self):
        lst = copied_tools()[0]
        d = tempfile.mkdtemp()
        try:
            for f in lst:
                shutil.copy(os.path.join(TOOLS, f), d)
            r = subprocess.run([sys.executable, os.path.join(d, "satellite_stage.py"), "--help"], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr[-800:])
            self.assertIn("--exclude-b", r.stdout)
        finally:
            shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
