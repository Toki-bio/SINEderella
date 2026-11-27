"""report_profile.profile on a plate where some windows are empty in every kept row (Atelerix e1-4 top100, run
ten11_nodash 2026-11-27): those windows have no median, and before the fix their NaN reached np.linalg.svd, which raised
"SVD did not converge" and stopped step 6. Now they are left out of the SVD and keep NaN in the mosaic track."""
import sys
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import report_profile as RP  # noqa: E402

PLATE = ROOT / "tests" / "data" / "aal_e1-4_top100_svd_nan.aln.fa"


def test_profile_survives_all_nan_windows():
    assert RP.profile(str(PLATE)) is not None


def test_svd_input_is_finite(monkeypatch):
    seen = []
    orig = np.linalg.svd

    def spy(a, *k, **kw):
        a = np.asarray(a)
        seen.append(bool(np.isfinite(a).all()))
        return orig(a, *k, **kw)

    monkeypatch.setattr(np.linalg, "svd", spy)
    RP.profile(str(PLATE))
    assert seen and all(seen)
