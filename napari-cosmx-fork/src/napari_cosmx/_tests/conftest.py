import os
from pathlib import Path


def pytest_configure(config):
    """Configure Qt for headless Linux test runs before qapp starts."""
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    os.environ.setdefault("QT_API", "pyqt5")

    runtime_dir = Path("/tmp/pytest-qt-runtime")
    runtime_dir.mkdir(parents=True, exist_ok=True)
    runtime_dir.chmod(0o700)
    os.environ.setdefault("XDG_RUNTIME_DIR", str(runtime_dir))