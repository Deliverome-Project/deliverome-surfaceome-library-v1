"""Environment loading for deliverome-analysis.

Reads ``.env`` at the repo root, with **shell-environment precedence** — a
variable already set in the shell (or by CI) wins over the same key in
``.env``. Matches the accessible-surfaceome convention so muscle memory
transfers.

Usage
-----
    from deliverome_analysis.env import load_env
    load_env()                  # call once at script / notebook start

Or pass a specific ``.env`` file path:

    load_env("/Users/me/.deliverome/.env")
"""

from __future__ import annotations

import os
from pathlib import Path

from dotenv import dotenv_values

_REPO_ROOT = Path(__file__).resolve().parents[2]
_DEFAULT_ENV_PATH = _REPO_ROOT / ".env"
_loaded = False


def load_env(path: str | Path | None = None, *, override: bool = False) -> dict[str, str]:
    """Load ``.env`` into ``os.environ``.

    Idempotent — calling it twice is a no-op unless ``override`` is set.

    Parameters
    ----------
    path
        Path to a ``.env`` file. Defaults to ``<repo_root>/.env``.
    override
        If True, re-load and overwrite ``os.environ`` values that already
        exist. Default False — shell env always wins.

    Returns
    -------
    dict[str, str]
        The variables that were actually applied to ``os.environ`` on this
        call (useful for assertion-style checks in notebooks).
    """
    global _loaded
    if _loaded and not override:
        return {}

    env_path = Path(path) if path else _DEFAULT_ENV_PATH
    if not env_path.exists():
        # Silent — many scripts (CI, container) supply env via shell instead.
        _loaded = True
        return {}

    applied: dict[str, str] = {}
    for key, value in dotenv_values(env_path).items():
        if value is None:
            continue
        if not override and key in os.environ:
            continue
        os.environ[key] = value
        applied[key] = value

    _loaded = True
    return applied


def require(*keys: str) -> None:
    """Raise RuntimeError if any of ``keys`` is missing from the environment.

    Use at the top of scripts that need a particular set of credentials:

        require("S3_ACCESS_KEY_ID", "S3_SECRET_ACCESS_KEY")
    """
    load_env()
    missing = [k for k in keys if not os.environ.get(k)]
    if missing:
        joined = ", ".join(missing)
        raise RuntimeError(
            f"Missing required environment variable(s): {joined}. "
            f"Set them in .env (see .env.example) or export them in the shell."
        )
