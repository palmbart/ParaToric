# python/paratoric/__init__.py
"""
Paratoric Python bindings.

Usage:
    from paratoric import get_sample, get_thermalization, get_hysteresis
"""

from __future__ import annotations

# Package version (when installed); fallback for dev checkouts.
try:
    from ._paratoric import __version__ 
except Exception:
    try:
        from importlib.metadata import version as _pkg_version
        __version__ = _pkg_version("paratoric")
    except Exception:
        __version__ = "0+local"

try:
    from . import _paratoric as _p
    extended_toric_code = _p.extended_toric_code
    __all__ = ["extended_toric_code", "__version__"]
except (ImportError, OSError) as _e:
    _CAUSE = _e
    # Helpful message if the extension isn't importable (not built or dlopen failure)
    _HINT = (
        "Could not import 'paratoric._paratoric'.\n"
        "Make sure the C++ extension is built and located at:\n"
        "  python/paratoric/_paratoric*.so (or .pyd on Windows)\n\n"
        "Quick fixes:\n"
        "  - Build with CMake so the .so is dropped into python/paratoric/\n"
        "  - Or `pip install -e python` inside your virtualenv\n"
        f"\nOriginal error:\n{_e}"
    )
    # Optional: expose a lazy attribute loader instead of failing immediately
    def __getattr__(name: str):
        if name in {"get_sample", "get_thermalization", "get_hysteresis"}:
            raise ImportError(_HINT) from _CAUSE
        raise AttributeError(name)

    __all__ = ["__version__"]  # still allow reading version

