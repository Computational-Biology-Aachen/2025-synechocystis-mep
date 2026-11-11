"""Compile sympy expressions into model functions.

Until mxlpy 0.68 this module worked around two mxlpy limitations (ADR 0003): generated
functions calling `math.exp` couldn't be parsed back (MxlPy#202), and long sums were
lowered in quadratic time (MxlPy#203). Both are fixed in mxlpy 0.69, so `expr_fn` is a
thin wrapper around `mxlpy.meta.sympy_tools.compile_fn` that also returns the args
(and sets the function's module, see `expr_fn`).
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

from mxlpy.meta.sympy_tools import compile_fn

if TYPE_CHECKING:
    from collections.abc import Callable

    import sympy


def expr_fn(expr: sympy.Expr, origin: str) -> tuple[Callable[..., float], list[str]]:
    """Compile `expr` (free symbols = model names) into a function and its args."""
    args = sorted(str(s) for s in expr.free_symbols)
    fn = compile_fn(expr, args, name="_" + re.sub(r"\W", "_", origin))
    # compile_fn leaves __module__ = None, so cloudpickle searched all of sys.modules for
    # each function (~10 ms each, every time a parallel scan pickles the model). Named
    # after a real module that doesn't hold it, it is pickled by value right away.
    fn.__module__ = __name__
    return fn, args


def clip_zero(x: float) -> float:
    """max(x, 0), applied to a surrogate output (legacy clipping, ADR 0002)."""
    return max(x, 0.0)
