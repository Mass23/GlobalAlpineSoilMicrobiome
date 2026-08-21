"""Export helpers – write covariate tables for downstream R workflows."""

from __future__ import annotations

from pathlib import Path
from typing import Union

import pandas as pd


def export_covariates(
    df: pd.DataFrame,
    path: Union[str, Path],
    fmt: str = "csv",
) -> Path:
    """Write a covariate DataFrame to disk.

    Parameters
    ----------
    df:
        DataFrame to export.
    path:
        Output file path.
    fmt:
        Output format: ``"csv"`` or ``"parquet"``.

    Returns
    -------
    Path
        Path to the written file.

    Raises
    ------
    ValueError
        If *fmt* is not ``"csv"`` or ``"parquet"``.
    """
    dest = Path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)

    if fmt == "csv":
        df.to_csv(dest, index=False)
    elif fmt == "parquet":
        df.to_parquet(dest, index=False)
    else:
        raise ValueError(f"fmt must be 'csv' or 'parquet', got {fmt!r}")

    return dest
