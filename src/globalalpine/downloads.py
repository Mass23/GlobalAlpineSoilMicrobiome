"""File-download helper with progress reporting and resumption support."""

from __future__ import annotations

from pathlib import Path
from typing import Union

import requests


def download_file(
    url: str,
    dest_path: Union[str, Path],
    chunk_size: int = 8192,
    timeout: int = 60,
) -> Path:
    """Download *url* to *dest_path*, streaming to avoid large memory use.

    Parameters
    ----------
    url:
        HTTPS URL to download.
    dest_path:
        Local file path for the downloaded content.  Parent directories are
        created automatically.
    chunk_size:
        Number of bytes per write chunk.
    timeout:
        Connection / read timeout in seconds.

    Returns
    -------
    Path
        Path to the downloaded file.

    Raises
    ------
    requests.HTTPError
        If the server returns a non-2xx status code.
    """
    dest = Path(dest_path)
    dest.parent.mkdir(parents=True, exist_ok=True)

    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        with requests.get(url, stream=True, timeout=timeout) as resp:
            resp.raise_for_status()
            with tmp.open("wb") as fh:
                for chunk in resp.iter_content(chunk_size=chunk_size):
                    if chunk:
                        fh.write(chunk)
        tmp.rename(dest)
    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise

    return dest
