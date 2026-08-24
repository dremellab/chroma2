#!/usr/bin/env python3
"""
Helpers shared by the MultiQC custom-content generator scripts
(make_*_mqc.py, aggregate_fragment_sizes_mqc.py).
"""

from __future__ import annotations

import os
from typing import IO, Mapping


def write_mqc_header(
    f: IO[str],
    id: str,
    section_name: str,
    plot_type: str,
    pconfig: Mapping[str, str],
) -> None:
    """Write the '# id:'/'# section_name:'/'# plot_type:'/'# pconfig:' comment
    header MultiQC custom-content expects, followed by one '#   key: value'
    line per pconfig entry (in the given order). Values are single-quoted to
    match the existing convention -- pass plain strings, not pre-quoted ones.
    """
    f.write(f"# id: '{id}'\n")
    f.write(f"# section_name: '{section_name}'\n")
    f.write(f"# plot_type: '{plot_type}'\n")
    f.write("# pconfig:\n")
    for key, value in pconfig.items():
        f.write(f"#   {key}: '{value}'\n")


def sample_name_from_path(path: str) -> str:
    """Sample name convention used across the *_mqc scripts: the first
    dot-separated token of the input file's basename."""
    return os.path.basename(path).split(".")[0]
