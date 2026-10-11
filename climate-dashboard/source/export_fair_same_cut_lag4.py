#!/usr/bin/env python3
"""Export matched current and four-year-lagged rows from each Oscar cutoff.

Run on Oscar, where the original cut1930 ... cut2025 trajectories are stored.
This reads existing FaIR outputs; it does not run FaIR or replace old exports.
Transfer the resulting single NPZ file for the local FaIR_comb_unB correction.
"""

import argparse
import json
from pathlib import Path

import numpy as np


def export_rows(input_dir, output_file, runs, start_year, end_year, ncols=841):
    years = np.arange(start_year, end_year + 1)
    if years.size == 0:
        raise ValueError("The end year must be at least the start year.")
    if output_file.exists():
        raise FileExistsError(f"Output already exists: {output_file}")

    payload = {
        "years": years,
        "current_boundary_years": years + 1,
        "lag4_boundary_years": years - 3,
        "runs": np.asarray(runs, dtype=int),
    }
    provenance = []
    for run in runs:
        current = np.full((len(years), ncols), np.nan)
        lag4 = np.full_like(current, np.nan)
        for i, year in enumerate(years):
            filename = f"all-headstails25r{run}cut{year}_temp_all.npy"
            source = input_dir / filename
            if not source.is_file():
                raise FileNotFoundError(
                    f"Missing cutoff trajectory: {source}. "
                    "Stopping rather than dropping a year or shifting labels."
                )
            arr = np.load(source, mmap_mode="r", allow_pickle=False)
            if arr.ndim != 2 or arr.shape[0] != 278 or not 0 < arr.shape[1] <= ncols:
                raise ValueError(f"Unexpected shape {arr.shape} in {source}; expected (278, <= {ncols}).")
            row_index = int(year) - 1750 + 1
            if not 4 <= row_index < arr.shape[0]:
                raise ValueError(f"Cutoff year {year} is outside the usable time bounds.")
            # Both rows come from THIS cutoff's posterior, with matching columns.
            current[i, :arr.shape[1]] = arr[row_index, :]
            lag4[i, :arr.shape[1]] = arr[row_index - 4, :]
            provenance.append({
                "file": filename,
                "current_row": row_index,
                "lag4_row": row_index - 4,
                "shape": list(arr.shape),
            })
        payload[f"r{run}_current_all"] = current
        payload[f"r{run}_lag4_all"] = lag4
        print(f"Prepared run {run}: {len(years)} matched cutoff/lag rows")

    payload["metadata_json"] = np.asarray(json.dumps({
        "input_dir": str(input_dir),
        "ensemble": "all-headstails25",
        "tag": "all",
        "timebounds_start": 1750,
        "lag_years": 4,
        "annual_label_convention": "following boundary: year - 1750 + 1",
        "sources": provenance,
    }))
    output_file.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_file, **payload)
    print(f"Saved {output_file}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path(
        "/users/jnickla1/data/jnickla1/fair-calibrate/output/all-headstails25"))
    parser.add_argument("--output", type=Path, default=Path("same_cut_all_lag4_headstails25.npz"))
    parser.add_argument("--first-run", type=int, default=0)
    parser.add_argument("--last-run", type=int, default=100)
    parser.add_argument("--start-year", type=int, default=1930)
    parser.add_argument("--end-year", type=int, default=2025)
    args = parser.parse_args()
    if args.first_run < 0 or args.last_run < args.first_run:
        parser.error("Use a non-negative, ascending run range.")
    if args.output.suffix != ".npz":
        parser.error("The output filename must end in .npz.")
    export_rows(args.input_dir, args.output,
                list(range(args.first_run, args.last_run + 1)),
                args.start_year, args.end_year)


if __name__ == "__main__":
    main()
