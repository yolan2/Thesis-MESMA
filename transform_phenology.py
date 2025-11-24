#!/usr/bin/env python3
"""
Transform and prepare phenology CSV files for R MESMA script.

This script reads a combined CSV (`all_locations_timeseries.csv`) from a results directory,
computes spectral indices and performs some filtering and outlier removal, then writes
`hls_phenology_data.csv` which is the required input for `fit_veg_mixture_mesma.R`.

Usage: python transform_phenology.py --results-dir /path/to/results --out-dir /path/to/phenology_results

This implements the behavior you provided in the `process_phenology_data` method.
"""
from __future__ import annotations

import argparse
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import sys

import numpy as np
import pandas as pd


def get_optimal_workers(mode: str = "mixed") -> int:
    """Small helper to pick a reasonable number of worker threads.
    Keep it conservative so it runs on a laptop without saturating RAM.
    """
    try:
        import multiprocessing
        n = max(1, multiprocessing.cpu_count() - 1)
    except Exception:
        n = 2
    # keep a modest limit
    return min(n, 8)


class PhenologyTransformer:
    def __init__(self, results_dir: Path | None = None, pheno_out_dir: Path | None = None):
        self.results_dir = Path(results_dir) if results_dir is not None else None
        if pheno_out_dir is None:
            # default to the same path used in the R script
            self.pheno_out_dir = Path("C:/Users/yolan/OneDrive/Documenten/UGENT/Master/masterproef/phenology_results")
        else:
            self.pheno_out_dir = Path(pheno_out_dir)
        self.pheno_out_dir.mkdir(parents=True, exist_ok=True)

    def process_phenology_data(self, combined_csv: Path | None = None, require_latlon: bool = False) -> None:
        """Step: read combined CSV, compute indices, filter, and write hls_phenology_data.csv"""
        print("\n=== Step: Processing Phenology Data ===")

        # Allow explicit combined_csv override, otherwise default to workspace's expected file
        if combined_csv is None:
            # Default combined CSV: use the incremental results produced by MAP
            combined_csv = Path(r"C:\MAP\incremental_results.csv")

        # If a results_dir was provided and combined_csv points to a directory style, prefer results_dir / filename
        if self.results_dir is not None and (isinstance(combined_csv, Path) and combined_csv.name == "all_locations_timeseries.csv"):
            combined_csv = self.results_dir / "all_locations_timeseries.csv"

        if not Path(combined_csv).exists():
            raise FileNotFoundError(f"Combined results CSV not found: {combined_csv}")

        df = pd.read_csv(combined_csv)
        print(f"✓ Loaded {len(df)} observations from {combined_csv}")

        # Normalize column names
        df.columns = df.columns.str.lower()

        # Ensure date column
        if "prediction_date" not in df.columns and "date" in df.columns:
            df["prediction_date"] = df["date"]
        if "prediction_date" not in df.columns:
            raise RuntimeError("No prediction_date column found in input CSV")

        # Parse dates; require y-m-d format (fail loudly if cannot parse)
        df["date"] = pd.to_datetime(df["prediction_date"], format="%Y-%m-%d", errors="raise")
        df["year"] = df["date"].dt.year
        df["doy"] = df["date"].dt.dayofyear

        lon_candidates = [c for c in df.columns if c in ("lon", "longitude", "x", "imagery_lon", "target_lon")]
        lat_candidates = [c for c in df.columns if c in ("lat", "latitude", "y", "imagery_lat", "target_lat")]
        lon_col = lon_candidates[0] if len(lon_candidates) > 0 else None
        lat_col = lat_candidates[0] if len(lat_candidates) > 0 else None

        if lon_col is not None:
            df[lon_col] = pd.to_numeric(df[lon_col], errors="coerce")
        if lat_col is not None:
            df[lat_col] = pd.to_numeric(df[lat_col], errors="coerce")

        if "location_id" not in df.columns:
            raise RuntimeError("Input CSV must provide a 'location_id' column; deriving it from lon/lat has been removed.")

        df["location_id"] = df["location_id"].astype(str)
        print("✓ Using existing 'location_id' column from incremental_results.csv")

        # Rename band columns if needed
        band_mapping = {
            "band_blue": "blue",
            "band_green": "green",
            "band_red": "red",
            "band_nir": "nir",
            "band_swir1": "swir1",
            "band_swir2": "swir2",
        }
        df = df.rename(columns=band_mapping)

        # Ensure bands numeric
        bands = ["blue", "green", "red", "nir", "swir1", "swir2"]
        for band in bands:
            if band in df.columns:
                df[band] = pd.to_numeric(df[band], errors="coerce")

        # Verify that we have at least the necessary bands to compute indices
        if not ("nir" in df.columns and "red" in df.columns and "green" in df.columns and "blue" in df.columns and "swir1" in df.columns and "swir2" in df.columns):
            missing = [b for b in bands if b not in df.columns]
            raise RuntimeError(f"Missing required spectral band columns to compute indices: {missing}")

        # Compute spectral indices (numeric-safe using eps)
        eps = 1e-9
        df["DVI"] = df["nir"] - df["red"]
        df["OSAVI"] = (df["nir"] - df["red"]) / (df["nir"] + df["red"] + 0.16)
        # MCARI: ((R - G) - 0.2*(R - B)) * (R / (G + eps))
        df["MCARI"] = ((df["red"] - df["green"]) - 0.2 * (df["red"] - df["blue"])) * (df["red"] / (df["green"] + eps))
        df["CRI"] = (1.0 / (df["green"] + eps)) - (1.0 / (df["red"] + eps))
        df["PRI"] = (df["green"] - df["red"]) / (df["green"] + df["red"] + eps)
        df["NIRv"] = df["nir"] * ((df["nir"] - df["red"]) / (df["nir"] + df["red"] + eps))
        df["PSRI"] = (df["red"] - df["blue"]) / (df["nir"] + eps)
        df["NBR"] = (df["nir"] - df["swir2"]) / (df["nir"] + df["swir2"] + eps)
        df["TCW"] = (df["swir1"] - df["swir2"]) / (df["swir1"] + df["swir2"] + eps)
        df["TCG"] = (df["green"] - df["red"]) / (df["green"] + df["red"] + eps)
        df["MNDWI"] = (df["green"] - df["swir1"]) / (df["green"] + df["swir1"] + eps)
        df["DUSTI"] = (df["red"] - df["blue"]) / (df["red"] + df["blue"] + eps)

        # MSAVI: stable formulation
        # MSAVI = (2*N + 1 - sqrt((2N + 1)^2 - 8(N - R))) / 2
        df["MSAVI"] = (2 * df["nir"] + 1 - np.sqrt(np.maximum(0.0, (2 * df["nir"] + 1) ** 2 - 8 * (df["nir"] - df["red"])))) / 2

        # Apply small correction to NIRv (as in your snippet)
        df["NIRv"] = df["NIRv"] * 1.3

        # DUSTI filtering
        dusti_threshold = 0.5
        n_before = len(df)
        df = df[df["DUSTI"] <= dusti_threshold]
        n_after = len(df)
        if n_before > n_after:
            print(f"✓ Applied DUSTI filter: removed {n_before - n_after} observations with DUSTI > {dusti_threshold}")

        # Outlier detection per (location_id, year) using MAD-based robust filtering
        optimal_indices = [
            "DVI",
            "OSAVI",
            "MCARI",
            "CRI",
            "PRI",
            "NIRv",
            "PSRI",
            "NBR",
            "TCW",
            "TCG",
            "MNDWI",
        ]

        def outlier_detection(group: pd.DataFrame) -> pd.DataFrame:
            if len(group) < 10:
                return group
            mask = pd.Series(True, index=group.index)
            for idx in optimal_indices:
                if idx in group.columns:
                    values = group[idx].dropna()
                    if len(values) >= 5:
                        med = values.median()
                        mad_val = (values - med).abs().median() * 1.4826
                        threshold = 4.0 * mad_val
                        # if mad_val is zero (no dispersion), do not drop
                        if not np.isfinite(threshold) or threshold <= 0:
                            continue
                        mask &= (group[idx] - med).abs() <= threshold
            return group[mask]

        # Group and run outlier detection in parallel for speed
        groups = list(df.groupby(["location_id", "year"]))
        if len(groups) > 0:
            workers = get_optimal_workers("mixed")
            results = []
            with ThreadPoolExecutor(max_workers=workers) as executor:
                futures = {executor.submit(outlier_detection, grp.copy()): key for key, grp in groups}
                for f in as_completed(futures):
                    try:
                        res = f.result()
                        if res is not None and len(res) > 0:
                            results.append(res)
                    except Exception as exc:
                        key = futures[f]
                        print(f"⚠ Outlier filtering failed for {key}: {exc}")

            if len(results) > 0:
                df = pd.concat(results).reset_index(drop=True)

        # Add required metadata columns
        # Choose the best available lon/lat source columns for imagery/target values
        def pick_column(preferred_candidates):
            for c in preferred_candidates:
                if c in df.columns:
                    return c
            return None

        imagery_lon_src = pick_column(["imagery_lon", "lon", "longitude", "target_lon"])
        imagery_lat_src = pick_column(["imagery_lat", "lat", "latitude", "target_lat"])
        target_lon_src = pick_column(["target_lon", "lon", "longitude", "imagery_lon"])
        target_lat_src = pick_column(["target_lat", "lat", "latitude", "imagery_lat"])

        # Use selected sources or fill with pd.NA if none available
        df["imagery_lon"] = df[imagery_lon_src] if imagery_lon_src is not None else pd.NA
        df["imagery_lat"] = df[imagery_lat_src] if imagery_lat_src is not None else pd.NA
        df["target_lon"] = df[target_lon_src] if target_lon_src is not None else pd.NA
        df["target_lat"] = df[target_lat_src] if target_lat_src is not None else pd.NA

        # Reorder columns to match R script expectation
        meta_cols = [
            "location_id",
            "lat",
            "lon",
            "imagery_lat",
            "imagery_lon",
            "target_lat",
            "target_lon",
            "date",
            "year",
            "doy",
        ]

        # But keep actual column names (we used lat/lon candidates)
        # ensure 'lat'/'lon' exist as canonical names for output
        if "lat" not in df.columns:
            if lat_col is not None and lat_col in df.columns:
                df["lat"] = df[lat_col]
            else:
                df["lat"] = pd.NA
        if "lon" not in df.columns:
            if lon_col is not None and lon_col in df.columns:
                df["lon"] = df[lon_col]
            else:
                df["lon"] = pd.NA

        idx_cols = optimal_indices + ["DUSTI", "MSAVI"]
        other_cols = [col for col in df.columns if col not in set(meta_cols + idx_cols)]
        out_cols = meta_cols + idx_cols + other_cols
        out_cols = [c for c in out_cols if c in df.columns]
        df_out = df.loc[:, out_cols]

        # Save phenology data - this is the CSV the R script expects
        out_file = self.pheno_out_dir / "hls_phenology_data.csv"
        df_out.to_csv(out_file, index=False)
        print(f"✓ Saved phenology data: {out_file} ({len(df_out)} rows)")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Transform combined HLS timeseries to phenology CSV for MESMA R script")
    parser.add_argument("--results-dir", required=False, help="Directory containing all_locations_timeseries.csv")
    parser.add_argument("--combined-csv", required=False, help=r"Path to combined csv file (overrides --results-dir). Default: C:\Users\yolan\OneDrive\Documenten\UGENT\Master\masterproef\phenology_results\hls_phenology_data.csv")
    parser.add_argument("--out-dir", required=False, help="Directory to write hls_phenology_data.csv (default: R script path)")
    parser.add_argument("--require-latlon", required=False, action="store_true", help="If set, require lon/lat even if location_id exists in input CSV")

    args = parser.parse_args(argv)
    results_dir = Path(args.results_dir) if args.results_dir else None
    combined_csv = Path(args.combined_csv) if args.combined_csv else None
    out_dir = Path(args.out_dir) if args.out_dir else None

    try:
        transformer = PhenologyTransformer(results_dir=results_dir, pheno_out_dir=out_dir)
        transformer.process_phenology_data(combined_csv=combined_csv, require_latlon=args.require_latlon)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
