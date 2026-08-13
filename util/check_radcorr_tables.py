#!/usr/bin/env python3
"""Validate the complete 3He radiative-weight archive or extracted tables."""

import argparse
import csv
import io
import math
import re
import sys
import tarfile
from collections import defaultdict
from pathlib import Path

ARCHIVE_NAME = "Newfit_20260710_fullxquad_15angles.tar.gz"
ANGLES = tuple(26.0 + 0.5 * index for index in range(15))
MODELS = tuple(f"SF{index}" for index in range(1, 6))
REQUIRED = ("Ebeam", "nu", "XSborn_unp", "XSrad_unp")
RAD_DNU_MEV = 1.0
RAD_DNU_TOL = 1.0e-6


def csv_sources(location: Path):
    if location.is_file():
        with tarfile.open(location, "r:gz") as archive:
            for member in archive.getmembers():
                if member.isfile() and member.name.endswith("_short.csv"):
                    stream = archive.extractfile(member)
                    yield member.name, io.TextIOWrapper(stream, encoding="utf-8")
        return
    for path in sorted(location.rglob("radiated_data_*deg_short.csv")):
        with path.open(newline="", encoding="utf-8") as stream:
            yield str(path), stream


def parse_name(name: str):
    model = re.search(r"/(SF[1-5])_G1F1cmplt_QE95/", name.replace("\\", "/"))
    angle = re.search(r"radiated_data_(\d+\.\d)deg_short\.csv$", name)
    if not model or not angle:
        raise ValueError(f"unexpected table path: {name}")
    return model.group(1), float(angle.group(1))


def validate_table(name, stream):
    reader = csv.DictReader(stream, skipinitialspace=True)
    if not reader.fieldnames:
        raise ValueError(f"missing header: {name}")
    fields = {field.strip() for field in reader.fieldnames}
    missing = set(REQUIRED) - fields
    if missing:
        raise ValueError(f"missing {sorted(missing)}: {name}")

    valid_nu = []
    saw_valid = False
    saw_invalid = False
    previous_nu = None
    row_count = 0
    for row in reader:
        row_count += 1
        values = {key.strip(): value for key, value in row.items()}
        ebeam = float(values["Ebeam"])
        nu = float(values["nu"])
        born = float(values["XSborn_unp"])
        radiated = float(values["XSrad_unp"])
        if abs(ebeam - 10380.0) > 0.1:
            raise ValueError(f"unexpected Ebeam={ebeam} MeV: {name}")
        valid = all(math.isfinite(value) and value > 0.0
                    for value in (born, radiated))
        if valid:
            if saw_invalid:
                raise ValueError(f"non-contiguous finite ratio domain: {name}")
            if previous_nu is not None:
                dnu = nu - previous_nu
                if dnu <= 0.0:
                    raise ValueError(f"non-monotonic nu grid: {name}")
                if abs(dnu - RAD_DNU_MEV) > RAD_DNU_TOL:
                    raise ValueError(
                        f"nonuniform finite nu grid (dnu={dnu} MeV): {name}"
                    )
            valid_nu.append(nu)
            previous_nu = nu
            saw_valid = True
        elif saw_valid:
            saw_invalid = True
    if len(valid_nu) < 2:
        raise ValueError(f"fewer than two finite rows: {name}")
    if valid_nu[-1] != 9910.0:
        raise ValueError(f"finite upper nu is not 9910 MeV: {name}")
    return row_count, valid_nu[0], valid_nu[-1], len(valid_nu), RAD_DNU_MEV


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "location", nargs="?",
        default=Path("src/interp") / ARCHIVE_NAME,
        type=Path,
        help="archive or extracted Newfit_20260710_fullxquad_15angles directory",
    )
    args = parser.parse_args()
    if not args.location.exists():
        parser.error(f"not found: {args.location}")

    tables = {}
    for name, stream in csv_sources(args.location):
        try:
            model, angle = parse_name(name)
            if (model, angle) in tables:
                raise ValueError(f"duplicate table for {(model, angle)}: {name}")
            tables[model, angle] = validate_table(name, stream)
        finally:
            stream.close()

    expected = {(model, angle) for model in MODELS for angle in ANGLES}
    missing = expected - set(tables)
    extra = set(tables) - expected
    if missing or extra:
        raise ValueError(f"missing={sorted(missing)}, extra={sorted(extra)}")

    by_angle = defaultdict(list)
    for (model, angle), summary in tables.items():
        by_angle[angle].append(summary)
    for angle, summaries in sorted(by_angle.items()):
        if len(set(summaries)) != 1:
            raise ValueError(f"SF grids differ at {angle:.1f} degrees")
        rows, nu_min, nu_max, valid_rows, dnu = summaries[0]
        print(f"{angle:4.1f} deg: nu_min={nu_min:.0f}, nu_max={nu_max:.0f} MeV, "
              f"N={valid_rows}, dnu={dnu:.0f} MeV, rows={rows}")
    print("PASS: 75 tables validated; finite ratios contiguous and 1 MeV spaced")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, tarfile.TarError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)
