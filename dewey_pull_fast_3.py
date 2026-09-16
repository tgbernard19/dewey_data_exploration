#!/usr/bin/env python3
"""
Dewey Data bulk puller.

Pulls three datasets -- work_visits, other_visits, home_visits -- for each
day in Apr 11-14 2025 (inclusive), into per-day subfolders, with hardened
retry logic. The Dewey CDN periodically returns HTTP 500s, so each partition
is fetched over multiple passes with backoff, then swept once more for files
that download but come back corrupt/unreadable.

Layout produced:
    <OUT_ROOT>\\work_visits\\2025-04-11\\part_0.parquet ... part_N.parquet
    <OUT_ROOT>\\work_visits\\2025-04-11\\_urls.txt        (manifest of signed URLs)
    <OUT_ROOT>\\other_visits\\2025-04-11\\...
    <OUT_ROOT>\\home_visits\\2025-04-11\\...
    ... one folder per day, per dataset.

Re-running is safe and resumable: files already present and readable are
skipped, only broken/missing parts are re-fetched.

Speed notes (vs. the ThreadPool-of-curl-processes version):
  * Each pass is ONE `curl --parallel` process -> connection reuse + HTTP/2
    multiplexing, no per-file process/TLS-handshake overhead.
  * Readability is cached: a file that once read cleanly is never re-checked.
  * Validation reads only the parquet footer (parquet_schema), which is both
    cheaper and the right signal for truncated downloads.

Target platform: Windows. Uses the bundled curl.exe (Win10 1803+, i.e.
curl >= 7.66 with --parallel). Paths are handled Windows-aware throughout.

Usage:
    python dewey_pull_fast_3.py

Requires: deweypy, duckdb (readability check), and curl (>=7.66) on PATH.
"""

import os
import sys
import time
import subprocess
from datetime import date, timedelta

from deweypy.auth import set_api_key
from deweypy.download.synchronous import get_dataset_files

# --------------------------------------------------------------------------
# CONFIG -- fill these in (same as you're swapping the API key)
# --------------------------------------------------------------------------
API_KEY = "PUT_YOUR_API_KEY_HERE"

# The datasets. Each is a different Dewey dataset/folder, so each gets its
# own ID.
DATASETS = {
    "work_visits":  "prj_xo9czjhu__fldr_gfv4qahxiwsd4dwy",
    "other_visits": "prj_xo9czjhu__fldr_8zme9bwbekydvezq",
    "home_visits":  "prj_xo9czjhu__fldr_d7cqgtcj3nyi4usp",
}

# Where everything lands. Point this at a real, non-synced local drive.
OUT_ROOT = r"E:\dewey-apr2025"   # use the actual letter DATA mounted as

# Inclusive date range: Apr 11-14 2025.
START_DATE = date(2025, 4, 11)
END_DATE   = date(2025, 4, 14)

# Retry / validation knobs.
MIN_BYTES        = 1024   # anything smaller is treated as a failed download
MAX_SIZE_PASSES  = 12     # passes for the size-based download loop
MAX_READ_PASSES  = 6      # extra passes to re-fetch corrupt-but-present files
CURL_TIMEOUT     = 300    # per-file curl max time (seconds)
NUM_WORKERS      = 8      # parallel transfers within one curl process

# --------------------------------------------------------------------------


def daterange(start, end):
    """Inclusive range of dates."""
    d = start
    while d <= end:
        yield d
        d += timedelta(days=1)


def partition_key(day):
    """
    Dewey partition key: zero-padded ISO 8601 (YYYY-MM-DD).

    This MUST match the stored partition column format exactly. String
    partition keys are compared lexically, so a non-padded '2025-4-11' sorts
    GREATER than '2025-04-11' (at char 5, '4' > '0') and would silently drop
    the partition from an `after` bound. isoformat() gives the padded form.
    """
    return day.isoformat()


def part_path(out_dir, i):
    return os.path.join(out_dir, f"part_{i}.parquet")


def duckdb_path(path):
    """Forward-slash form of a path for embedding in a DuckDB SQL string."""
    return path.replace("\\", "/")


def size_broken(path):
    """True if the file is missing or suspiciously small."""
    return (not os.path.exists(path)) or os.path.getsize(path) < MIN_BYTES


def download_batch(urls, out_dir, todo):
    """
    Download the given part indices in ONE curl process.

    A single `curl --parallel` invocation reuses connections to the origin and
    multiplexes transfers, instead of spawning one process (and one TLS
    handshake) per file. The url/output pairs are fed via a config file on
    stdin because 256 pairs would overflow the Windows command-line length
    limit (~32K chars). Output paths are written with forward slashes: curl's
    config parser treats backslashes in quoted values as escape sequences and
    would corrupt Windows paths (E:\\dewey... -> E:dewey...).
    """
    cfg = []
    for i in todo:
        dest = part_path(out_dir, i).replace("\\", "/")
        cfg.append(f'url = "{urls[i]}"')
        cfg.append(f'output = "{dest}"')
    subprocess.run(
        ["curl", "--parallel", "--parallel-max", str(NUM_WORKERS),
         "-L", "--fail", "--retry", "3", "--retry-delay", "2",
         "--max-time", str(CURL_TIMEOUT), "--config", "-"],
        input="\n".join(cfg).encode(), capture_output=True,
    )


# ---- signed-URL resolution ------------------------------------------------

def urls_for(data_id, pk):
    """Resolve fresh signed URLs for one dataset/partition."""
    return get_dataset_files(
        data_id,
        partition_key_after=pk,
        partition_key_before=pk,
        to_list=True,
    )


# ---- readability check (catches files that download but are corrupt) ------

def _load_duckdb():
    try:
        import duckdb
        return duckdb.connect()
    except Exception as e:  # noqa: BLE001
        print(f"  [warn] duckdb unavailable ({e}); skipping readability check.")
        return None


def unreadable_parts(con, out_dir, n, verified):
    """
    Indices of present files that DuckDB cannot read (corrupt/stub).

    `verified` is a set of indices already confirmed readable; those are
    skipped (a valid parquet does not spontaneously rot), and newly confirmed
    files are added to it. Validation reads only the file footer via
    parquet_schema -- no row-group decompression -- which is cheaper and is
    the right signal for a truncated download (its footer is missing).
    """
    if con is None:
        return []
    bad = []
    for i in range(n):
        if i in verified:
            continue
        p = part_path(out_dir, i)
        if size_broken(p):
            continue  # handled by the size loop
        try:
            con.sql(f"SELECT * FROM parquet_schema('{duckdb_path(p)}') LIMIT 1").fetchall()
            verified.add(i)
        except Exception:  # noqa: BLE001
            bad.append(i)
    return bad


# ---- per-day driver -------------------------------------------------------

def download_day(name, data_id, day, con):
    """Fetch one dataset's partition for one day. Returns (n, n_missing)."""
    pk = partition_key(day)
    out_dir = os.path.join(OUT_ROOT, name, day.isoformat())
    os.makedirs(out_dir, exist_ok=True)
    verified = set()  # part indices confirmed readable this run (cache)

    print(f"\n=== {name} | {day.isoformat()} (partition_key={pk}) ===")

    # 1. Resolve signed URLs. Fail LOUDLY if this comes back empty.
    try:
        urls = urls_for(data_id, pk)
    except Exception as e:  # noqa: BLE001
        raise RuntimeError(
            f"get_dataset_files FAILED for {name} {day.isoformat()} "
            f"(data_id={data_id}): {e}"
        ) from e

    if not urls:
        raise RuntimeError(
            f"NO URLS returned for {name} {day.isoformat()} "
            f"(data_id={data_id}, partition_key={pk}). "
            f"Check the data_id, date, partition-key FORMAT, and that your "
            f"API key has access."
        )

    n = len(urls)
    print(f"  {n} URLs -> {out_dir}  ({NUM_WORKERS} parallel transfers)")

    # Manifest, so the exact signed URLs used are recoverable.
    with open(os.path.join(out_dir, "_urls.txt"), "w") as f:
        f.write("\n".join(urls) + "\n")

    # 2. Size-based download passes with backoff.
    #    Signed URLs expire (typically 15-60 min). A big national partition can
    #    take longer than the TTL to grind through, after which every curl 403s
    #    and the loop burns its budget on dead links. So: if a pass makes ZERO
    #    progress, assume expiry and re-resolve fresh URLs before continuing.
    prev_todo = None
    for attempt in range(MAX_SIZE_PASSES):
        todo = [i for i in range(n) if size_broken(part_path(out_dir, i))]
        if not todo:
            print(f"  all {n} parts present after {attempt} pass(es).")
            break
        if todo == prev_todo:
            print(f"  no progress last pass -> re-resolving signed URLs "
                  f"(likely expired).")
            new_urls = urls_for(data_id, pk)
            if len(new_urls) != n:
                # Ordering/count changed under us; safest to abort this day
                # rather than map URLs to the wrong part indices.
                raise RuntimeError(
                    f"URL count changed on re-resolve for {name} "
                    f"{day.isoformat()}: was {n}, now {len(new_urls)}. "
                    f"Aborting day to avoid index/URL mismatch."
                )
            urls = new_urls
            with open(os.path.join(out_dir, "_urls.txt"), "w") as f:
                f.write("\n".join(urls) + "\n")
        prev_todo = todo
        print(f"  pass {attempt + 1}: fetching {len(todo)} part(s)")
        download_batch(urls, out_dir, todo)
        time.sleep(2 * (attempt + 1))  # back off more each round
    else:
        still = [i for i in range(n) if size_broken(part_path(out_dir, i))]
        print(f"  [warn] size loop gave up; {len(still)} still missing: {still}")

    # 3. Readability sweep: re-fetch files that are present but corrupt.
    #    Same expiry logic: refresh URLs if a sweep makes no progress.
    prev_bad = None
    for attempt in range(MAX_READ_PASSES):
        bad = [i for i in range(n) if size_broken(part_path(out_dir, i))]
        bad += unreadable_parts(con, out_dir, n, verified)
        bad = sorted(set(bad))
        if not bad:
            break
        if bad == prev_bad:
            print(f"  read-sweep no progress -> re-resolving signed URLs.")
            new_urls = urls_for(data_id, pk)
            if len(new_urls) != n:
                raise RuntimeError(
                    f"URL count changed on re-resolve for {name} "
                    f"{day.isoformat()}: was {n}, now {len(new_urls)}."
                )
            urls = new_urls
        prev_bad = bad
        print(f"  read-sweep {attempt + 1}: re-fetching {len(bad)} corrupt/missing: {bad}")
        download_batch(urls, out_dir, bad)
        time.sleep(2 * (attempt + 1))

    # 4. Final tally for this day. Cheap: unreadable_parts skips everything in
    #    `verified`, so this only stats files and checks any stragglers.
    missing = [i for i in range(n) if size_broken(part_path(out_dir, i))]
    missing += unreadable_parts(con, out_dir, n, verified)
    missing = sorted(set(missing))
    if missing:
        print(f"  [FAIL] {name} {day.isoformat()}: {len(missing)}/{n} bad -> {missing}")
    else:
        print(f"  [OK]   {name} {day.isoformat()}: {n}/{n} present & readable")
    return n, len(missing)


def main():
    # Guard against running with the placeholder key still in place.
    if API_KEY in ("", "filler") or "PUT_" in API_KEY:
        sys.exit(
            "Refusing to run: set a real API_KEY at the top of this script "
            "(still set to the 'filler' placeholder)."
        )
    if any("PUT_" in v for v in DATASETS.values()):
        sys.exit("Refusing to run: fill in all dataset IDs first.")

    set_api_key(API_KEY)
    os.makedirs(OUT_ROOT, exist_ok=True)
    con = _load_duckdb()

    days = list(daterange(START_DATE, END_DATE))
    print(f"Pulling {len(DATASETS)} dataset(s) x {len(days)} day(s) into {OUT_ROOT}")

    summary = {}  # (name, day) -> (n, n_missing)
    for name, data_id in DATASETS.items():
        for day in days:
            n, n_missing = download_day(name, data_id, day, con)
            summary[(name, day.isoformat())] = (n, n_missing)

    # ---- final report ----
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    total_bad = 0
    for (name, day), (n, n_missing) in summary.items():
        total_bad += n_missing
        flag = "OK  " if n_missing == 0 else "FAIL"
        print(f"  [{flag}] {name:13s} {day}  {n - n_missing:>4}/{n} good")
    if total_bad:
        print(f"\n{total_bad} part(s) still bad across all days. "
              f"Re-run the script to retry only those (it resumes).")
        sys.exit(1)
    print("\nAll partitions complete, present, and readable.")


if __name__ == "__main__":
    main()
