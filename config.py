"""
Shared configuration for the Dewey Data pulling scripts.

Non-secret settings (dataset IDs, output paths, date ranges, retry knobs)
live here so they can be edited without touching script logic.

The Dewey API key is deliberately NOT stored in this file. Set it via the
DEWEY_API_KEY environment variable instead, e.g.:

    export DEWEY_API_KEY="akv1_..."          # macOS/Linux
    setx DEWEY_API_KEY "akv1_..."             # Windows (new shells)

Never commit real API keys to source control.
"""

import os
from datetime import date

# --------------------------------------------------------------------------
# Credentials
# --------------------------------------------------------------------------
# Read from the environment; left blank if unset so callers can decide how
# to handle a missing key (e.g. fail loudly with a helpful message).
API_KEY = os.environ.get("DEWEY_API_KEY", "")

# --------------------------------------------------------------------------
# Datasets
# --------------------------------------------------------------------------
# Each is a different Dewey dataset/folder, so each gets its own ID.
DATASETS = {
    "work_visits":  "prj_xo9czjhu__fldr_gfv4qahxiwsd4dwy",
    "other_visits": "prj_xo9czjhu__fldr_8zme9bwbekydvezq",
    "home_visits":  "prj_xo9czjhu__fldr_d7cqgtcj3nyi4usp",
}

# --------------------------------------------------------------------------
# Output location
# --------------------------------------------------------------------------
# Where everything lands. Point this at a real, non-synced local drive.
# Can be overridden with the DEWEY_OUT_ROOT environment variable.
OUT_ROOT = os.environ.get("DEWEY_OUT_ROOT", r"E:\dewey-apr2025")

# --------------------------------------------------------------------------
# Date range (inclusive)
# --------------------------------------------------------------------------
START_DATE = date(2025, 4, 11)
END_DATE = date(2025, 4, 14)

# --------------------------------------------------------------------------
# Retry / validation knobs
# --------------------------------------------------------------------------
MIN_BYTES = 1024          # anything smaller is treated as a failed download
MAX_SIZE_PASSES = 12      # passes for the size-based download loop
MAX_READ_PASSES = 6       # extra passes to re-fetch corrupt-but-present files
CURL_TIMEOUT = 300        # per-file curl max time (seconds)
NUM_WORKERS = 8           # parallel transfers within one curl process
