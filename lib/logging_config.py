"""
Standalone logging configuration for the MOA Catalog HTTP client.
"""

import os
from pathlib import Path

# Default log file lives under lib/logs/, overridable via MOA_CLIENT_LOG
# (e.g. to point somewhere else when running from a notebook or CI).
_DEFAULT_LOG_PATH = Path(__file__).parent / "logs" / "client.log"

DEBUG_LOG = os.environ.get("MOA_CLIENT_LOG", str(_DEFAULT_LOG_PATH))

# Ensure the log directory exists before any FileHandler opens it.
Path(DEBUG_LOG).parent.mkdir(parents=True, exist_ok=True)
