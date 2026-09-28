"""Write public browser settings; secret and legacy service-role keys are rejected."""
import argparse
import json
import os
from pathlib import Path
from urllib.parse import urlparse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    url = os.environ.get("SUPABASE_URL", "").rstrip("/")
    key = os.environ.get("SUPABASE_PUBLISHABLE_KEY", "")
    if url and key:
        parsed = urlparse(url)
        if parsed.scheme != "https" or not parsed.hostname or parsed.username or parsed.query or parsed.fragment:
            raise SystemExit("SUPABASE_URL must be an HTTPS project URL.")
        if not key.startswith("sb_publishable_"):
            raise SystemExit("Only a Supabase publishable key may be included in the public gallery.")
        settings = {"supabaseUrl": url, "publishableKey": key,
                    "turnstileSiteKey": os.environ.get("TURNSTILE_SITE_KEY", "")}
        (args.directory / "config.js").write_text("window.SUNCET_VOTING = " + json.dumps(settings, indent=2) + ";\n")
    elif url or key:
        raise SystemExit("Provide both SUPABASE_URL and SUPABASE_PUBLISHABLE_KEY.")
    elif not (args.directory / "config.js").exists():
        raise SystemExit("No public voting configuration is available.")


if __name__ == "__main__":
    main()
