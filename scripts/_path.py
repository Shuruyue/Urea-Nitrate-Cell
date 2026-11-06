"""Let scripts/ run without installing the package (adds the project root to sys.path); print in UTF-8."""
import os
import sys

# When output is redirected to a file or pipe, Windows uses the locale encoding (cp950 on a
# Traditional-Chinese system), which cannot encode characters such as µ and −
for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        _stream.reconfigure(encoding="utf-8")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
OUTPUTS = os.path.join(ROOT, "outputs")
