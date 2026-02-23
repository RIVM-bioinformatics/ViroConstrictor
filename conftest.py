import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parents[1]
root_str = str(repo_root)
if root_str not in sys.path:
    sys.path.insert(0, root_str)