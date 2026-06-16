"""Run the container platform CLI when invoked as a module or package path."""

import sys

from container_manager.src.version import REPO_ROOT

try:
    from container_manager.src.cli import main
except ModuleNotFoundError:
    # Support direct invocation: python container_manager/src/
    sys.path.insert(0, str(REPO_ROOT))
    from container_manager.src.cli import main


if __name__ == "__main__":
    main()
