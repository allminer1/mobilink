"""linkstat: chance-corrected, read-length-robust co-localization of gene classes on long reads and contigs."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("linkstat")
except PackageNotFoundError:  # running from a source tree without installation
    __version__ = "0.0.0+unknown"

__all__ = ["__version__"]
