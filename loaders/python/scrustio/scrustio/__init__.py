"""scrustio: load scrust's outputs into AnnData objects.

Reads the plain files scRust's command-line tool (scrust) writes, through the manifest in each
output's summary.json (schema ``scrust-files``, version 1; docs/output_manifest.md in the scRust
repository), and builds AnnData objects. Each matrix goes only into a slot meant for its scale:
the raw counts into ``layers["counts"]``; SCTransform's corrected counts, their log1p and the
residuals only into slots for normalised data.
"""

from .anndata_loader import load_anndata, write_h5ad
from .manifest import SCHEMA, SCHEMA_VERSION, Manifest, ScrustError, read_manifest, verify_sha256
from .slots import SLOTS, check_slot, read_for_slot, slot_table

__version__ = "0.1.0"
TESTED = {"anndata": "0.13.4"}

__all__ = [
    "SCHEMA", "SCHEMA_VERSION", "SLOTS", "TESTED", "Manifest", "ScrustError", "check_slot", "load_anndata",
    "read_for_slot", "read_manifest", "slot_table", "verify_sha256", "write_h5ad",
]
