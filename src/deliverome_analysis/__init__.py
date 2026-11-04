"""deliverome-analysis — scientist-driven analyses for the Deliverome project.

Common imports:

    from deliverome_analysis.plotting import setup_plotting_style, save_figure
    from deliverome_analysis.s3 import read_tsv, read_parquet, download, list_keys
    from deliverome_analysis.env import load_env
    from deliverome_analysis.surfaceome import surface_accessions, canonical_accessions
    from deliverome_analysis.seqtools import domesticate_gene, recode_gene, generate_barcode
    from deliverome_analysis import cellosaurus, orcs  # public-dataset readers (CVCL, BioGRID ORCS)

See ``docs/`` for the full guide.
"""

__version__ = "0.1.0"
