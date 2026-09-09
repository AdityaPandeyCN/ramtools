# Sphinx configuration for the RAMTools documentation.
# Build locally with:  python3 -m sphinx -b html docs/source build/docs/html

project = "RAMTools"
copyright = "2025, compiler-research"
author = "compiler-research"

extensions = ["sphinx.ext.todo"]
templates_path = ["_templates"]
exclude_patterns = []

html_theme = "alabaster"
html_static_path = ["_static"]
html_theme_options = {
    "github_user": "compiler-research",
    "github_repo": "ramtools",
    "github_banner": True,
    "fixed_sidebar": True,
    "description": "Genomic alignments in ROOT's RNTuple format",
}

highlight_language = "bash"
