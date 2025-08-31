"""
Configuration file for the Sphinx documentation builder.
"""

# pylint: disable=wrong-import-position, invalid-name

import os

import tomli
from docutils.parsers.null import Parser as NullParser
from sphinx.application import Sphinx


def setup(app: Sphinx):
    app.add_source_parser(NullParser)


THIS_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(THIS_DIR)

# -- Project information -----------------------------------------------------


project = "art_modern"
author = "YU Zhejian"
copyright_string = f"2024-2025, {author}"  # FIXME
release = "1.1.5"  # FIXME

# -- General configuration ---------------------------------------------------
html_theme = "furo"
extensions = [
    "sphinx.ext.todo",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "myst_parser",
    "sphinx_copybutton",
    # "sphinx_design",
]

myst_enable_extensions = ["deflist", "dollarmath"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", ".virtualenv/**"]

# html_static_path = ['_static']

# Source code suffixes
source_suffix = {".rst": "restructuredtext", ".md": "myst"}


with open(os.path.join(THIS_DIR, "latex_preamble.tex"), "r+") as r:
    PREAMBLE = r.read()

latex_elements = {
    # Additional stuff for the LaTeX preamble.
    "preamble": PREAMBLE,
}
