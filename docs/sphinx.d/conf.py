"""
Configuration file for the Sphinx documentation builder.
"""

# pylint: disable=wrong-import-position, invalid-name

import os
import datetime
import sys

import jinja2
from sphinx import __version__ as sphinx_version
from myst_parser import __version__ as myst_parser_version

from docutils.parsers.null import Parser as NullParser
from sphinx.application import Sphinx


def setup(app: Sphinx):
    app.add_source_parser(NullParser)


THIS_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(THIS_DIR)

# -- Project information -----------------------------------------------------


project = "art_modern"
author = "YU Zhejian"
copyright_string = f"2024-{datetime.datetime.now().year}, {author}"
release = os.environ.get("PACKAGE_VERSION")

# -- General configuration ---------------------------------------------------
html_theme = "furo"
extensions = [
    "sphinx.ext.todo",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "myst_parser",
    "sphinx_copybutton",
    # Replaces "sphinx.ext.imgconverter",
    "sphinxcontrib.cairosvgconverter",
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

template_path = os.path.join(THIS_DIR, "versions.md.jinja2")
context = {
    "sphinx_version": sphinx_version,
    "python_version": sys.version,
    "myst_parser_version": myst_parser_version,
}

with open(template_path, encoding="UTF-8") as f:
    template = jinja2.Template(f.read())

with open(os.path.join(THIS_DIR, "src", "versions.md"), "w", encoding="UTF-8") as f:
    f.write(template.render(context))
