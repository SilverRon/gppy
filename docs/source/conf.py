# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'gpPy'
copyright = '2024, Gregory S.H. Paek'
author = 'Gregory S.H. Paek'
release = '0.9'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# extensions = []

templates_path = ['_templates']
exclude_patterns = []

html_baseurl = 'https://silverron.github.io/gppy/'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'alabaster'
extensions = ['sphinx.ext.autodoc']
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_favicon = 'favicon.ico'

import os
import sys
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('../../config'))
sys.path.insert(0, os.path.abspath('../../phot'))
sys.path.insert(0, os.path.abspath('../../preprocess'))
sys.path.insert(0, os.path.abspath('../../util'))
