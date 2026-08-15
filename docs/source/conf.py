# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'RabbitSketch'
copyright = '2026, RabbitSketch contributors'
author = 'RabbitSketch contributors'
release = '2.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

try:
    import sphinx_rtd_theme  # noqa: F401
except ImportError:
    # Keep local/source-distribution documentation buildable with Sphinx's
    # bundled theme; docs/requirements.txt installs the preferred RTD theme.
    html_theme = 'alabaster'
else:
    html_theme = 'sphinx_rtd_theme'
html_static_path = []
