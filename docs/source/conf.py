# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import os
import sys
sys.path.insert(0, os.path.abspath('../../'))
# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'MetaNorm'
copyright = '2023, Jackson Barth'
author = 'Jackson Barth, Yuqiu Yang, Guanghua Xiao, Xinlei Wang'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx_rtd_theme',
    # 'sphinx.ext.autodoc',
    # View source code 
    # 'sphinx.ext.viewcode',
    # Check building time
    # 'sphinx.ext.duration',
    # Generate CLI documents 
    # 'sphinxcontrib.autoprogram',
    #'sphinx.ext.doctest',
    #'sphinx.ext.coverage',
    # 'sphinx.ext.napoleon',
    # Generate github history
    # 'sphinx_git',
]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_logo = "../../assets/logo.png"
html_theme_options = {
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': True,
    'logo_only': True,
    'display_version': False
}
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']