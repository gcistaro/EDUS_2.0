import os
import sys
sys.path.insert(0, os.path.abspath('.'))

project = "EDUS"
author = "Your Name"
release = "1.0"

extensions = [
    "myst_parser",
    "sphinx_rtd_theme",
    "sphinx_copybutton",
]

myst_enable_extensions = [
    "colon_fence",    # allows ::: blocks
    "dollarmath",     # allows $...$ and $$...$$
    "amsmath",        # equation environments
]

templates_path = ["_templates"]
html_static_path = ["_static"]

# Wannier90-like theme
html_theme = "sphinx_rtd_theme"
html_logo = "_static/logo.png"   # optional
html_favicon = "_static/favicon.ico"  # optional

html_theme_options = {
    "collapse_navigation": True,
    "navigation_depth": 4,
    "sticky_navigation": True,
    "titles_only": False,
}

# IMPORTANT: defines the homepage
master_doc = "index"
