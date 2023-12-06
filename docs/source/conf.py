# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import inspect
import sys
import os

import django
from django.utils.html import strip_tags
from django.utils.encoding import force_text

sys.path.insert(0, os.path.abspath('../..'))
os.environ['DJANGO_SETTINGS_MODULE'] = 'fragalysis.settings'
django.setup()


# -- Project information -----------------------------------------------------

project = 'Fragalysis-Backend'
project_copyright = '2020, Rachael Skyner'
author = 'Rachael Skyner'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.coverage', 'sphinx.ext.napoleon']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

master_doc = 'index'


# A function that extracts documentation from the Django model classes.
# Field documentation is extracted from the model's help_text or verbose_name
# along with automated docs for the type.
#
# Inspired by https://www.djangosnippets.org/snippets/2533/
def process_docstring(app, what, name, obj, options, lines):
    # Unused arguments
    del app, what, name, options

    # This causes import errors if left outside the function
    from django.db import models

    # Only look at objects that inherit from Django's base model class
    if inspect.isclass(obj) and issubclass(obj, models.Model):
        for field in obj._meta.get_fields():  # pylint: disable=protected-access
            # Try the column's help text or the verbose name
            text = ""
            if hasattr(field, "help_text"):
                text = strip_tags(force_text(field.help_text))
            if not text and hasattr(field, "verbose_name"):
                text = force_text(field.verbose_name).capitalize()
            if text:
                lines.append(f':param {field.attname}: {text}')
                lines.append(f':type {field.attname}: {type(field).__name__}')

    return lines


# Register the docstring processor with sphinx
def setup(app):
    app.connect('autodoc-process-docstring', process_docstring)
