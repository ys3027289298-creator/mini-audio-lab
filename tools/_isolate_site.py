"""Bootstrap for `python3 -I`: -I drops PYTHONPATH and user site-packages,
but the script's own directory stays on sys.path. This module lives next to
the analysis scripts and explicitly re-adds the interpreter's user
site-packages, without ever adding the current working directory."""
import os
import site
import sys

_user = site.getusersitepackages()
if os.path.isdir(_user) and _user not in sys.path:
    sys.path.append(_user)
