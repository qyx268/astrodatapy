# Licensed under a 3-clause BSD style license - see LICENSE.rst
from . import misc_tools
from . import number_density
from . import correlation
from . import clustering


def _installed_version():
    try:
        from importlib.metadata import version  # Python >= 3.8
        return version("astrodatapy")
    except ImportError:
        pass
    except Exception:
        return ""
    try:
        from pkg_resources import get_distribution  # Python < 3.8
        return get_distribution("astrodatapy").version
    except Exception:
        return ""


__version__ = _installed_version()
