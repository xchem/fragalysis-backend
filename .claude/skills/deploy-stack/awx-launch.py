#!/usr/bin/env python3
"""Run the awxkit ``awx`` CLI with two runtime compatibility shims applied.

awxkit's console script is unmaintained against recent Pythons and crashes
before it ever talks to AWX. This launcher patches the two breakages
in-process, then hands straight over to ``awxkit.cli.run()`` with the same
arguments, so it can be used exactly like the ``awx`` command:

    python3 awx-launch.py job_templates launch "<template>" \
        --extra_vars "..." --monitor --wait

The shims (each guarded so it is a no-op when the underlying Python/awxkit is
already fine):

1. pkg_resources -- awxkit does ``import pkg_resources`` only to read its own
   version string. setuptools 81+ removed that API, so on a modern setuptools
   the import raises ModuleNotFoundError. We install a minimal stand-in backed
   by the stdlib ``importlib.metadata`` before awxkit is imported.

2. argparse -- awxkit's HelpfulArgumentParser overrides the private
   ``_parse_known_args(self, args, ns)``. Python 3.13 added a third positional
   argument (``intermixed``) to that call, so the two-argument override raises
   ``TypeError: ... takes 3 positional arguments but 4 were given``. We replace
   it with a version that accepts and forwards any extra positionals.
"""

import sys
import types


def _install_pkg_resources_shim():
    """Make ``import pkg_resources`` succeed without setuptools < 81."""
    try:
        import pkg_resources  # noqa: F401  (succeeds on setuptools < 81)

        return
    except ModuleNotFoundError:
        pass

    import importlib.metadata as metadata

    shim = types.ModuleType("pkg_resources")

    class _Distribution:
        def __init__(self, version):
            self.version = version

    def get_distribution(name):
        # awxkit only reads ``.version`` off the result.
        return _Distribution(metadata.version(name))

    shim.get_distribution = get_distribution
    sys.modules["pkg_resources"] = shim


def _patch_argument_parser():
    """Let HelpfulArgumentParser tolerate Python 3.13's extra argparse arg."""
    from argparse import ArgumentParser

    from awxkit.cli.utils import HelpfulArgumentParser

    def _parse_known_args(self, args, namespace, *rest):
        for help_flag in ("-h", "--help"):
            # Mirror awxkit's original behaviour: drop -h/--help so the CLI
            # prints usage info instead of argparse's terse help.
            if help_flag in args:
                args.remove(help_flag)
        return ArgumentParser._parse_known_args(self, args, namespace, *rest)

    HelpfulArgumentParser._parse_known_args = _parse_known_args


def main():
    _install_pkg_resources_shim()
    _patch_argument_parser()

    from awxkit.cli import run

    # Present as "awx" so any usage/error text reads naturally.
    sys.argv[0] = "awx"
    return run()


if __name__ == "__main__":
    sys.exit(main())
