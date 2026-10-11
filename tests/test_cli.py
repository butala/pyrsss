"""Tests for pyrsss.util.cli (console-script logging options)."""
import logging

from pyrsss.util.cli import add_log_options, configure_logging
from argparse import ArgumentParser


def _args(argv):
    parser = ArgumentParser()
    add_log_options(parser)
    return parser.parse_args(argv)


def test_default_is_info():
    assert configure_logging(_args([])) == logging.INFO


def test_quiet_and_verbose_levels():
    assert configure_logging(_args(['--quiet'])) == logging.WARNING
    assert configure_logging(_args(['--verbose'])) == logging.DEBUG
    assert configure_logging(_args(['-q'])) == logging.WARNING
    assert configure_logging(_args(['-v'])) == logging.DEBUG


def test_quiet_and_verbose_are_mutually_exclusive():
    import pytest

    with pytest.raises(SystemExit):
        _args(['--quiet', '--verbose'])


def test_short_flags_can_be_left_for_passthrough():
    import pytest

    parser = ArgumentParser()
    add_log_options(parser, short=False)
    args = parser.parse_args(['--quiet'])
    assert configure_logging(args) == logging.WARNING
    with pytest.raises(SystemExit):
        parser.parse_args(['-q'])
