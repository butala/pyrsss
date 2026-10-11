"""Shared CLI plumbing for the console scripts.

The scripts' loggers only say something if logging is configured on the
entry point -- ``logging.basicConfig`` under ``if __name__ == '__main__'``
never runs when a console script calls :func:`main` directly. Every CLI
therefore calls :func:`add_log_options` on its parser and
:func:`configure_logging` on the parsed arguments:

    parser = ArgumentParser(...)
    ...
    add_log_options(parser)
    args = parser.parse_args(argv[1:])
    configure_logging(args)

Default is INFO (progress is visible); ``--quiet`` keeps only warnings
and errors; ``--verbose`` turns on debug. Use ``short=False`` for CLIs
that forward unknown options to another program (so ``-q``/``-v`` stay
available to it).
"""

import logging

FORMAT = '%(asctime)s %(levelname)s %(name)s: %(message)s'


def add_log_options(parser, short=True):
    """Add ``--quiet``/``--verbose`` (and ``-q``/``-v`` when *short*)."""
    group = parser.add_mutually_exclusive_group()
    quiet = ['--quiet']
    verbose = ['--verbose']
    if short:
        quiet.insert(0, '-q')
        verbose.insert(0, '-v')
    group.add_argument(*quiet, action='store_true',
                       help='only warnings and errors')
    group.add_argument(*verbose, action='store_true',
                       help='debug logging')


def configure_logging(args):
    """
    Configure logging from the options of :func:`add_log_options`
    (missing attributes fall back to INFO).
    """
    if getattr(args, 'quiet', False):
        level = logging.WARNING
    elif getattr(args, 'verbose', False):
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(level=level, format=FORMAT, force=True)
    return level
