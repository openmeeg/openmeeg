from contextlib import contextmanager

from ._openmeeg_wrapper import DEBUG, ERROR, INFORMATION, WARNING, Logger

_warn_map = dict(
    error=ERROR,
    warning=WARNING,
    info=INFORMATION,
    debug=DEBUG,
)


def set_log_level(level):
    """Set the logging level.

    Parameters
    ----------
    level : str
        Can be ``'error'``, ``'warning'``, ``'info'``, or ``'debug'``.
    """
    if not isinstance(level, str) or level not in _warn_map:
        raise ValueError(
            f"Unknown level {repr(level)}, must be one of " f"{list(_warn_map)}"
        )
    Logger.logger().set_info_level(_warn_map[level])


def get_log_level():
    """Get the current logging level.

    Returns
    -------
    level : str
        See :func:`set_log_level`.
    """
    rev_map = {val: key for key, val in _warn_map.items()}
    return rev_map[Logger.logger().get_info_level()]


@contextmanager
def use_log_level(level):
    """Context manager for logging level.

    Parameters
    ----------
    level : str
        See :func:`set_log_level`.
    """
    old_level = get_log_level()
    set_log_level(level)
    try:
        yield
    finally:
        set_log_level(old_level)
