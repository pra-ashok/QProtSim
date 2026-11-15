"""Allow strict documentation builds from source archives without Git history."""

import logging

from mkdocs.plugins import event_priority


timestamp_logger = logging.getLogger("mkdocs.plugins.mixins/git")


def _ignore_unused_timestamp_warning(record):
    return record.getMessage() != "mixins/git: No git commit found."


@event_priority(100)
def on_config(config):
    # Shadcn checks Git history even when its timestamp display is disabled.
    if not config.theme.get("show_datetime", False):
        timestamp_logger.addFilter(_ignore_unused_timestamp_warning)


def on_post_build(config):
    timestamp_logger.removeFilter(_ignore_unused_timestamp_warning)


def on_build_error(error):
    timestamp_logger.removeFilter(_ignore_unused_timestamp_warning)
