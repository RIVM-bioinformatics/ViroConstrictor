import logging

from rich.console import Console
from snakemake_interface_logger_plugins.base import LogHandlerBase

import ViroConstrictor
import ViroConstrictor.logging
from snakemake_logger_plugin_viroconstrictor.handler import ViroConstrictorLogHandler


class LogHandler(LogHandlerBase, ViroConstrictorLogHandler):
    def __post_init__(self) -> None:
        snakemake_logger = logging.getLogger("snakemake")
        snakemake_logger.propagate = False
        snakemake_logger.setLevel(logging.CRITICAL)
        console = Console(
            stderr=False,
        )

        ViroConstrictorLogHandler.__init__(
            self, console=console, settings=self.common_settings
        )
        # self.setFormatter(RichFormatter(console, self.common_settings.printshellcmds))
        # self.addFilter(RichFilter())

    @property
    def writes_to_stream(self) -> bool:
        """
        Whether this plugin writes to stderr/stdout
        """
        return False

    @property
    def writes_to_file(self) -> bool:
        """
        Whether this plugin writes to a file
        """
        return False

    @property
    def has_filter(self) -> bool:
        """
        Whether this plugin attaches its own filter
        """
        return True

    @property
    def has_formatter(self) -> bool:
        """
        Whether this plugin attaches its own formatter
        """
        return True

    @property
    def needs_rulegraph(self) -> bool:
        """
        Whether this plugin requires the DAG rulegraph.
        """
        return True

    def close(self) -> None:
        """
        Close the log handler.
        """
        super().close()

        # Close the console
        self.console = None

        # De-initialize the ViroConstrictorBaseLogHandler
        ViroConstrictor.logging.ViroConstrictorBaseLogHandler.close(self)
