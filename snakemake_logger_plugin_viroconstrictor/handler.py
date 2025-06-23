from rich.console import Console
from snakemake_interface_logger_plugins.settings import OutputSettingsLoggerInterface

import ViroConstrictor
import ViroConstrictor.logging
from ViroConstrictor.logging import ViroConstrictorBaseLogHandler


class ViroConstrictorLogHandler(ViroConstrictorBaseLogHandler):
    """
    Custom log handler for ViroConstrictor, extends ViroConstrictorBaseLogHandler.
    All custom log handling is done in the ViroConstrictorBaseLogHandler, this is just the hook to initialize it as a snakemake logger plugin.

    Attributes
    ----------
    console : Console
        Rich Console instance for enhanced terminal output.
    printshellcmds : bool
        Whether to print shell commands.
    nocolor : bool
        Whether to disable colored output.
    quiet : bool
        Whether to suppress output.
    debug_dag : bool
        Whether to enable DAG debugging output.
    verbose : bool
        Whether to enable verbose output.
    show_failed_logs : bool
        Whether to show logs for failed rules.
    stdout : bool
        Whether to redirect output to stdout.
    dryrun : bool
        Whether it is a dry run.
    Methods
    -------
    close()
        Closes the log handler.
    """

    def __init__(
        self, console: Console, settings: OutputSettingsLoggerInterface, *args, **kwargs
    ) -> None:
        self.console = console

        super().__init__(logfile_path=ViroConstrictor.logging.logfile, *args, **kwargs)

        self.printshellcmds = settings.printshellcmds
        self.nocolor = settings.nocolor
        self.quiet = settings.quiet
        self.debug_dag = settings.debug_dag
        self.verbose = settings.verbose
        self.show_failed_logs = settings.show_failed_logs
        self.stdout = settings.stdout
        self.dryrun = settings.dryrun

    def close(self) -> None:
        super().close()
