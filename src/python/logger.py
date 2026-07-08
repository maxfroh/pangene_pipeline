import logging

SUCCESS_LEVEL = 21

logging.addLevelName(SUCCESS_LEVEL, "SUCCESS")


class ANSIColors:
    GRAY = "\033[90m"
    GREEN = "\033[92m"
    BOLD = "\033[1m"
    OFF = "\033[0m"


class CustomLogger(logging.getLoggerClass()):
    def debug(self, message, *args, **kwargs):
        message = ANSIColors.GRAY + message + ANSIColors.OFF
        super().debug(message, *args, **kwargs)

    def success(self, message, *args, **kwargs):
        if self.isEnabledFor(SUCCESS_LEVEL):
            self._log(SUCCESS_LEVEL, message, args, **kwargs, stacklevel=2)


logging.setLoggerClass(CustomLogger)
