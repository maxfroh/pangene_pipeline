import logging
import re

from colorama import Fore, Back, Style


SUCCESS_LEVEL = 21

logging.addLevelName(SUCCESS_LEVEL, f"{Fore.GREEN}SUCCESS{Fore.RESET}")


class CustomLogger(logging.getLoggerClass()):
    def debug(self, message, *args, **kwargs):
        message = Style.DIM + message + Style.RESET_ALL
        super().debug(message, *args, **kwargs, stacklevel=2)

    def success(self, message, *args, **kwargs):
        message = re.sub(r"success\w*", lambda m: f"{Fore.GREEN}{m.group(0)}{Fore.RESET}", message, flags=re.IGNORECASE)
        if self.isEnabledFor(SUCCESS_LEVEL):
            self._log(SUCCESS_LEVEL, message, args, **kwargs, stacklevel=2)
        else:
            self._log(logging.INFO, message, args, **kwargs, stacklevel=2)


logging.setLoggerClass(CustomLogger)
