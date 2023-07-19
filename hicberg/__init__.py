import sys
import logging
import logging.config

# logging.basicConfig(
#     level=logging.INFO,
#     format="%(asctime)s [%(levelname)s] %(message)s",
#     handlers=[
#         logging.FileHandler("hicberg.log"),
#         logging.StreamHandler()
#         ], stream=sys.stdout
#     )



loggers = logging.getLogger("__name__")
loggers.setLevel(logging.DEBUG)

# Create the Handler for logging data to a file
logger_handler = logging.FileHandler(filename='hicberg.log')
logger_handler.setLevel(logging.INFO)

# Create a Formatter for formatting the log messages
logger_formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')

# Add the Formatter to the Handler
logger_handler.setFormatter(logger_formatter)

# Add the Handler to the Logger
loggers.addHandler(logger_handler)
# loggers.info('Completed configuring logger()!') 




