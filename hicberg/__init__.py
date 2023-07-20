# import sys
# import logging
# import logging.config

# logging.basicConfig(
#     level=logging.INFO,
#     format="%(asctime)s [%(levelname)s] %(message)s",
#     handlers=[
#         logging.FileHandler("hicberg.log"),
#         logging.StreamHandler()
#         ], stream=sys.stdout
#     )



# loggers = logging.getLogger("__name__")
# loggers.setLevel(logging.DEBUG)


# # Create the Handler for logging data to a file
# logger_handler = logging.FileHandler(filename='hicberg.log')
# logger_handler.setLevel(logging.INFO)

# loggers.

# # Create a Formatter for formatting the log messages
# logger_formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')

# # Add the Formatter to the Handler
# logger_handler.setFormatter(logger_formatter)

# # Add the Handler to the Logger
# loggers.addHandler(logger_handler)
# # loggers.info('Completed configuring logger()!') 


# logger = logging.getLogger(__name__)
# logger.setLevel(logging.DEBUG)
# # terminal logger
# stream_handler = logging.StreamHandler(sys.stdout)
# stream_handler.setLevel(logging.CRITICAL)
# logger.addHandler(stream_handler)
# # file logger
# file_handler = logging.FileHandler(f'{__name__}.log', 'w')
# file_handler.setLevel(logging.DEBUG)
# logger.addHandler(file_handler)


# logger = logging.getLogger(__name__)

# log_format = logging.Formatter('%(asctime)s %(levelname)s %(message)s')

# console_handler = logging.StreamHandler(sys.stdout)
# console_handler.setLevel(logging.INFO)
# console_handler.setFormatter(log_format)

# logger.addHandler(console_handler)

# file_handler = logging.FileHandler('hicberg.log', 'a')
# file_handler.setLevel(logging.INFO)
# file_handler.setFormatter(log_format)

# logger.addHandler(file_handler)
# logger.setLevel(logging.INFO)

# logging.basicConfig(
#     level=logging.INFO,
#     format="%(asctime)s [%(levelname)s] %(message)s",
#     handlers=[
#         logging.FileHandler("debug.log"),
#         logging.StreamHandler()
#     ]
# )



