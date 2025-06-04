import time
import logging
import logging.config
import pysam

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Create handlers
console_handler = logging.StreamHandler()
file_handler = logging.FileHandler(f'hicberg_{time.strftime("%Y_%m_%d_%H_%M_%S")}.log')
console_handler.setLevel(logging.INFO)
file_handler.setLevel(logging.INFO)

# Create formatters and add it to handlers
console_format = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s', "%Y-%m-%d -- %H:%M:%S")
file_format = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s', "%Y-%m-%d --  %H:%M:%S")
console_handler.setFormatter(console_format)
file_handler.setFormatter(file_format)

# Add handlers to the logger
logger.addHandler(console_handler)
logger.addHandler(file_handler)
logger.propagate = False


save = pysam.set_verbosity(0)

