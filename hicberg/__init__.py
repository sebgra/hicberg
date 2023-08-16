import logging
import logging.config
import pysam

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Create handlers
console_handler = logging.StreamHandler()
file_handler = logging.FileHandler('hicberg.log')
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

# Benchmark logger

benchmark_logger = logging.getLogger(__name__)
benchmark_logger.setLevel(logging.INFO)

# Create handlers
benchmark_console_handler = logging.StreamHandler()
benchmark_file_handler = logging.FileHandler('hicberg_benchmark.log')
benchmark_console_handler.setLevel(logging.INFO)
benchmark_file_handler.setLevel(logging.INFO)

# Create formatters and add it to handlers
console_format = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s', "%Y-%m-%d -- %H:%M:%S")
file_format = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s', "%Y-%m-%d --  %H:%M:%S")
benchmark_console_handler.setFormatter(console_format)
benchmark_file_handler.setFormatter(file_format)

# Add handlers to the logger
benchmark_logger.addHandler(benchmark_console_handler)
benchmark_logger.addHandler(benchmark_file_handler)
benchmark_logger.propagate = False

save = pysam.set_verbosity(0)

