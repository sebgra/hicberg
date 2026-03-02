from abc import ABC
from typing import List
import subprocess
import logging

logger = logging.getLogger(__name__)

class BaseCommand(ABC):
    def __init__(self, cpus: int = 1, dry_run: bool = False) -> None:
        self.cpus : int = cpus
        self. dry_run : bool = dry_run
        
    def run(self, command: List[str], **kwargs)-> None:
        """
        TO BE completed
        """
        cmd_str: str = " ".join(map(str, command))
        if self.dry_run:
            logger.info(f"[DRY RUN]: {cmd_str}")
            return None
        logger.debug(f"Executing command: {cmd_str}")
        return subprocess.run(command, check = True, capture_output = True, text = True, **kwargs)
        