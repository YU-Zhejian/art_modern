"""

"""
import sys
from enum import IntEnum
from typing import Union


class FileType(IntEnum):
    CONST_COV = 0
    COV_TSV = 1
    PBSIM3_TRANSCRIPT = 2

def validate(in_fq_path: str, in_cov: Union[float, str], file_type: FileType) -> None:
    ...

if __name__ == "__main__":
    sys.exit(0)
