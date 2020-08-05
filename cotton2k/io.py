import os
from locale import atof, atoi
from pathlib import Path
from datetime import datetime

from .utils import get_line_data


def parse_profile(content):
    lines = content.splitlines()
    return {
        'description': lines[0][20:].strip(),
        'dateEmerge': datetime.strptime(lines[1][:11], "%d-%b-%Y").date(),
        'dateSimStart': datetime.strptime(lines[1][15:26], "%d-%b-%Y").date(),
        'dateSimEnd': datetime.strptime(lines[1][30:41], "%d-%b-%Y").date(),
    }
