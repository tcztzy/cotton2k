#!/usr/bin/env python
import os
import sys

from argparse import ArgumentParser
from datetime import date
from difflib import context_diff
from itertools import product
from pathlib import Path
import re

today = date.today().strftime("%A, %B %d, %Y").ljust(30)
extensions = ["B01", "F01", "PLM", "PLT", "S01", "SMP", "WA2", "WAT", "LWP", "CHB", "RUT", "TM1", "TM2", "TMS", "NB0", "NB1", "NB2", "NB3", "NB4"]
fixtures = Path(__file__).parent / "fixtures"
parser = ArgumentParser()
parser.add_argument("output_dir", type=Path, default=Path('.'))
parser.add_argument("--profiles", metavar="P", nargs="+", default=["hakb1", "hava98"], help="Profiles to be compared")

def replace_zero(lines: list[str]):
    result = []
    for line in lines:
        if re.match('^(  [1-9]| [1-3]\d| 40)', line):
            line = line[:3] + line[3:].replace('0', ' ', -1)
        result.append(line)
    return result

def main(output_dir: Path, profiles: list[str]):
    for pro, ext in product(profiles, extensions):
        output = output_dir / f"{pro}.{ext}"
        fixture = fixtures / f"{pro}.{ext}"

        if output.exists():
            if not fixture.exists():
                sys.stderr.write(f"Fixture {fixture} doesn't exist.{os.linesep}")
                continue
            print(f"Diff {output} and {fixture}")
            a = output.read_text().splitlines(True)
            b = fixture.read_text().format(today).splitlines(True)
            if ext == "SMP":
                a = replace_zero(a)
                b = replace_zero(b)
            sys.stdout.writelines(context_diff(a, b))


if __name__ == "__main__":
    main(**vars(parser.parse_args()))
