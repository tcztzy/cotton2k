import argparse
import pathlib

from cotton2k import run

parser = argparse.ArgumentParser()
parser.add_argument("profile", help="Profile file path")
parser.add_argument(
    "--output",
    action="store_true",
    help="Dump simulation result to user data directory.",
)
args = parser.parse_args()
profile = args.profile
if (path := pathlib.Path(profile)).exists():
    run(path, output=args.output)
else:
    raise NotImplementedError(
        f"File {path} cannot be found.\n"
        "Currently, cotton2k only support run simulation from profile file."
    )
