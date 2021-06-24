import argparse
import pathlib

from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from cotton2k import ENGINE_URL, run
from cotton2k.models import Profile

parser = argparse.ArgumentParser()
parser.add_argument("profile", help="Profile file path or name")
parser.add_argument(
    "--output",
    action="store_true",
    help="Dump simulation result to user data directory.",
)
args = parser.parse_args()
if (path := pathlib.Path(args.profile)).exists():
    run(path, output=args.output)
else:
    session = Session(bind=create_engine(ENGINE_URL))
    profile = session.query(Profile).where(Profile.name == args.profile).one_or_none()
    if profile is None:
        raise FileNotFoundError(f"Cannot found the profile {args.profile}")
    run(profile, output=args.output)
