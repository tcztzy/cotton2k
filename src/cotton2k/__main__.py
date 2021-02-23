import sys

from cotton2k import run

if len(sys.argv) != 2:
    sys.exit("Usage:\ncotton2k <profile name>\n")
run(sys.argv[1])
