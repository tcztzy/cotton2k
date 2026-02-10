#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

echo "[1/5] cargo fmt --all --check"
cargo fmt --all --check

echo "[2/5] cargo check --all-targets"
cargo check --all-targets

echo "[3/5] cargo test --all-targets"
cargo test --all-targets

echo "[4/5] verify legacy C++ artifacts are removed"
for legacy_file in build.rs global.h CottonSimulation.h GeneralFunctions.h; do
    if [[ -e "$legacy_file" ]]; then
        echo "Found unexpected legacy file: $legacy_file" >&2
        exit 1
    fi
done

echo "[5/5] verify stale migration comments are removed"
legacy_hits="$(rg -n "(\\.cpp|global\\.h|CottonSimulation\\.h|GeneralFunctions\\.h|transition from C\\+\\+ to Rust|translation to C\\+\\+)" src README.md Cargo.toml || true)"
if [[ -n "$legacy_hits" ]]; then
    echo "Found stale migration references:" >&2
    echo "$legacy_hits" >&2
    exit 1
fi

echo "Pure-Rust regression checks passed."
