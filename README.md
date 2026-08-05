# Cotton2K Simulation Model

[![Rust](https://img.shields.io/badge/Rust-1.70%2B-blue)](https://www.rust-lang.org)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://www.python.org)
[![Crates.io](https://img.shields.io/crates/v/cotton2k)](https://crates.io/crates/cotton2k)
[![License](https://img.shields.io/badge/License-AGPL--3.0--or--later-green)](LICENSE)

Cotton2K is a cotton simulation model specially adapted for irrigated cotton production in arid regions. It was originally written by [Prof. Avishalom Marani][marani] and has been rewritten in Rust with Python bindings.

## Features

- High-performance cotton growth simulation
- Native Rust implementation with Python bindings
- TOML/JSON configuration format
- CSV output format
- Cross-platform support

## Workspace Layout

The Rust implementation is a Cargo workspace. The simulation engine crate
(`cotton2k`) coordinates these independent model crates:

- `cotton2k-core`: shared profile, state, legacy globals, and numerical helpers
- `cotton2k-atmosphere`: weather, radiation, humidity, wind, and evapotranspiration
- `cotton2k-soil`: soil hydrology, nitrogen, and temperature
- `cotton2k-plant`: plant growth, phenology, roots, abscission, and nitrogen

Build or test all workspace members with:

```bash
cargo test --workspace
```

## Usage

### As a Rust Library

Add to your `Cargo.toml`:
```toml
[dependencies]
cotton2k = "0.1.0"
```

Or install the command-line binaries from crates.io:
```bash
cargo install cotton2k
```

### As a Python Package

Install from PyPI:
```bash
pip install cotton2k
```

Or install locally:
```bash
cd bindings/python
pip install .
```

Example usage:
```python
import cotton2k as c2k

# Run simulation with profile file
c2k.run("path/to/profile.toml")
```

### As a Worker Process (for GUI/Batch)

The local worker binary runs one simulation job per process and emits JSONL progress events.

```bash
cargo run --bin cotton2k-worker -- \
  run \
  --profile /path/to/profile.toml \
  --run-dir /path/to/runs/job-001 \
  --run-id job-001
```

Artifacts in `run-dir`:
- `output.csv`
- `meta.json`
- `request.json`
- `input/profile.toml` or `input/profile.json`

### As a Batch Scheduler (parallel worker processes)

Prepare a jobs file (`jobs.json`):

```json
[
  { "profile": "/path/to/case-a/profile.toml", "run_id": "case-a" },
  { "profile": "/path/to/case-b/profile.toml", "run_id": "case-b" }
]
```

Run in parallel (default max parallelism is `min(4, physical_cpu_cores)`):

```bash
cargo run --bin cotton2k-batch -- \
  run \
  --jobs /path/to/jobs.json \
  --runs-root /path/to/runs \
  --max-parallel 2
```

## Requirements

- Rust 1.70+ (for Rust usage)
- Python 3.10+ (for Python bindings)
- No libclang/LLVM toolchain is required for Python wheel builds

## Building

Build Rust library:
```bash
cargo build --release
```

Build worker binary:
```bash
cargo build --release --bin cotton2k-worker
```

Build batch binary:
```bash
cargo build --release --bin cotton2k-batch
```

Build and install Python bindings:
```bash
cd bindings/python
maturin develop --release
```

Run the pure-Rust regression check suite:
```bash
./scripts/regression_pure_rust.sh
```

Run the repository checks manually:
```bash
pre-commit run --all-files
```

Install the checks as a Git pre-commit hook:
```bash
pre-commit install
```

### Build with Nix

This repository now provides a flake with:
- `packages.<system>.default` / `packages.<system>.cotton2k`
- `apps.<system>.default`
- `devShells.<system>.default`

Build the package:
```bash
nix build .#cotton2k
```

Run directly:
```bash
nix run .#cotton2k -- /path/to/profile.toml
```

Enter the development shell:
```bash
nix develop
```

## Roadmap

- [ ] Complete Rust migration
- [ ] GUI interface
- [ ] Improved documentation
- [ ] SWAP model integration

## License

This project is licensed under:
- [AGPLv3+](LICENSE) for the Rust implementation
- Original Cotton2K by Marani remains under [GPLv2+](https://plantscience.agri.huji.ac.il/avishalom-marani/cotton2k_source)

### Correction
I used to think this project could be licensed under MIT license, but finally I found that my translate-to-Rust implementation should also follow GPL, so I relicensed all code in Rust & Python under AGPLv3+.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.

[marani]: https://plantscience.agri.huji.ac.il/avishalom-marani
