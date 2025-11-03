# Cotton2K Simulation Model

[![Rust](https://img.shields.io/badge/Rust-1.70%2B-blue)](https://www.rust-lang.org)
[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://www.python.org)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)

Cotton2K is a cotton simulation model specially adapted for irrigated cotton production in arid regions. It was originally written by [Prof. Avishalom Marani][marani] and has been rewritten in Rust with Python bindings.

## Features

- High-performance cotton growth simulation
- Native Rust implementation with Python bindings
- TOML/JSON configuration format
- CSV output format
- Cross-platform support

## Usage

### As a Rust Library

Add to your `Cargo.toml`:
```toml
[dependencies]
cotton2k = { git = "https://github.com/tcztzy/cotton2k" }
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

## Requirements

- Rust 1.70+ (for Rust usage)
- Python 3.10+ (for Python bindings)
- libclang (temporary, during transition from C++ to Rust)

## Building

Build Rust library:
```bash
cargo build --release
```

Build and install Python bindings:
```bash
cd bindings/python
maturin develop --release
```

## Roadmap

- [ ] Complete Rust migration
- [ ] GUI interface
- [ ] Improved documentation
- [ ] SWAP model integration

## License

This project is licensed under:
- [MIT License](LICENSE) for the Rust implementation
- Original Cotton2K by Marani remains under [GPLv2](https://plantscience.agri.huji.ac.il/avishalom-marani/cotton2k_source)

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.

[marani]: https://plantscience.agri.huji.ac.il/avishalom-marani
