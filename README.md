# Cotton2K simulation model

Cotton2K is a cotton simulation model specially adapted for irrigated cotton production in arid regions.
It was originally written by [Prof. Avishalom Marani][marani] and distributed under GPL, and rewrite in Rust by [Tang Ziya][tang].

These packages are distributed under [GPL 3.0 or later](https://www.gnu.org/licenses/gpl-3.0.en.html). This version has been upgraded to run on modern OS.

## Migrate from version 4.0

I use `TOML` as input file format and `CSV` as output file format. Documents coming soon.

## What about Python?

I wrote a [Python version](https://github.com/tcztzy/cotton2k-core) before this. It works, but I have serious performance issues, it toke about 2min one single simulation, compared with 3s in C++ version. I tried to improve, but it is too difficult for me.

I provided Python API via PyO3! see `bindings` directory in the root of this repository.

## Roadmap

1. GUI
2. Fix calculation errors
3. Update formulas, grant SWAP model
4. Refactoring all in Rust
4. Maybe I will change model structure

## Requirements

C++ compiler (clang and msvc tested), Rust toolchain (nightly tested), FORTRAN compiler and Meson build.

Additionally, you need to notice that [`LIBCLANG_PATH`](https://rust-lang.github.io/rust-bindgen/requirements.html#clang) should be satisfied. This dependency will be eventually removed after C++ codes are all refactored into Rust.

## Build

```
cargo build
```

## Licenses

* [Cotton2K by Marani: GPLv2](https://plantscience.agri.huji.ac.il/avishalom-marani/cotton2k_source)
<!-- * [GOSSYM: CC0](https://data.nal.usda.gov/dataset/gossym)
* [SWAP: GPLv2](https://www.swap.alterra.nl/DownloadRecent/swap4.0.1/Swap4.0.1.htm) -->
* [THIS REPO: GPLv3 or later](https://github.com/tcztzy/cotton2k)

## Contributing

If you have any problem, feel free to open an issue.

[marani]: https://plantscience.agri.huji.ac.il/avishalom-marani
[tang]: https://github.com/tcztzy
