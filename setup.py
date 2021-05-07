import os
import subprocess
from collections import defaultdict
from glob import glob

try:
    from multiprocessing import cpu_count

    from Cython.Build import build_ext, cythonize

    extensions = cythonize("src/_cotton2k/*.pyx", nthreads=cpu_count())
except ImportError:
    from pathlib import Path

    from setuptools import Extension
    from setuptools.command.build_ext import build_ext

    extensions = list(
        map(
            lambda source: Extension(
                f"{source.parent.name}.{source.stem}",
                sources=[str(source).replace(".pyx", ".cpp")],
            ),
            Path("src/_cotton2k").glob("*.pyx"),
        )
    )
from pyproject_toml import setup

extra_compile_args = defaultdict(lambda: ["-std=c++20"])
extra_compile_args["msvc"] = ["/std:c++latest"]
libraries = defaultdict(lambda: ["cotton2k"])
libraries["nt"].extend(["ws2_32", "userenv", "advapi32"])


class cotton2k_build_ext(build_ext):
    def build_extensions(self):
        subprocess.call(["cargo", "build", "--release"])
        args = extra_compile_args[self.compiler.compiler_type]
        for extension in self.extensions:
            extension.sources = glob("src/_cotton2k/*.cpp")
            extension.libraries = libraries[os.name]
            extension.library_dirs = ["target/release"]
            extension.extra_compile_args = args
        super().build_extensions()


setup(
    packages=["cotton2k"],
    package_dir={"": "src"},
    package_data={"cotton2k": ["*.json", "*.csv"]},
    ext_modules=extensions,
    cmdclass={"build_ext": cotton2k_build_ext},
)
