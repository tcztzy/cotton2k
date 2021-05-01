import os
import subprocess
from collections import defaultdict
from glob import glob

from Cython.Build import cythonize
from pyproject_toml import setup
from setuptools import Extension
from setuptools.command.build_ext import build_ext

extra_compile_args = defaultdict(lambda: ["-std=c++20"])
extra_compile_args["msvc"] = ["/std:c++latest"]
libraries = defaultdict(lambda: ["cotton2k"])
extra_compile_args["nt"].extend(["ws2_32", "userenv", "advapi32"])


class cotton2k_build_ext(build_ext):
    def build_extension(self, ext):
        subprocess.call(["cargo", "build", "--release"])
        args = extra_compile_args[self.compiler.compiler_type]
        for extension in self.extensions:
            extension.libraries = libraries[os.name]
            extension.extra_compile_args = args
        super().build_extension(ext)


setup(
    packages=["cotton2k"],
    package_dir={"": "src"},
    package_data={"cotton2k": ["*.json", "*.csv"]},
    ext_modules=cythonize(
        [
            Extension(
                "_cotton2k",
                sources=[
                    "src/_cotton2k/__init__.pyx",
                    *filter(
                        lambda p: "__init__" not in p,
                        glob("src/_cotton2k/*.cpp"),
                    ),
                ],
                library_dirs=["target/release"],
            )
        ]
    ),
    cmdclass={"build_ext": cotton2k_build_ext},
)
