import logging
import os
import subprocess
from collections import defaultdict
from glob import glob

import setuptools
from pyproject_toml import setup
from setuptools.command.develop import develop

log = logging.getLogger("COTTON2K")
try:
    from multiprocessing import cpu_count

    from Cython.Build import build_ext, cythonize

    extensions = cythonize(
        "src/_cotton2k/*.pyx", nthreads=cpu_count() if os.name != "nt" else 0
    )
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

extra_compile_args = defaultdict(lambda: ["-std=c++20"])
extra_compile_args["msvc"] = ["/std:c++latest"]
libraries = defaultdict(lambda: ["cotton2k"])
libraries["nt"].extend(["ws2_32", "userenv", "advapi32"])


class cotton2k_build_ext(build_ext):
    def build_extensions(self):
        cargo_build = ["cargo", "build"]
        if self.debug:
            self.cython_directives = {"linetrace": True}
        else:
            cargo_build.append("--release")
        subprocess.call(cargo_build)
        args = extra_compile_args[self.compiler.compiler_type]
        for extension in self.extensions:
            extension.sources = glob("src/_cotton2k/*.cpp")
            extension.libraries = libraries[os.name]
            extension.library_dirs = [
                "target/" + ("debug" if self.debug else "release")
            ]
            extension.extra_compile_args = args
            if self.debug:
                extension.define_macros = [("CYTHON_TRACE", 1)]
        super().build_extensions()


class cotton2k_develop(develop):
    user_options = develop.user_options + [
        ("debug", "g", "compile/link with debugging information")
    ]
    boolean_options = develop.boolean_options + ["debug"]

    def initialize_options(self):
        super().initialize_options()
        self.debug = 0

    def install_for_development(self):
        # Without 2to3 inplace works fine:
        self.run_command("egg_info")

        # Build extensions in-place
        self.reinitialize_command("build_ext", inplace=1, debug=self.debug)
        self.run_command("build_ext")

        if setuptools.bootstrap_install_from:
            self.easy_install(setuptools.bootstrap_install_from)
            setuptools.bootstrap_install_from = None

        self.install_namespaces()

        # create an .egg-link in the installation dir, pointing to our egg
        log.info("Creating %s (link to %s)", self.egg_link, self.egg_base)
        if not self.dry_run:
            with open(self.egg_link, "w") as f:
                f.write(self.egg_path + "\n" + self.setup_path)
        # postprocess the installed distro, fixing up .pth, installing scripts,
        # and handling requirements
        self.process_distribution(None, self.dist, not self.no_deps)


setup(
    packages=["cotton2k", "_cotton2k"],
    package_dir={"": "src"},
    package_data={"cotton2k": ["*.json", "*.csv"]},
    ext_modules=extensions,
    cmdclass={"build_ext": cotton2k_build_ext, "develop": cotton2k_develop},
)
