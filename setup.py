import subprocess
from pathlib import Path

import toml
from Cython.Build import cythonize
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

with open("pyproject.toml") as pyproject_toml:
    project = toml.load(pyproject_toml)["project"]

with open(project.pop("readme"), encoding="utf-8") as readme:
    project["long_description"] = readme.read()
    project["long_description_content_type"] = "text/x-rst"

author = project.pop("authors")[0]
project["author"] = author["name"]
project["author_email"] = author["email"]
project["python_requires"] = project.pop("requires-python")
project["install_requires"] = project.pop("dependencies", None)
urls = project.pop("urls")
project["url"] = urls.pop("homepage")
project["project_urls"] = urls

extensions = [
    Extension(
        "_cotton2k",
        sources=[
            "src/_cotton2k/__init__.pyx",
            *map(
                str,
                filter(
                    lambda p: p.stem != "__init__", Path("src/_cotton2k").glob("*.cpp")
                ),
            ),
        ],
        library_dirs=["target/release"],
        libraries=["cotton2k"],
        extra_compile_args=["-std=c++20"],
    )
]


class cotton2k_build_ext(build_ext):
    def build_extension(self, ext):
        subprocess.call(["cargo", "build", "--release"])
        super().build_extension(ext)


setup(
    packages=["cotton2k"],
    package_dir={"": "src"},
    package_data={"cotton2k": ["*.json", "*.csv"]},
    ext_modules=cythonize(extensions),
    cmdclass={"build_ext": cotton2k_build_ext},
    **project
)
