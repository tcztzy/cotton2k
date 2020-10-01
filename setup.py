"""Setup file for cotton2k."""
from setuptools import find_packages, setup

CLASSIFIERS = """\
Development Status :: 2 - Pre-Alpha
Intended Audience :: Science/Research
License :: OSI Approved :: GNU Affero General Public License v3 or later \
(AGPLv3+)
Topic :: Scientific/Engineering :: Bio-Informatics
Programming Language :: Python :: 3
Programming Language :: Python :: 3.8"""


with open("README.md") as f:
    README = f.read()


setup(
    name="cotton2k",
    version="2020.10.1",
    author="Tang Ziya",
    author_email="tcztzy@hotmail.com",
    description="Reimplementation for Cotton2k simulation model",
    long_description=README,
    url="https://github.com/tcztzy/cotton2k",
    license="AGPLv3+",
    keywords="cotton simulation model",
    packages=find_packages(),
    install_requires=["PySimpleGUI>=4.29.0"],
    python_requires=">=3.8",
    classifiers=CLASSIFIERS.splitlines(),
)
