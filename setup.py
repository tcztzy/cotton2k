from setuptools import setup, find_packages


CLASSIFIERS="""\
Development Status :: 2 - Pre-Alpha
Intended Audience :: Science/Research
License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)
Topic :: Scientific/Engineering :: Bio-Informatics
Programming Language :: Python :: 3
Programming Language :: Python :: 3.8"""


setup(
    name="cotton2k",
    version="2020.08.14",
    author="Tang Ziya",
    author_email="tcztzy@gmail.com",
    description="Reimplementation for Cotton2k simulation model",
    license="AGPLv3+",
    keywords="cotton simulation model",
    packages=find_packages(),
    install_requires=["appdirs>=1.4.4"],
    test_requires=["pytest", "pytest-cov"],
    python_requires=">=3.8",
    classifiers=CLASSIFIERS.splitlines(),
)
