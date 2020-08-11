from setuptools import setup, find_packages


setup(
    name='cotton2k',
    version="2020.08.11",
    author="Tang Ziya",
    author_email="tcztzy@gmail.com",
    description="Reimplementation for Cotton2k simulation model",
    keywords="cotton simulation model",
    packages=find_packages(),
    install_requires=["appdirs>=1.4.4"]
)
