import toml
from skbuild import setup

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
print(project)
setup(
    packages=["cotton2k"],
    package_dir={"": "src"},
    package_data={"cotton2k": ["*.json", "*.csv"]},
    **project
)
