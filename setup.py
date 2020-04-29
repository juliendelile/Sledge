import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

install_requires = [
    'pandas>=1.0.3',
    'python-igraph>=0.8.1',
]

setuptools.setup(
    name="sledge-no-username",
    version="0.1.0",
    author="Julien Delile",
    author_email="julien.delile@gmail",
    description="Synthetic LinEage Dataset GEnerator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/juliendelile/Sledge",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=install_requires,
)