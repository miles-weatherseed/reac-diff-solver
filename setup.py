from setuptools import setup, find_packages

long_description = ""
with open("README.md", "r") as fi:
    long_description = fi.read()

setup(
    name="reac-diff-solver",
    version="0.0.1",
    author="Barnum Swannell, Danail Stoychev, Miles Weatherseed, Muriel van der Laan",
    author_email="barnum.swannell@worc.ox.ac.uk, danail.stoychev@exeter.ox.ac.uk, miles.weatherseed@st-annes.ox.ac.uk, muriel.vanderlaan@wolfson.ox.ac.uk",
    description="A module for solving general reaction-diffusion systems",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/miles-weatherseed/scientific_computing_assessment",
    package_dir={'': "reac_diff_solver"},
    packages=find_packages(where="reac_diff_solver"),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=["numpy", "scipy"]
)
