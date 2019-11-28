from setuptools import setup

setup(
   name='reac-diff-solver',
   version='1.0',
   description='A module for solving general reaction-diffusion systems',
   author='Miles Weatherseed, Muriel van der Laan, Barnum Swannell, Danail Stoychev',
   packages=['reac-diff-solver'],  #same as name
   install_requires=['numpy', 'scipy'], #external packages as dependencies
)