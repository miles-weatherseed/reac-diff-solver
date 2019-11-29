# Reaction-Diffusion PDE System Solver
### Barnum Swannell, Danail Stoychev, Miles Weatherseed, Muriel van der Laan
### 29th November 2019

This repository can solve any reaction-diffusion system in 2D using finite element and conjugate gradient methods. 

There is also the ability to take a time series and use inference to recover the diffusion parameters from the underlying system.

## Installation

First, clone the repository by opening the terminal in the desired directory and typing

```
git clone https://github.com/miles-weatherseed/reac-diff-solver.git
```

Then install the module and all dependencies by typing


```
pip install .
```

## PDE system solver

## Parameter inference

## Unit Testing
To run the unit tests, first navigate to the src directory and then use the `unittest` module like this:

```
cd src
python -m unittest discover -s ../tst
```
