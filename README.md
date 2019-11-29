# Reaction-Diffusion PDE System Solver
### Barnum Swannell, Danail Stoychev, Miles Weatherseed, Muriel van der Laan
### 29th November 2019

This module can solve any reaction-diffusion system of the form

$$\frac{du}{dt} = 1$$

using finite element and conjugate gradient methods. 

Furthermore, the user has the ability to take some observed data and recover the parameters from the original model using inference methods.

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
