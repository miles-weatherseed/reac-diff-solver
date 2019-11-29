# Reaction-Diffusion PDE System Solver
### Barnum Swannell, Danail Stoychev, Miles Weatherseed, Muriel van der Laan
### 29th November 2019

This module can solve any reaction-diffusion system of the form

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;u}{\partial&space;t}&space;=&space;D_u&space;\nabla^2u&space;&plus;&space;f(u,v)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\frac{\partial&space;u}{\partial&space;t}&space;=&space;D_u&space;\nabla^2u&space;&plus;&space;f(u,v)" title="\frac{\partial u}{\partial t} = D_u \nabla^2u + f(u,v)," /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;v}{\partial&space;t}&space;=&space;D_v&space;\nabla^2u&space;&plus;&space;g(u,v)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\frac{\partial&space;v}{\partial&space;t}&space;=&space;D_v&space;\nabla^2u&space;&plus;&space;g(u,v)" title="\frac{\partial v}{\partial t} = D_v \nabla^2u + g(u,v)." /></a>

The system is solved numerically using finite element and conjugate gradient methods. 

In addition to this, the user has the ability to take some observed data and recover the parameters from the original model using inference methods.

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
