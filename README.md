# Reaction-Diffusion PDE System Solver
### Barnum Swannell, Danail Stoychev, Miles Weatherseed, Muriel van der Laan
### 29th November 2019

This module can solve any reaction-diffusion system of the form

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;u}{\partial&space;t}&space;=&space;D_u&space;\nabla^2u&space;&plus;&space;f(u,v)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\frac{\partial&space;u}{\partial&space;t}&space;=&space;D_u&space;\nabla^2u&space;&plus;&space;f(u,v)" title="\frac{\partial u}{\partial t} = D_u \nabla^2u + f(u,v)," /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;v}{\partial&space;t}&space;=&space;D_v&space;\nabla^2v&space;&plus;&space;g(u,v)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\frac{\partial&space;v}{\partial&space;t}&space;=&space;D_v&space;\nabla^2v&space;&plus;&space;g(u,v)" title="\frac{\partial v}{\partial t} = D_v \nabla^2v + g(u,v)" /></a>

The system is solved numerically using finite element and conjugate gradient methods. 

In addition to this, the user has the ability to take some observed data and recover the parameters from the original model using inference methods.

## Installation

First, clone the repository by opening the terminal in the desired directory and typing

```
git clone https://github.com/miles-weatherseed/reac-diff-solver.git
```

Then navigate to the newly created directory and install the module and all dependencies by typing

```
cd reac-diff-solver
pip install .
```

## API Reference

The API reference can be found at the following location:

https://miles-weatherseed.github.io/reac-diff-solver/

## PDE system solver

By default, the solver class expects a system of the format:


<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;u}{\partial&space;t}&space;=&space;D_u&space;\nabla^2u" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\frac{\partial&space;u}{\partial&space;t}&space;=&space;D_u&space;\nabla^2u" title="\frac{\partial u}{\partial t} = D_u \nabla^2u" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;v}{\partial&space;t}&space;=&space;D_v&space;\nabla^2v" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\frac{\partial&space;v}{\partial&space;t}&space;=&space;D_v&space;\nabla^2v" title="\frac{\partial v}{\partial t} = D_v \nabla^2v" /></a>

Initial conditions are defined via a function created by the user. For example, to provide the initial conditions for the Gray-Scott Equations, the user defines the function

```
def initial_conditions(X,Y):

    u = np.ones_like(X)
    v = np.zeros_like(X)

    # add perturbation in corner
    u[:20,:20] = 0.5
    v[:20,:20] = 0.25

    # add small amount of random noise to break symmetry
    u += 0.1*np.random.standard_normal(X.shape)
    v += 0.1*np.random.standard_normal(X.shape)

    return [u,v]
    
```

The user also has the option to add a nonlinear reaction function to the system. Again, this is done by defining a function. Continuing to look at the Gray-Scott Equations, the reaction functions are defined by

```
def Gray_Scott_reaction_terms(u,v, parameters = [1.0, 1.0]):
    k = parameters[0]
    F = parameters[1]
    return [-u*v**2 + F*(1-u), u*v**2 - (k+F)*v] [u,v]
    
```
This information, along with the values of parameters and the size and resolution of the finite differences grid, are then passed to the solver class. The output provides the values of u and v at each timestep at each point in the grid. This can then be visualized using in-built animation features, producing visually appealing outputs like the gifs below.

![](examples/GrayScottEquations_Spots.gif) ![](examples/GrayScottEquations_Stripes1.gif)

## Parameter inference
The user also has the option to tackle to problem in reverse, by providing some observed time series of u and v and recovering the values of the underlying system's parameters using inference. 

The system is set up in much the same way as the solving problem. First, we provide the inference class with the observed data

```
inference = Inference(u_data, v_data, times)
```

Next, we set the boundaries for the domain under consideration and provide a reaction function if we are expecting one


```
inference.set_model(xBounds, yBounds)
inference.set_reaction_function(function)
```

Once this is setup, we simply call the function *fit_params* which uses the Nelder-Mead algorithm to minimize an L2 error cost function. Note, this is not a trivial process and it will likely take a long time even on a fast machine. The user has to provide an initial estimate of the parameters in the form of *x0* which is a list *[D_u, D_v, ...]*.


```
print(inference.fit_params(x0))
```

## Unit Testing
To run the unit tests, first navigate to the src directory and then use the `unittest` module like this:

```
cd src
python -m unittest discover -s ../tst
```
