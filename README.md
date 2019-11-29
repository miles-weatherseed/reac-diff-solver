# Reaction-Diffusion System Solver
This repository can solve any reaction-diffusion system in 2D using finite element and conjugate gradient methods. 

There is also the ability to take a time series and use inference to recover the diffusion parameters from the underlying system.

## Unit Testing
To run the unit tests, first navigate to the src directory and then use the `unittest` module like this:

```
cd src
python -m unittest discover -s ../tst
```
