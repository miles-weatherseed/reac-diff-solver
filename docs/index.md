<div class="related" role="navigation" aria-label="related navigation">

### Navigation

*   [index](genindex.html "General Index")
*   [modules](py-modindex.html "Python Module Index")
*   [reac-diff-solver 0.0.1 documentation](#) »

</div>

<div class="document">

<div class="documentwrapper">

<div class="bodywrapper">

<div class="body" role="main">

<div class="section" id="introduction">

# Introduction[¶](#introduction "Permalink to this headline")

This Python module provides a library for solving reaction-diffusion systems, using finite element and conjugate gradient methods. The general form of the systems is:

<div class="math notranslate nohighlight">\[\begin{split}\begin{align*} \frac{\partial u}{\partial t} &= \text{D}_u \nabla^{2} u + f(u,v) \\ \frac{\partial v}{\partial t} &= \text{D}_v \nabla^{2} v + g(u,v) \end{align*}\end{split}\]</div>

<div class="section" id="installation">

## Installation[¶](#installation "Permalink to this headline")

Clone the repo and install via pip:

<div class="highlight-default notranslate">

<div class="highlight">

<pre><span></span><span class="n">git</span> <span class="n">clone</span> <span class="n">https</span><span class="p">:</span><span class="o">//</span><span class="n">github</span><span class="o">.</span><span class="n">com</span><span class="o">/</span><span class="n">miles</span><span class="o">-</span><span class="n">weatherseed</span><span class="o">/</span><span class="n">reac</span><span class="o">-</span><span class="n">diff</span><span class="o">-</span><span class="n">solver</span><span class="o">.</span><span class="n">git</span>
<span class="n">cd</span> <span class="n">reac_diff_solver</span>
<span class="n">pip</span> <span class="n">install</span> <span class="o">.</span>
</pre>

</div>

</div>

</div>

<div class="section" id="example">

## Example[¶](#example "Permalink to this headline")

See README for a detailed example on how to use the library.

</div>

</div>

<div class="section" id="module-reac_diff_solver.conjgrad"><span id="api-reference"></span>

# API Reference[¶](#module-reac_diff_solver.conjgrad "Permalink to this headline")

<dl class="function">

<dt id="reac_diff_solver.conjgrad.conjugate_gradients">`reac_diff_solver.conjgrad.``conjugate_gradients`<span class="sig-paren">(</span>_A_, _b_, _x0_, _tol=0.0001_, _nmax=100_<span class="sig-paren">)</span>[¶](#reac_diff_solver.conjgrad.conjugate_gradients "Permalink to this definition")</dt>

<dd>

Uses the conjugate gradient algorithm to solve the system Ax = b.

<dl class="field-list simple">

<dt class="field-odd">Parameters</dt>

<dd class="field-odd">

*   **A** (_2d numpy array_ _or_ _scipy sparse matrix_) – symmetric, positive definite matrix

*   **b** (_1d numpy array_) – right-hand side of the system

*   **x0** (_1d numpy array_) – starting point of the iteration

*   **tol** (_float_) – tolerance, iteration terminates when the norm of the residual has been reduced by the tolerance

*   **nmax** (_int_) – maximum number of iterations

</dd>

<dt class="field-even">Raises</dt>

<dd class="field-even">

*   **TypeError** – the matrix A is not square

*   **TypeError** – the dimensions of A & b don’t match

*   **TypeError** – the dimensions of A & x0 don’t match

*   **Exception** – the iteration has failed to converge in nmax iterations

</dd>

<dt class="field-odd">Returns</dt>

<dd class="field-odd">

solution of the equation, list of the norm of the residual after each iteration

</dd>

<dt class="field-even">Return type</dt>

<dd class="field-even">

tuple of 1d numpy array and a list of floats

</dd>

</dl>

</dd>

</dl>

<span class="target" id="module-reac_diff_solver.inference"></span>

<dl class="class">

<dt id="reac_diff_solver.inference.Inference">_class_ `reac_diff_solver.inference.``Inference`<span class="sig-paren">(</span>_u_data_, _v_data_, _times_<span class="sig-paren">)</span>[¶](#reac_diff_solver.inference.Inference "Permalink to this definition")</dt>

<dd>

Bases: [`<span class="pre">reac_diff_solver.solver.Solver</span>`](#reac_diff_solver.solver.Solver "reac_diff_solver.solver.Solver")

<dl class="method">

<dt id="reac_diff_solver.inference.Inference.cost_func">`cost_func`<span class="sig-paren">(</span>_parametervalues_<span class="sig-paren">)</span>[¶](#reac_diff_solver.inference.Inference.cost_func "Permalink to this definition")</dt>

<dd>

Takes newly proposed set of parameter values and returns the L2 error between the values of u,v given by the solver with these values and the values of u,v from the input

<dl class="field-list simple">

<dt class="field-odd">Parameters</dt>

<dd class="field-odd">

**parametervalues** – the proposed set of parameter values at this point in the optimization

</dd>

<dt class="field-even">Returns</dt>

<dd class="field-even">

the L2 error with these parameters

</dd>

</dl>

</dd>

</dl>

<dl class="method">

<dt id="reac_diff_solver.inference.Inference.fit_params">`fit_params`<span class="sig-paren">(</span>_x0_<span class="sig-paren">)</span>[¶](#reac_diff_solver.inference.Inference.fit_params "Permalink to this definition")</dt>

<dd>

The master function to call. Fits parameters from some initial estimate by minimizing the cost function using the Nelder-Mead method.

<dl class="field-list simple">

<dt class="field-odd">Parameters</dt>

<dd class="field-odd">

**x0** – the initial estimate of parameters

</dd>

<dt class="field-even">Return best_fit_parameters</dt>

<dd class="field-even">

a vector of the fitted parameters

</dd>

</dl>

</dd>

</dl>

<dl class="method">

<dt id="reac_diff_solver.inference.Inference.set_model">`set_model`<span class="sig-paren">(</span>_xBounds_, _yBounds_, _initial_conditions_function_<span class="sig-paren">)</span>[¶](#reac_diff_solver.inference.Inference.set_model "Permalink to this definition")</dt>

<dd>

Sets up the reaction diffusion model we wish to consider. We determine the x and y boundaries for consideration and provide some initial conditions for our system.

<dl class="field-list simple">

<dt class="field-odd">Parameters</dt>

<dd class="field-odd">

*   **xBounds** – x-range of the problem

*   **yBounds** – y-range of the problem

*   **initial_condition_function** (_function that takes two 2d numpy arrays and returns a list of two 2d numpy arrays_) – calculates the values of u and v at t=0 at each gridpoint

</dd>

</dl>

</dd>

</dl>

<dl class="method">

<dt id="reac_diff_solver.inference.Inference.set_reaction_function">`set_reaction_function`<span class="sig-paren">(</span>_function_<span class="sig-paren">)</span>[¶](#reac_diff_solver.inference.Inference.set_reaction_function "Permalink to this definition")</dt>

<dd>

This is optional. We can choose to provide a nonlinear term to our system.

<dl class="field-list simple">

<dt class="field-odd">Parameters</dt>

<dd class="field-odd">

**function** (_function that takes two numpy arrays_ _(__containing values of u and v__)_ _and a list of parameters and returns a list of two numpy arrays_) – calculates the value of the reaction terms at the given u and v

</dd>

</dl>

</dd>

</dl>

</dd>

</dl>

<span class="target" id="module-reac_diff_solver.solver"></span>

<dl class="class">

<dt id="reac_diff_solver.solver.Solver">_class_ `reac_diff_solver.solver.``Solver`<span class="sig-paren">(</span>_xBounds_, _yBounds_, _gridSize_, _initial_condition_function_<span class="sig-paren">)</span>[¶](#reac_diff_solver.solver.Solver "Permalink to this definition")</dt>

<dd>

Bases: `<span class="pre">object</span>`

<dl class="method">

<dt id="reac_diff_solver.solver.Solver.set_grid">`set_grid`<span class="sig-paren">(</span>_xBounds_, _yBounds_, _gridSize_<span class="sig-paren">)</span>[¶](#reac_diff_solver.solver.Solver.set_grid "Permalink to this definition")</dt>

<dd>

Sets the grid parameters.

<dl class="field-list simple">

<dt class="field-odd">Parameters</dt>

<dd class="field-odd">

*   **xBounds** (_list of two floats_) – x-range of the problem

*   **yBounds** (_list of two floats_) – y-range of the problem

*   **gridSize** (_int_) – number of gridpoints to use (in both the x- & y-direction)

</dd>

<dt class="field-even">Returns</dt>

</dl>

</dd>

</dl>

<dl class="method">

<dt id="reac_diff_solver.solver.Solver.set_initialConditions">`set_initialConditions`<span class="sig-paren">(</span>_initial_condition_function_<span class="sig-paren">)</span>[¶](#reac_diff_solver.solver.Solver.set_initialConditions "Permalink to this definition")</dt>

<dd>

Set the initial conditions used by the solver.

<dl class="field-list simple">

<dt class="field-odd">Parameters</dt>

<dd class="field-odd">

**initial_condition_function** (_function that takes two 2d numpy arrays and returns a list of two 2d numpy arrays_) – calculates the values of u and v at t=0 at each gridpoint

</dd>

<dt class="field-even">Returns</dt>

</dl>

</dd>

</dl>

<dl class="method">

<dt id="reac_diff_solver.solver.Solver.set_reactionFunction">`set_reactionFunction`<span class="sig-paren">(</span>_function_<span class="sig-paren">)</span>[¶](#reac_diff_solver.solver.Solver.set_reactionFunction "Permalink to this definition")</dt>

<dd>

Set the reaction term of the equation (defaults to zero).

<dl class="field-list simple">

<dt class="field-odd">Parameters</dt>

<dd class="field-odd">

**function** (_function that takes two numpy arrays_ _(__containing values of u and v__)_ _and a list of parameters and returns a list of two numpy arrays_) – calculates the value of the reaction terms at the given u and v

</dd>

<dt class="field-even">Returns</dt>

</dl>

</dd>

</dl>

<dl class="method">

<dt id="reac_diff_solver.solver.Solver.set_timeStepLength">`set_timeStepLength`<span class="sig-paren">(</span>_length_<span class="sig-paren">)</span>[¶](#reac_diff_solver.solver.Solver.set_timeStepLength "Permalink to this definition")</dt>

<dd>

Set the size of the timestep used by the solver.

<dl class="field-list simple">

<dt class="field-odd">Parameters</dt>

<dd class="field-odd">

**length** (_float_) – timestep

</dd>

<dt class="field-even">Returns</dt>

</dl>

</dd>

</dl>

<dl class="method">

<dt id="reac_diff_solver.solver.Solver.solve">`solve`<span class="sig-paren">(</span>_times_, _parameters_, _printProgress=False_<span class="sig-paren">)</span>[¶](#reac_diff_solver.solver.Solver.solve "Permalink to this definition")</dt>

<dd>

Solves the equation at the given times.

<dl class="field-list simple">

<dt class="field-odd">Parameters</dt>

<dd class="field-odd">

*   **times** (_list of floats_) – times at which the solution is desired.

*   **parameters** (_list_) – parameters to give to the reaction function

*   **printProgress** – whether or not to print progress updates

</dd>

<dt class="field-even">Returns</dt>

<dd class="field-even">

the solution of u and v at the given times

</dd>

<dt class="field-odd">Return type</dt>

<dd class="field-odd">

list of two 2d numpy arrays (self.gridSize by self.gridSize)

</dd>

</dl>

</dd>

</dl>

</dd>

</dl>

<span class="target" id="module-reac_diff_solver.animation"></span>

<dl class="function">

<dt id="reac_diff_solver.animation.animate">`reac_diff_solver.animation.``animate`<span class="sig-paren">(</span>_u_, _total_time_, _filename=''_<span class="sig-paren">)</span>[¶](#reac_diff_solver.animation.animate "Permalink to this definition")</dt>

<dd>

Creates an animation of u

<dl class="field-list simple">

<dt class="field-odd">Parameters</dt>

<dd class="field-odd">

*   **u** (_list of 2d numpy arrays_) – frames to show

*   **total_time** (_float_) – total time over which the animation should run

*   **filename** (_string_) – file to which to save the animation, if no filename is provided, the animation will not be saved

</dd>

<dt class="field-even">Returns</dt>

</dl>

</dd>

</dl>

</div>

<div class="section" id="indices-and-tables">

# Indices and tables[¶](#indices-and-tables "Permalink to this headline")

*   [<span class="std std-ref">Index</span>](genindex.html)

*   [<span class="std std-ref">Module Index</span>](py-modindex.html)

*   [<span class="std std-ref">Search Page</span>](search.html)

</div>

</div>

</div>

</div>

<div class="sphinxsidebar" role="navigation" aria-label="main navigation">

<div class="sphinxsidebarwrapper">

### [Table of Contents](#)

*   [Introduction](#)
    *   [Installation](#installation)
    *   [Example](#example)
*   [API Reference](#module-reac_diff_solver.conjgrad)
*   [Indices and tables](#indices-and-tables)

<div role="note" aria-label="source link">

### This Page

*   [Show Source](_sources/index.rst.txt)

</div>

<div id="searchbox" style="display: none" role="search">

### Quick search

<div class="searchformwrapper">

<form class="search" action="search.html" method="get"><input type="text" name="q" aria-labelledby="searchlabel"> <input type="submit" value="Go"></form>

</div>

</div>

<script type="text/javascript">$('#searchbox').show(0);</script></div>

</div>

</div>

<div class="related" role="navigation" aria-label="related navigation">

### Navigation

*   [index](genindex.html "General Index")
*   [modules](py-modindex.html "Python Module Index") |
*   [reac-diff-solver 0.0.1 documentation](#) »

</div>

<div class="footer" role="contentinfo">© Copyright 2019, Barnum Swannell, Danail Stoychev, Miles Weatherseed, Muriel van der Laan. Created using [Sphinx](http://sphinx-doc.org/) 2.2.0.</div>
