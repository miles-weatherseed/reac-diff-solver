.. toctree::
   :maxdepth: 2
   :caption: Table of Contents

Introduction
============

This Python module provides a library for solving reaction-diffusion systems, using finite element and conjugate gradient methods. The general form of the systems is:

.. math::

    \begin{align*}
    \frac{\partial u}{\partial t} &= \text{D}_u \nabla^{2} u + f(u,v) \\
    \frac{\partial v}{\partial t} &= \text{D}_v \nabla^{2} v + g(u,v)
    \end{align*}

Installation
------------

Clone the repo and install via pip:

.. code-block::

    git clone https://github.com/miles-weatherseed/reac-diff-solver.git
    cd reac_diff_solver
    pip install .

Example
-------

See README for a detailed example on how to use the library.

API Reference
=============

.. automodule:: reac_diff_solver.conjgrad
   :members:
   :undoc-members:
   :show-inheritance:


.. automodule:: reac_diff_solver.inference
   :members:
   :undoc-members:
   :show-inheritance:


.. automodule:: reac_diff_solver.solver
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: reac_diff_solver.animation
   :members:
   :undoc-members:
   :show-inheritance:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
