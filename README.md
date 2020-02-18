# Super-resolution of near-colliding point sources
This repository holds the code reproducing the numerical results from the paper D. Batenkov, G. Goldman, and Y. Yomdin, “Super-resolution of near-colliding point sources,” arXiv:1904.09186 [math], Apr. 2019, to appear in *Information and Inference: a Journal of the IMA*.



## Prerequisites

- [Julia](http://julialang.org) >= 1.2
- Julia packages: **`Polynomials`, `Combinatorics`,` LinearAlgebra`, `Random`, `Statistics`**
- Plotting packages: **`Plots`, `PyPlot`**
- In order to run the notebooks, you should also install **`IJulia`**.

## Instructions

In order to reproduce the figures, please execute the Jupyter notebooks:

- [Figs3and4.ipynb](./Figs3and4.ipynb) for Figures 3 and 4;
- [Figs5and6.ipynb](./Figs5and6.ipynb) for Figures 5 and 6.

The above notebooks call the code from the corresponding Julia source files: [fig3.jl](./fig3.jl), [fig4.jl](./fig4.jl), [fig51.jl](./fig51.jl), [fig52.jl](./fig52.jl), and [fig6.jl](./fig6.jl). These can be executed directly from the Julia prompt, in which case you do not need **`IJulia`**:

`include("./fig3.jl")`

and so on. The figure files are written to the  `figures` subfolder.

