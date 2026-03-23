# Welcome to [Chibli, Genet & Imperiale. Stability analysis of a new curl-based full field reconstruction method in 2D isotropic nearly-incompressible elasticity. Submitted.]'s demos!

Static and interactive demos can be found at [https://mgenet.github.io/Curl-Based-Reconstruction-Paper-Demos](https://mgenet.github.io/Curl-Based-Reconstruction-Paper-Demos), or directly on [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mgenet/Curl-Based-Reconstruction-Paper-Demos/master?urlpath=lab/tree/./demos).

## Installation

A working installation of **FEniCS** (version 2019.1.0, including the Dolfin Python interface) and **VTK** (including the Python interface) is required to run the code.  
The simplest way to set up your system is using **Conda**:

```bash
conda create -y -c conda-forge -n fenics fenics=2019.1.0 matplotlib=3.5 mpi4py=3.1.3 numpy=1.24 scipy=1.10 pandas=1.3 pip python=3.10  
conda activate fenics
