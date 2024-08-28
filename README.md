[![Paper](https://img.shields.io/badge/paper-arXiv%3A2404.17676-B31B1B.svg)](https://arxiv.org/abs/2404.17676)

# Toward a 2D Local Implementation of Quantum LDPC Codes

[Noah Berthusen](https://noahberthusen.github.io), Dhruv Devulapalli, Eddie Schoute, Andrew M. Childs, Michael J. Gullans, Alexey V. Gorshkov, Daniel Gottesman

### Abstract
Geometric locality is an important theoretical and practical factor for quantum low-density parity-check (qLDPC) codes which affects code performance and ease of physical realization. For device architectures restricted to 2D local gates, naively implementing the high-rate codes suitable for low-overhead fault-tolerant quantum computing incurs prohibitive overhead. In this work, we present an error correction protocol built on a bilayer architecture that aims to reduce operational overheads when restricted to 2D local gates by measuring some generators less frequently than others. We investigate the family of bivariate bicycle qLDPC codes and show that they are well suited for a parallel syndrome measurement scheme using fast routing with local operations and classical communication (LOCC). Through circuit-level simulations, we find that in some parameter regimes bivariate bicycle codes implemented with this protocol have logical error rates comparable to the surface code while using fewer physical qubits.

### Description
This repository includes information, code, and data to generate the figures in the paper.

### Figures
All the codes used to create the figures in the paper are found in the `/figures/figure_scripts` folder. They are all written in Python, and use the matplotlib library. Files used to draw certain figures in the paper can be found in the `/figures/figure_svgs` folder.
- `bell_fidelity.py` Generates Fig. 6(a).
- `no_idle_sim_results.py` Generates Fig. 12.
- `potential_savings.py` Generates Fig. 8.
- `routing_algo.py` Generates Fig. 3.
- `sim_results.py` Generates Fig. 11.

### Simulations

#### Circuit-level simulations
Code used to perform the circuit-level simulations of the bivariate bicycle codes can be found in the `/src_py/bivariate_bicycle` directory. Listed below are the files and their functions:
- `code_distance.g` [GAP](https://www.gap-system.org/) program which calculates the distance of a CSS code using the [QDistRnd](https://github.com/QEC-pages/QDistRnd) package. The files referenced, `QX.mtx` and `QZ.mtx` are the X and Z parity check matrices of the code in the [Matrix Market](https://networkrepository.com/mtx-matrix-market-format.html) file format.
- `code_exploration.py` Performs a computer search for BB codes which satisfy a number of user defined parameters. The user can enter $\ell$ and $m$, and the program will find BB codes with a toric embedding encoding a nonzero number of logical qubits. For the purposes of the paper, valid codes then went through preliminary simulations to determine if they were good candidates for the scheme as a whole.
- `general.py` Main driver file for the circuit-level simulations as defined in Section IV of the paper. To simulate a BB code, the user inputs its $A$ and $B$ polynomials, embedding parameters, and (short- and long-range) routing depth. The latter of which can be determine with the code in `/routing`. Circuit-level simulations powered by [Stim](https://github.com/quantumlib/Stim) are then performed for a user-defined number of rounds, with the long-range generators being measured every `lr_time` rounds.
- `mec.py` Calculates the minimal enclosing circle of a set of $(x,y)$ points.
- `result.py` Handles the format and saving of the simulation results.
- `teleportation.py` Performs circuit-level simulations of the Bell pair purification and application of the long-range CNOT gate. The user can specify the length of the long-range CNOT to perform, and then the purification process as described in Section IVB is performed. This file was used to generate the data presented in Fig. 6(a).

#### Routing simulations

Python code used to perform the simulations of Fig. 3 are located in `/src_py/routing`, while C++ code with the same functionality is located in `/src_cpp/routing`.
- `Grid.py` Defines a `Grid` class for routing a set of generators under the assumptions and conditions described in Section IIC. Given a $M\times N$ grid of "qubits", and a set of generators given as a list of $(x,y)$ locations of the qubits which they contain, the program attempt to efficiently route each of the generators. The algorithm used is described in Algorithm 1. The output is the number of steps required to route all of the generators.
- `mec.py` Calculates the minimal enclosing circle of a set of $(x,y)$ points.
- `route.py` Driver code for gathering the data used to create Fig. 3. Random BB code instances were generated, and then their generators were routed using `Grid.py` and the greedy routing algorithm. Given a set of generators defined as a list of lists of $(x,y)$ coordinates, the algorithm can be executed manually.

The C++ code located in `/src_cpp/routing` is a more efficient implementation of the Python code and can be used for large-scale investigations into the performance of the greedy routing algorithm. The functionality is identical and is driven from `main.cpp`. The user can define a `vector<vector<Point>>` to describe the $(x,y)$ locations of qubits in generators. It can then be passed into the `grid.cpp` class which also implements the routing algorithm. The output is again the number of steps required to route all of the provided generators.
