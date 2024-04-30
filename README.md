[![Paper](https://img.shields.io/badge/paper-arXiv%3A2404.17676-B31B1B.svg)](https://arxiv.org/abs/2404.17676)

# Toward a 2D Local Implementation of Quantum LDPC Codes

[Noah F. Berthusen](https://noahberthusen.github.io), Dhruv Devulapalli, Eddie Schoute, Andrew M. Childs, Michael J. Gullans, Alexey V. Gorshkov, Daniel Gottesman

### Abstract
Geometric locality is an important theoretical and practical factor for quantum low-density parity-check (qLDPC) codes which affects code performance and ease of physical realization. For device architectures restricted to 2D local gates, naively implementing the high-rate codes suitable for low-overhead fault-tolerant quantum computing incurs prohibitive overhead. In this work, we present an error correction protocol built on a bilayer architecture that aims to reduce operational overheads when restricted to 2D local gates by measuring some generators less frequently than others. We investigate the family of bivariate bicycle qLDPC codes and show that they are well suited for a parallel syndrome measurement scheme using fast routing with local operations and classical communication (LOCC). Through circuit-level simulations, we find that in some parameter regimes bivariate bicycle codes implemented with this protocol have logical error rates comparable to the surface code while using fewer physical qubits.

### Description
This repository includes information, code, and data to generate the figures in the paper.

### Figures
All the codes used to create the figures in the paper are found in the **figure_scripts** folder. They are all written in Python, and use the matplotlib library.
