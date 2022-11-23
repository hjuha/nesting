# Nesting Decompositions

This repository contains the supplementary materials for the AAAI 2023 paper (to be published) "A Faster Practical Approximation Scheme for the Permanent" by Juha Harviainen and Mikko Koivisto.

## deepar.cpp

The file contains the implementations of rejection sampling schemes by Harviainen et al. (2021) ([Github](https://github.com/Kalakuh/deepar/tree/main/DeepAR)) and our implementation of ExtHuber. Compilation requires C++ library Boost, which comes automatically with many Linux distributions, as well as the file `Breal.hpp`. Before compilation, the file can be configured: 

`BOUND`: The used scheme – (AdaPart = 0, HL = 1, and ExtHuber = 2). HL requires that the matrix has entries in [0, 1].

`TASSA`: If true, entries that are not part of any permutation of positive weight are replaced by zeros.

`SHARPEN`: The level of preprocessing – 0: nothing, 1: matrix is made nearly doubly stochastic, 2: our sharpening scheme is applied

An (epsilon, delta)-approximation of an n × n matrix A is obtained by using a depth-d bound with the following input:

```
epsilon delta d time_limit
n
A_11  ..  ..  A_1n
 ..   ..       ..
 ..       ..   ..
A_n1  ..  ..  A_nn
```

The output contains four values: logarithm of the estimate, preprocessing time, sampling time, and total running time.

## sis.cpp

A slightly modified sequential importance sampler of Alimohammadi et al. (2021) ([Github](https://github.com/mohammadroghani/SIS)). Produces improving estimates of the permanent over time. The parameter `d` is ignored and and used only by rejection samplers.

Input format:

```
n d time_limit
A_11  ..  ..  A_1n
 ..   ..       ..
 ..       ..   ..
A_n1  ..  ..  A_nn
```

## is_plot.cpp

Our implementation of Smith & Dawkins' (2001) importance sampler PPS (R=1). Produces improving estimates of the permanent over time. The parameter `d` is ignored and and used only by rejection samplers.

Input format:
```
n d time_limit
A_11  ..  ..  A_1n
 ..   ..       ..
 ..       ..   ..
A_n1  ..  ..  A_nn
```

## generator.cpp, mixed_gen.py, generator.py, and gen.py

The scripts generate different kinds of matrices. They are not well-documented, but should be simple to modify and use. `generator.cpp` generates Bernoulli(p) instances with positive permanent, `generator.py` instance classes `Block Diagonal` and `Random Permutations`, `mixed_gen.py` the `Mixed` instance, and `gen.py` the rest.

## instances/

To enable replication of the results, this folder contains the instances used for evaluating the performance of the importance samplers. In addition, there is a script `visualizer.py` for visualizing the matrices.

# Example

Try, for example, the following script that generates a 10×10 matrix with Bernoulli(0.2)-distributed entries and gets an (0.1, 0.05)-approximation with ExtHuber-S using a depth-5 bound:
```
g++ generator.cpp -O2 -o bernoulli_generator
g++ deepar.cpp -O2 -o exthuber_sharp
echo "10 0.2 5" | ./bernoulli_generator
echo "10 0.2 5" | ./bernoulli_generator | ./exthuber_sharp
```
