# Code for paper on DSIRP

## Getting started
Firstly, the C++ code has to be compilated. To this end, use one of the following commands, depending your version of make: 
```
make or mingw32-make
```
Then, the command to launch the code uses the following pattern:
```
./irp path-to-instance -seed <seed> -veh <number-of-vehicle> -scenarios <number-of-scenarios> -iter <max-iteration> -iterNonProd <max-iter-without-improve> -time <max-time> -traces <verbosity>
```
The different parameters are:
- `seed`: the seed for reproductibility (without specifying, the seed is random).
- `veh`: the numbers of available vehicles. In our experiments and to compare to coelho paper, only one vehicle is available.
- `scenarios`: the number of scenarios to generate to take our decision.
- `iter`: the maximum number of iterations allowed to each step of the horizon.
- `iterNonProd`: the maximum number of iterations without improving the best known solution (again, for each step of the horizon).
- `time`: the maximum amount of time allow to each step of the horizon.
- `traces`: the verbosity of the optimization. If 0, nothing appears. If 1, the full trace of resolution is shown.

-Concrete example: `./irp dsirp/standard/dirp-5-5-1.dat -seed 42 -veh 1 -scenarios 5 -iter 20000 -iterNonProd 20000 -time 10 -traces 1`

## Detailed Files Description
Here, you can find the description of every different file in this repository:
- The `data_dsirp folder` folder contains all the instances (standard, seasonal, correlated) we use for our experiences. They come from http://www.leandro-coelho.com/instances/stochastic-and-dynamic-inventory-routing/. A README file is also here to remind the structure of these instances.
- The `dsirp_results` folder contains all the experiences results. A README file is also here to describe their structure.
- `expes.py`
- `commandline.h` and `commandline.cpp`
- `Client.h` and `Client.cpp` contain the definition of the Client class, related to customers (and also supplier).
- `Genetic.h` and `Genetic.cpp`
- `Individu.h` and `Individu.cpp`
- `LinearPiece.h` and `LinearPiece.cpp`
- `LotSizingSolver.h` and `LotSizingSolver.cpp`
- `main.cpp`
- `Mutations.cpp`
- `Noeud.h` and `Noeud.cpp`
- `Params.h` and `Params.cpp`
- `PLFunction.h` and `PLFunction.cpp`
- `Population.h` and `Population.cpp`
- `Rng.h` and `Rng.cpp`
- `Route.h` and `Route.cpp`
- `Vehicle.h` and `Vehicle.cpp`


### How to launch experiments
