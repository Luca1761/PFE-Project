# Code for paper on DSIRP

## Getting started
Firstly, the C++ code has to be compilated. To this end, use one of the following commands, depending your version of make: 
```
make or mingw32-make
```

You also can clear it using the `clean` command.

Then, the command to launch the code uses the following pattern:
```
./irp path-to-instance -seed <seed> -cores <number-of-cores> -endDay <0-or-1> -trueDemand1<0-or-1> -deterministic <0-or-1> -veh <number-of-vehicle> -scenarios <number-of-scenarios> -iter <max-iteration> -iterNonProd <max-iter-without-improve> -time <max-time> -traces <verbosity>
```
The different parameters are:
- `seed`: the seed for reproductibility (without specifying, the seed is random).
- `cores`: the number of cores available. One thread (per scenario) will be assigned to each core.
- `endDay`: boolean that indicates when we limit inventories. If true, $I^{t-1}_i+q_i^t-d_i^t\leq U_i$. Otherwise, $I^{t-1}_i+q_i^t\leq U_i$.
- `trueDemand1`: boolean that indicates if we know demand of day 1. If true, we know $d_i^1$. Otherwise, we don't and it depends of scenario $\omega$: $d_i^{1,\omega}$.
- `deterministic`: boolean that indicates if we should apply _HGS-IRP_ or _HGS-DSIRP_ (to compare the stochastic version to the deterministic one).
- `veh`: the numbers of available vehicles. In our experiments and to compare to coelho paper, only one vehicle is available.
- `scenarios`: the number of scenarios to generate to take our decision. In our experiments, we set this number to 10.
- `iter`: the maximum number of iterations allowed to each step of the horizon. In our experiments, we set this number to 20000.
- `iterNonProd`: the maximum number of iterations without improving the best known solution (again, for each step of the horizon). In our experiments, we set this number to 20000
- `time`: the maximum amount of time allow to each step of the horizon. In our experiments, we set this number to 1200.
- `traces`: the verbosity of the optimization. If 0, nothing appears. If 1, the full trace of resolution is shown.

Concrete example: `./irp instances/seasonal/dirp-5-5-1.dat -seed 42 -cores 4 -endDay 0 -trueDemand1 0 -deterministic 0 -veh 1 -scenarios 5 -iter 20000 -iterNonProd 20000 -time 10 -traces 1`

## Detailed Files Description
Here, you can find the description of every different file in this repository:
- The `instances` folder folder contains all the instances (standard, seasonal, correlated) we use for our experiences. They come from http://www.leandro-coelho.com/instances/stochastic-and-dynamic-inventory-routing/. A README file is also here to remind the structure of these instances.
- The `results` folder contains all the experiences results, Figures and solutions. A README file is also here to describe their structure. Currently, we only did all the experiences for standard dataset. Experiences for seasonal and correlated should be done by adapting scenarios generation.
- `expes.py` is used to in `job.sh` to launch experiences on compute canada with given expe id. Parameters for experiences can be modified in this file.
- `commandline.h` and `commandline.cpp`
- `Client.h` and `Client.cpp` contain the definition of the Client class, related to customers (and also supplier).
- `Genetic.h` and `Genetic.cpp` contain the definition of the Genetic class, used to evolve our population.
- `Individual.h` and `Individual.cpp` contain the definition of the Individual, used to represent a solution with the tour and delivered quantities.
- `LinearPiece.h` and `LinearPiece.cpp` contain the definition of linear pieces, used to build piecewise linear functions.
- `PLFunction.h` and `PLFunction.cpp` contain the definition of piecewise linear functions to compute backward dynamic programming.
- `LotSizingSolver.h` and `LotSizingSolver.cpp` containg the definition of our Dynamic Programming Operator.
- `main.cpp` launch the full _HGS-DSIRP_ resolution for a given instance and given parameters.
- `Mutations.cpp` contain every mutation from _HGS_.
- `Node.h` and `Node.cpp`
- `Params.h` and `Params.cpp` contain the definition of the Params class, containing every instance and solve parameters.
- `Population.h` and `Population.cpp` contain the definition of the Population class, containing every solution we use for the genetic algorithm.
- `Rng.h` and `Rng.cpp` contain the definition of a random generator class.
- `Route.h` and `Route.cpp` contain the definition of the Route class used to define our tours.
- `Vehicle.h` and `Vehicle.cpp` contain the definition of the Vehicle class also use to define our tours.


### How to launch experiments
To Launch experiments, one must modify the parameters (and the type of dataset: standard, seasonal or correlated) from the `expes.py` file and use `job.sh` to launch them on a computational structure as compute canada. All the experiences on standard dataset have been done (deterministic- HGS-DSIRP and HGS-DSIRPR (endDay = true)) but correlated and seasonal need some changes in the scenario generation (because they have some patterns that need to be take into account for generation). Particularly, we use in our experiments `seed=42`, `cores=10`, `trueDemand1=0`, `-veh=1`, `scenarios=10`, `iter=20000`, `iterNonProd=20000`, `time=1200`, `traces=0`.  
