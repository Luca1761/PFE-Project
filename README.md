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
This repository is organized as follows:
- The 

### How to launch experiments
