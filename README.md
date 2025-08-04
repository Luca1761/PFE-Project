# Code for paper on Inventory Routing Problem

## Setup CPlex
- Updating...

## Build
```
make
```
## How to run
```
./irp path-to-instance -seed <seed> -veh <number-of-vehicle> -scenarios <nbScenarios> -iter <maxIter> -iterNonProd <maxIterNonProd> -time <maxtime> -traces <verbosity>
```

-Example: `./irp dsirp/standard/dirp-5-5-1.dat -seed 42 -veh 1 -scenarios 1 -iter 20000 -iterNonProd 20000 -time 3 -traces 1`

## stockout test4
./irp Data/Small/Istanze0105h3/abs1n10_1.dat -seed 1000 -veh 2 -stock 100