# Code for paper on DSIRP

## Getting started

```
make
```

```
./irp path-to-instance -seed <seed> -veh <number-of-vehicle> -scenarios <nbScenarios> -iter <maxIter> -iterNonProd <maxIterNonProd> -time <maxtime> -traces <verbosity>
```

-Example: `./irp dsirp/standard/dirp-5-5-1.dat -seed 42 -veh 1 -scenarios 1 -iter 20000 -iterNonProd 20000 -time 3 -traces 1`

## Detailed Files Description
This repository is organized as follows:

### How to launch experiments