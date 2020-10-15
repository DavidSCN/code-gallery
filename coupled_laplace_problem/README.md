Black-box coupled Laplace equation
------------------------------------------
## Overview

What is preCICE? What is black-box coupling? What is this tutorial about?

## Requirements
* Version `9.2` or greater of `deal.II`

* [preCICE](https://github.com/precice/precice/wiki#1-get-precice)

## Compiling and running
Similar to the example programs, run
```
cmake -DDEAL_II_DIR=/path/to/deal.II .
```
in this directory to configure the problem.
You can switch between debug and release mode by calling either
```
make debug
```
or
```
make release
```
This command will generate two executables: one for the `coupled_laplace_problem` and one for the `dummy` participant.
```
make run
```
excutes the `coupled_laplace_problem`. In order to run the coupled simulation, execute
```
./dummy
```
in the same directory from another terminal window.

## Results

Presenting some results? Essentially a single picture?

## Further reading

* about [preCICE](https://github.com/precice/precice/) in general
* preCICE coupled deal.II codes in the [dealii-adapter](https://github.com/precice/dealii-adapter)
* some literature ?
