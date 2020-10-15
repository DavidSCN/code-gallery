Laplace equation coupled to an external simulation program
------------------------------------------
## Overview

preCICE allows to couple deal.II to external simulation software, such as OpenFOAM, SU2, or CalculiX. To keep dependencies of this example minimal we couple deal.II to an external c++ program, which provides a time varying boundary condition. The deal.II code consists mainly of the `step-4` tutorial program, where a simple Laplace problem is solved. Coupling with preCICE is usually carried out along surfaces in order to apply a Dirichlet-Neumann coupling between two domains (volume coupling is also possible). For the sake of simplicity, we couple here one side of our quadrilateral domain unidirectional with a c++ program, which generates a parabolic boundary profile with a time varying amplitude. The boundary values are consequntly used in the Laplace solver as Dirichlet boundary condition.

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
This command will generate two executables: one for the `coupled_laplace_problem` and one for the `fancy_boundary_condition` participant.
```
make run
```
excutes the `coupled_laplace_problem`. In order to run the coupled simulation, execute
```
./fancy_boundary_condition
```
in the same directory from another terminal window.

## Results

TODO: Add some results

## Further reading

* A complete overview over the preCICE project can be found on the [official preCICE webpage](https://www.precice.org/)
* The [source code](https://github.com/precice/precice/) of preCICE is hosted on Github
* If you are in particular interested in preCICE coupled deal.II codes (currently dealing with solid mechanics), have a look in the [dealii-adapter repository](https://github.com/precice/dealii-adapter)
* [The preCICE reference paper](https://www.sciencedirect.com/science/article/abs/pii/S0045793016300974)
