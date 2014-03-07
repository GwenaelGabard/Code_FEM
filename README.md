### Code FEM

This Matlab program provides a platform to facilitate the implementation of
Finite Element and Boundary Element formulations. It is intended as a research
tool.


## Structure

Three levels can be considered.

1. Using an existing finite element formulation to solve a given problem.

2. Implementation and validating a new finite element formulation. This is the
main purpose of this code. It is designed to make it (relatively) easy to
go from a variational formulation on paper to a working implementation that can
be validated and assessed. The different finite element formulations, or boundary
element formulations, available are located in the Library directory.

3. The low level and generic tasks required by all finite element models, such
as allocating degrees of freedom, the assembly of element matrices in a global
system. All these aspects are provided by functions and scripts in the Core
directory.


## Documentation

Further details will be made available in the Documentation directory.


## Examples

See the Examples directory.


## Some background

The first version of this code was originally developed at the Acoustics Group
at the University of Technology of Compiegne, France, back in 2000. Various
revisions were further develop, and the code was used for several research
projects. It has since been used in other research groups, leading to several
different version being developed independently. This repository is an attempt
to bring some of these features together.
