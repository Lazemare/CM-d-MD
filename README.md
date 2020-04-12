# CMdMD: Commute Map directed Molecular Dynamics


## What's this?

CM-d-MD is a package used for performing the so-called Commute Map directed Molecular Dynamics accelerated sampling method. This idea is based on the  Diffusion-Map-Directed Molecular Dynamics (https://dx.doi.org/10.1021/jp401911h , https://dx.doi.org/10.1039/c3cp54520b) and the Commute Map (https://dx.doi.org/10.1021/acs.jctc.6b00762).



## Dependences

Julia packages:

`Julia > 1.2.0` `PyCall`

Others:

`NAMD > 2.9` `MDtraj = 1.9.3`

## Installation

CM-d-MD could only run under Linux operating system.

1. Install NAMD (http://www.ks.uiuc.edu/Research/namd/) and MDtraj using pip or conda:

   ```
   pip install mdtraj
   ```

   or

   ```
   conda install mdtraj
   ```

   Make sure when typing `namd2` command under terminal, you could start the NAMD program.

2. Install Julia and PyCall package. You could make a check by typing this command in the Julia REPL, which will give you the version of MDtraj:

   ```julia
   using PyCall
   md = pyimport("mdtraj")
   println(md.version.version)
   ```

3. Clone this package to folder called `CMdMD`, then `cd` into `CMdMD`, type this command:

   ```julia
   include("CMdMD.jl")
   ```

   This will load the CM-d-MD codes and show the CM-d-MD version.

## Examples

In the example folder we provide some examples, check the `README.md` there for more information about how to work with CM-d-MD.

## Run

1. How to run CM-d-MD

   ```julia
   include("/path/to/CM-d-MD/CMdMD.jl")
   cmdmdrun("input_file")
   ```

2. Configure file

## Problems

The Koopman reweighting scheme do not work by now, please use default reweighting scheme.
