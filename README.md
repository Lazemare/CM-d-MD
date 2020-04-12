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

   The CM-d-MD configure file contains two parts basically. 
   The first part defines some MD simulation parameters:

   `NAMDConfigureFile`: name of NAMD configure file.
   
   `Iterations`: the iteration number of the CM-d-MD simulation.
   
   `Steps`: MD steps per iteration.
   
   `DCDFreq`: frequency of writing DCD file.
   
   `Nproc`: number of CPU cores to use in the simulation.
   
   `GPU`: GPU id to use in the simulation. -1 for do not use GPU.

   The second part defines parameters about the commute map:

   `Selection`: the complete structure contains the region you would like to perform accelerated sampling. All selections in CM-d-MD use MDtraj syntax. For example, if we have a protein solvated by explicit water, and we would like to sample the conformation change of the protein, the `Selection` key word should be set as `"all and protein"`.
   
   `OrderParameter`: the input data type of the commute map calculation. Acceptable values are: `RMSD`, `Distance` and `Position`. 
   
   `AlignReference`: use which structure to align frames. It should be part of the structure selected by `Selection`. For example, `protein and resid 0 to 10`.
   
   `ReferenceFrame`: align frames to which frame in the trajectory.
   
   `StretchPrecision`: stretch the structure by transforming coordinates of each atom or each residue. Acceptable values are: `atoms` and `residues`. 
   
   `StretchFactor`: the strength of structure stretching, must be larger than 1.
   
   `RMSDRegionSelection`, `DisRegionSelection`, `PosRegionSelection`: calculate RMSD/distance/position of which region in structure selected by `Selection`. It should be part of the structure selected by `Selection`. For example, `protein and resid 10 to 20`.
   
   `RMSDStructureSelection`, `DisStructureSelection`, `PosStructureSelection`: calculate RMSD/distance/position of which structure in selected region. Acceptable values are: `all`, `backbone`, `sidechain`, or `nonh`.
   
   `RMSDPrecision`, `DisPrecision`, `PosPrecision`: calculate RMSD/distance/position of atoms or residues in selected structure. Acceptable values are: `atoms` and `residues`.
   
   `CommuteMap`: perform dimension reduction with commute map algorithm or TICA algorithm.
   
   `LagTime`: lag time of commute map or TICA calculation.
   
   `Frequency`: frequency of loading the trajectories.

## Problems

The Koopman reweighting scheme do not work by now, please use default reweighting scheme.
