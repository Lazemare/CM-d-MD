"""
# CMdMD: Commute Map directed Molecular Dynamics
"""
module CMdMD

using PyCall
using Printf
using Statistics
using LinearAlgebra
using DelimitedFiles

# CMdMD infomations.
__version__ = "0.0.2"
__author__ = "Lazemare"
__path__ = "$(@__DIR__)"
__file__ = "$(@__DIR__)/CMdMD.jl"

# PyObject python module mdtraj
md = pyimport("mdtraj")

"""
# Parameter of CM-d-MD simulation system 

## MD parameter
- `namdconfigurefile::String`: NAMD configure file name.
- `initcoord::String`: PDB file with initial coordniate.
- `topo::String`: Topology file.
- `iter::Int=1`: Iteration number. 
- `steps::Int=1`: MD steps each iteration. 
- `dcdfreq::Int=1`: How offen to write DCD File. 
- `nproc::Int=1`: Run NAMD with nproc cpu cores. 
- `gpu::Union{String,Int}=-1`: Use GPU or not, `-1` for not use, 
or like `"0,1,2"`. 

## Order parameters
- `orderparam::String`: Use which order parameter as the input of commute map 
calculation. Acceptable values are: `rmsd`, `distance` and `position`. 
- `selec::String="all"`: Must be a complete structure, for example: 
- `all and not water`. When loading trajectory, only this part will be loaded. 
- `alignref::String="all"`: Used to align the whole structure, must be part of 
- `selec`, for example: `protein and resid 1 to 20` 
- `refframe::Int=0`: Use which frame as the reference frame to calculate RMSD 
or other order parameters. 

### RMSD
Note that if you use this option, you must define `alignref` and `refframe` too,
to align structures before calculating positions.
- `rmsdregionselection::String="all"`: We only calculate RMSD of atoms in this 
region, must be part of `select`, for example: `protein and resid 1 to 2` 
- `rmsdstructureselection::String="all"`: We only calculate RMSD of atoms in this
selection, acceptable values are: `all`, `backbone`, `sidechain` and `nonh`. 
- `rmsdprecision::String="atoms"`: Calculate in which level. Acceptable values are: 
`atoms` and `residues`. 

### Distance
- `disprecision::String="atoms"`: Use center coordniates of which structure as the
order parameter. Acceptable values are: `atoms` and `residues`. 
- `disregionselection::String="all"`: We only calculate position of structures in
this region, must be part of `select`, for example: `protein and resid 1 to 2`. 
- `disstructureselection::String="all"`: We only calculate position of structures 
in this selection, acceptable values are: `all`, `backbone`, `sidechain` and `nonh`. 

### Position
Note that if you use this option, you must define `alignref` and `refframe` too,
to align structures before calculating positions.
- `posprecision::String="atoms"`: Use center coordniates of which structure as the
order parameter. Acceptable values are: `atoms` and `residues`. Note that if your
structure is large, do not use this option.
- `posregionselection::String="all"`: We only calculate distances between structures 
in this region, must be part of `select`, for example: `protein and resid 1 to 2`. 
- `posstructureselection::String="all"`: We only calculate distances between structures
in this selection, acceptable values are: `all`, `backbone`, `sidechain` and `nonh`. 

## Stretch parameter
- `stretchprecision::String="atoms"`: When stretching, change coordniates of which 
structure together. Acceptable values are: `atoms` and `residues`.
- `stretchfactor::Float64=1.15`: Strength of stretching, should be larger than 1.0.

## Map parameter
- `freq::Int=1`: How offen we load the trajectory. Our tests suggest that 
`steps/dcdfreq/freq > 10000` is a good choise. 
- `commutemap::Bool=true`: Use commute map or normal TICA method. 
- `lag::Int=1`: Lag time of commute map calculation. 
"""
mutable struct MDsystem
    # MD param
    namdconfigurefile::String
    initcoord::String
    topo::String
    iter::Int
    steps::Int
    dcdfreq::Int
    nproc::Int
    gpu::Union{String,Int}
    # Order parameters
    orderparam::String
    selec::String
    refframe::Int
    alignref::String
    # RMSD parameter
    rmsdregionselection::String
    rmsdstructureselection::String
    rmsdprecision::String
    # Distance parameter
    disprecision::String
    disregionselection::String
    disstructureselection::String
    # Position parameter
    posprecision::String
    posregionselection::String
    posstructureselection::String
    # Stretch parameter
    stretchprecision::String
    stretchfactor::Float64
    # Map parameter
    method::String
    freq::Int
    commutemap::Bool
    lag::Int
    # Initializer
    MDsystem() = new()
end

# Include files
# IOs
include("$(@__DIR__)/io/traj.jl")
include("$(@__DIR__)/io/parser.jl")
include("$(@__DIR__)/io/namd_conf.jl")
# Core functions
include("$(@__DIR__)/core/cmap.jl")
include("$(@__DIR__)/core/pick.jl")
include("$(@__DIR__)/core/rmsd.jl")
include("$(@__DIR__)/core/position.jl")
include("$(@__DIR__)/core/distance.jl")
include("$(@__DIR__)/core/stretch.jl")
# Main calculation loop
include("$(@__DIR__)/md/run.jl")
# Helper functions
include("$(@__DIR__)/core/matrix.jl")

# export functions
export
    # main function
    runcmdmd, 
    # IOs
    loadtraj,
    convtraj,
    convtraj_all,
    inputparser,
    read_conf_namd,
    write_conf_namd,
    runnamd,
    # Core functions
    cmap,
    get_rmsd_org,
    get_rmsd_qcp,
    get_rmsd_by_residues,
    get_rmsd_by_atoms,
    get_distance_by_atoms,
    get_distance_by_residues,
    get_position_by_atoms,
    get_position_by_residues,
    pick_minmax,
    pick_pmf,
    pick_pmf2,
    pick_maxpmf,
    pick_maxpmf2,
    pick_pmf2plus,
    stretch_by_atoms,
    stretch_by_residues
end

using Main.CMdMD
println("CMdMD version: $(CMdMD.__version__)")
