"""
    stretch_by_atoms(t::PyObject, t_all::PyObject, selec::String, basin::Int,
    endpoint::Int, factor::Float64, alignref::String="all")

Stretch coordinates of all atoms towards given direction.

# Arguments:
- `t::PyObject`: trajectory contains structure selected by `MDsystem.select`.
- `t_all::PyObject`: origin trajectory.
- `selec::String`: must be the same as `MDsystem.select`.
- `basin::Int`: frame number of basin conformation (start from 1).
- `endpoint::Int`: frame number of endpoint conformation (start from 1).
- `factor::Float64`: stretching factor, must be larger than 1.
- `alignref::String="all"`: must be the same as `MDsystem.alignref`.

See also: [`stretch_by_residues`](@ref)
"""
function stretch_by_atoms(t::PyObject,t_all::PyObject,selec::String,basin::Int,
         endpoint::Int,factor::Float64,alignref::String="all")

    # Get topology
    topology = t.topology
    # Align
    t.center_coordinates()
    if alignref == "all"
        t.superpose(t,frame=0)
    else
        align_selection = topology.select("all and $alignref")
        t.superpose(t,frame=0,atom_indices=align_selection)
    end
    # t
    t_basin = get(t,basin,nothing)
    basin_xyz = t_basin.xyz[1,:,:]
    t_endpoint = get(t,endpoint,nothing)
    endpoint_xyz = t_endpoint.xyz[1,:,:]
    transM = (endpoint_xyz .- basin_xyz) * (factor - 1)
    transM_ = reshape(transM,(1,t.n_atoms,3))
    t_endpoint.xyz = t_endpoint.xyz .+ transM_
    # t_all
    # Align t_basin to its position in t_all
    t_tmp = get(t_all,endpoint,nothing)
    selection_tmp = topology.select("all and $selec")
    t_tmp.restrict_atoms(selection_tmp)
    t_endpoint.superpose(t_tmp)
    # change coordinates in t_all
    topology_all = t_all.topology
    selection = topology_all.select("$selec") .+ 1
    t_endpoint_all = get(t_all,endpoint,nothing)
    endpoint_all_xyz = copy(t_endpoint_all.xyz)
    endpoint_all_xyz[1,selection,:] = t_endpoint.xyz[1,:,:]
    t_endpoint_all.xyz = endpoint_all_xyz

return t_endpoint, t_endpoint_all
end

"""
    stretch_by_residues(t::PyObject, t_all::PyObject, selec::String, basin::Int, 
    endpoint::Int, factor::Float64, alignref::String="all")

Stretch coordinates of residues towards given direction, 
by the center coordinates of each residue.

# Arguments:
- `t::PyObject`: trajectory contains structure selected by `MDsystem.select`.
- `t_all::PyObject`: origin trajectory.
- `selec::String`: must be the same as `MDsystem.select`.
- `basin::Int`: frame number of basin conformation (start from 1).
- `endpoint::Int`: frame number of endpoint conformation (start from 1).
- `factor::Float64`: stretching factor, must be larger than 1.
- `alignref::String="all"`: must be the same as `MDsystem.alignref`.

See also: [`stretch_by_atoms`](@ref)
"""
function stretch_by_residues(t::PyObject,t_all::PyObject,selec::String,basin::Int,
         endpoint::Int,factor::Float64,alignref::String="all")

    # Get topology
    topology = t.topology
    # Align
    t.center_coordinates()
    if alignref == "all"
        t.superpose(t,frame=0)
    else
        align_selection = topology.select("all and $alignref")
        t.superpose(t,frame=0,atom_indices=align_selection)
    end
    # t
    t_basin = get(t,basin,nothing)
    basin_xyz = t_basin.xyz[1,:,:]
    basin_pos = zeros(1,3)
    t_endpoint = get(t,endpoint,nothing)
    endpoint_xyz = t_endpoint.xyz[1,:,:]
    endpoint_pos = zeros(1,3)
    transM = zeros(t.n_atoms,3)
    # Stretch.
    flag = 1
    for residue in topology.residues
        index = residue.index
        selection = topology.select("all and resid $index")
        basin_pos = sum(basin_xyz[(selection.+1),:],dims=1) ./ size(selection,1)
        endpoint_pos = sum(endpoint_xyz[(selection.+1),:],dims=1) ./ size(selection,1)
        (transM[(selection.+1),:] = 
            repeat(((endpoint_pos .- basin_pos) .* (factor - 1)),size(selection,1),1))
        flag = flag + 1
    end
    transM_ = reshape(transM,(1,t.n_atoms,3))
    t_endpoint.xyz = t_endpoint.xyz .+ transM_
    # t_all
    # Align t_basin to its position in t_all
    t_tmp = get(t_all,endpoint,nothing)
    selection_tmp = topology.select("all and $selec")
    t_tmp.restrict_atoms(selection_tmp)
    t_endpoint.superpose(t_tmp)
    # change coordinates in t_all
    topology_all = t_all.topology
    selection = topology_all.select("$selec") .+ 1
    t_endpoint_all = get(t_all,endpoint,nothing)
    endpoint_all_xyz = copy(t_endpoint_all.xyz)
    endpoint_all_xyz[1,selection,:] = t_endpoint.xyz[1,:,:]
    t_endpoint_all.xyz = endpoint_all_xyz

return t_endpoint, t_endpoint_all
end
