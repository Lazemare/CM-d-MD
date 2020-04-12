"""
    get_position_by_atoms(t::PyObject, refframe::Int=0, alignref::String="all",
    regselect::String="all", strucselect::String="all"; mp=true)

Get positions of atoms in the trajectory of the selected structure. 
Return a 3n×m array, where n is the frame number, m is the atom number:

`atom_1_x_frame_0 atom_1_y_frame_0 ... atom_n_x_frame_0   ... `
`atom_1_x_frame_1 atom_1_y_frame_1 ... atom_n+1_x_frame_0 ... `

# Arguments:
- `t::PyObject`: the trajectory.
- `refframe::Int=0`: the reference frame used to align frames.
- `alignref::String="all"`: structure used to align frames.
- `regselect::String="all"`: region of the system to get positions.
- `strucselect::String="all"`: structure in given region (all, backbone, 
sidechain, or nonh) to get positions.
- `mp::Bool=true`: if use multi-processing.

See also: [`get_position_by_residues`](@ref)
"""
function get_position_by_atoms(t::PyObject,refframe::Int=0,alignref::String="all",
         regselect::String="all",strucselect::String="all";mp=true)

    # Get topology
    topology = t.topology
    # Align
    t.center_coordinates()
    if alignref == "all"
        t.superpose(t,frame=refframe)
    else
        align_selection = topology.select("all and $alignref")
        t.superpose(t,frame=refframe,atom_indices=align_selection)
    end
    # select
    if strucselect == "nonh"
        strucselect = "not element H"
    end
    selection = topology.select("all and $strucselect and $regselect")
    t.restrict_atoms(selection)
    # update topology
    topology = t.topology
    # Get position
    n_frames = copy(t.n_frames)
    n_atoms  = copy(t.n_atoms)
    pos = zeros(n_frames,(n_atoms*3))
    xyz = copy(t.xyz)
    if mp
        Threads.@threads for i = 0:(n_atoms-1)
            pos[:,(3*i+1):(3*i+3)] = xyz[:,(i+1),:]
        end
    else 
        for i = 0:(n_atoms-1)
            pos[:,(3*i+1):(3*i+3)] = xyz[:,(i+1),:]
        end
    end
    pos = pos .* 10

return pos
end

"""
    get_position_by_residues(t::PyObject, refframe::Int=0, alignref::String="all",
    regselect::String="all", strucselect::String="all"; mp=true)

Get positions of residues in the trajectory of the selected structure. 
Return a 3n×m array, where n is the frame number, m is the residue number:

`res_1_x_frame_0 res_1_y_frame_0 ... res_n_x_frame_0   ... `
`res_1_x_frame_1 res_1_y_frame_1 ... res_n+1_x_frame_0 ... `

# Arguments:
- `t::PyObject`: the trajectory.
- `refframe::Int=0`: the reference frame used to align frames.
- `alignref::String="all"`: structure used to align frames.
- `regselect::String="all"`: region of the system to get positions.
- `strucselect::String="all"`: structure in given region (all, backbone, 
sidechain, or nonh) to get positions.
- `mp::Bool=true`: if use multi-processing.

See also: [`get_position_by_atoms`](@ref)
"""
function get_position_by_residues(t::PyObject,refframe::Int=0,alignref::String="all",
         regselect::String="all",strucselect::String="all";mp=true)

    # Get topology
    topology = t.topology
    # Align
    t.center_coordinates()
    if alignref == "all"
        t.superpose(t,frame=refframe)
    else
        align_selection = topology.select("all and $alignref")
        t.superpose(t,frame=refframe,atom_indices=align_selection)
    end
    # select
    if strucselect == "nonh"
        strucselect = "not element H"
    end
    selection = topology.select("all and $regselect")
    t.restrict_atoms(selection)
    # update topology
    topology = t.topology
    # Calculate distance
    n_frames = copy(t.n_frames)
    n_residues = copy(t.n_residues)
    dis = zeros(n_frames,sum(range(1,stop=(n_residues-1))))
    xyz = copy(t.xyz)
    pos3 = zeros(n_frames,n_residues,3)
    flag = 1
    # Get center coordinates of each residue.
    for residue in topology.residues
        index = residue.index
        selection = topology.select("all and resid $index and $strucselect")
        pos3[:,flag,:] = sum(xyz[:,(selection.+1),:],dims=2) ./ size(selection,1)
        flag = flag + 1
    end
    # re-formate
    pos = zeros(n_frames,(n_residues*3))
    if mp
        Threads.@threads for i = 0:(n_residues-1)
            pos[:,(3*i+1):(3*i+3)] = pos3[:,(i+1),:]
        end
    else 
        for i = 0:(n_residues-1)
            pos[:,(3*i+1):(3*i+3)] = pos3[:,(i+1),:]
        end
    end
    pos = pos .* 10

return pos
end
