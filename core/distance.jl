"""
    get_distance_by_atoms(t::PyObject, regselect::String="all",
    strucselect::String="all"; mp=true)

Get distance between atoms belong to the selected structure in each frame of the
trajectory,. Return a n×m array, where n is the frame number, m is the atom number.

# Arguments:
- `t::PyObject`: the trajectory.
- `regselect::String="all"`: region of the system to calculate RMSD.
- `strucselect::String="all"`: structure in given region (all, backbone,
sidechain, or nonh) to calculate RMSD.
- `mp::Bool=true`: if use multi-processing.

NOTE: THIS COULD BE REAL SLOW AND REQUIRE VERY MUCH MEMORY!!!

See also: [`get_distance_by_residues`](@ref)
"""
function get_distance_by_atoms(t::PyObject,regselect::String="all",
         strucselect::String="all";mp=true)

    # Get topology
    topology = t.topology
    # Align
    t.center_coordinates()
    # select
    if strucselect == "nonh"
        strucselect = "not element H"
    end
    selection = topology.select("all and $strucselect and $regselect")
    t.restrict_atoms(selection)
    # update topology
    topology = t.topology
    # Calculate distance
    n_frames = copy(t.n_frames)
    n_atoms  = copy(t.n_atoms)
    order = zeros(sum(range(1,stop=(n_atoms-1))),2)
    dis = zeros(n_frames,sum(range(1,stop=(n_atoms-1))))
    xyz = copy(t.xyz)
    if mp
        Threads.@threads for i = 1:n_atoms
            for j = (i+1):n_atoms
                order[sign(i-1)*sum(range(((n_atoms+1)-i),stop=(n_atoms-1))) + (j - i),1] = i
                order[sign(i-1)*sum(range(((n_atoms+1)-i),stop=(n_atoms-1))) + (j - i),2] = j
                dis[:,sign(i-1)*sum(range(((n_atoms+1)-i),stop=(n_atoms-1))) + (j - i)] =
                sum((xyz[:,i,:].-xyz[:,j,:]).^2,dims=2).^0.5
            end
        end
    else
        for i = 1:n_atoms
            for j = (i+1):n_atoms
                order[sign(i-1)*sum(range(((n_atoms+1)-i),stop=(n_atoms-1))) + (j - i),1] = i
                order[sign(i-1)*sum(range(((n_atoms+1)-i),stop=(n_atoms-1))) + (j - i),2] = j
                column = sign(i-1)*sum(range(((n_atoms+1)-i),stop=(n_atoms-1))) + (j - i)
                dis[:,column] = sum((xyz[:,i,:].-xyz[:,j,:]).^2,dims=2).^0.5
            end
        end
    end
    dis = dis .* 10

return dis, Int.(order)
end

"""
    get_distance_by_residues(t::PyObject, regselect::String="all",
    strucselect::String="all"; mp=true)

Get distance between residues belong to the selected structure in each frame of the
trajectory,. Return a n×m array, where n is the frame number, m is the residue number.

# Arguments:
- `t::PyObject`: the trajectory.
- `regselect::String="all"`: region of the system to calculate RMSD.
- `strucselect::String="all"`: structure in given region (all, backbone,
sidechain, or nonh) to calculate RMSD.
- `mp::Bool=true`: if use multi-processing.

See also: [`get_distance_by_atoms`](@ref)
"""
function get_distance_by_residues(t::PyObject,regselect::String="all",
         strucselect::String="all";mp=true)

    # Get topology
    topology = t.topology
    # Align
    t.center_coordinates()
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
    order = zeros(sum(range(1,stop=(n_residues-1))),2)
    dis = zeros(n_frames,sum(range(1,stop=(n_residues-1))))
    xyz = copy(t.xyz)
    pos = zeros(n_frames,n_residues,3)
    flag = 1
    # Get center coordinates of each residue.
    for residue in topology.residues
        index = residue.index
        selection = topology.select("all and resid $index and $strucselect")
        pos[:,flag,:] = sum(xyz[:,(selection.+1),:],dims=2) ./ size(selection,1)
        flag = flag + 1
    end
    # Get distance between each residue.
    if mp
        Threads.@threads for i = 1:n_residues
            for j = (i+1):n_residues
                order[sign(i-1)*sum(range(((n_residues+1)-i),stop=(n_residues-1))) + (j - i),1] = i
                order[sign(i-1)*sum(range(((n_residues+1)-i),stop=(n_residues-1))) + (j - i),2] = j 
                dis[:,sign(i-1)*sum(range(((n_residues+1)-i),stop=(n_residues-1))) + (j - i)] = 
                sum((pos[:,i,:].-pos[:,j,:]).^2,dims=2).^0.5
            end
        end
    else
        for i = 1:n_residues
            for j = (i+1):n_residues
                order[sign(i-1)*sum(range(((n_residues+1)-i),stop=(n_residues-1))) + (j - i),1] = i
                order[sign(i-1)*sum(range(((n_residues+1)-i),stop=(n_residues-1))) + (j - i),2] = j 
                column = sign(i-1)*sum(range(((n_residues+1)-i),stop=(n_residues-1))) + (j - i)
                dis[:,column] = sum((pos[:,i,:].-pos[:,j,:]).^2,dims=2).^0.5
            end
        end
    end
    dis = dis .* 10

return dis, Int.(order)
end
