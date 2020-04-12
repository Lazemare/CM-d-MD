"""
    get_rmsd_org(t::PyObject, alignref::String="all")

Get RMSD between each frame of the trajectory of the whole structure, 
using MDtraj RMSD engine. Return a n×n array, where n is the frame number.

# Arguments:
- `t::PyObject`: the trajectory.
- `alignref::String="all"`: structure used to align frames.
"""
function get_rmsd_org(t::PyObject, alignref::String="all")

    # RMSDs
    # To avoid the influence of large conformation change,
    # we align all frames to the first frame; and do not 
    # apply t.center_coordinates() to each frame, to reduce
    # the computational cost. NOTE that these methods have not 
    # been test carefully, and may be problematic.
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
    # calculate RMSD
    RMSD = zeros(t.n_frames,t.n_frames)
    for i = 0:(t.n_frames-1)
        RMSD[(i+1),:] = md.rmsd(t,t,i,atom_indices=nothing,precentered=true)
    end
    # Convert nm to Å and reduce the numerical error.
    RMSD = (RMSD + RMSD') * 5

return RMSD
end

"""
    get_rmsd_direct(t::PyObject,alignref::String="all",
    regselect::String="all",strucselect::String="all")

Get RMSD between each frame of the trajectory of the selected structure,
using direct RMSD algorithm.  Return a n×n array, where n is the frame number.

# Arguments:
- `t::PyObject`: the trajectory.
- `alignref::String="all"`: structure used to align frames.
- `regselect::String="all"`: region of the system to calculate RMSD.
- `strucselect::String="all"`: structure in given region (all, backbone, 
sidechain, or nonh) to calculate RMSD.

See also: [`get_rmsd_qcp`](@ref)
"""
function get_rmsd_direct(t::PyObject,alignref::String="all",
         regselect::String="all",strucselect::String="all")

    # calculate RMSD of part of the trajectory.
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
    # select
    regselection = topology.select("$regselect")
    t.restrict_atoms(regselection)
    # update topology
    topology = t.topology
    if strucselect == "nonh"
        strucselect = "not element H"
    end
    selection = topology.select("all and $strucselect")
    # calculate RMSD from coordinates
    coord = copy(t.xyz)
    RMSD = zeros(t.n_frames,t.n_frames)
    for i = 1:t.n_frames
        b = reshape(coord[i,:,:],(1,t.n_atoms,3))
        c = repeat(b,t.n_frames,1,1)
        xyz = coord .- c
        xyz = xyz .^ 2
        RMSD[i,:] = sum(sum(xyz[:,(selection.+1),:], dims=2), dims=3) ./ size(selection,1)
    end
    RMSD = (RMSD .^ 0.5 + RMSD'.^ 0.5) * 5

return RMSD
end

" Calculate RMSD with the QCP mrthod "
function _rmsd(A::Array, B::Array)

    G_A = sum(A.^2)
    G_B = sum(B.^2)
    M = B' * A
    K = zeros(4,4)
    K[1,1] = M[1,1] + M[2,2] + M[3,3]
    K[2,2] = M[1,1] - M[2,2] - M[3,3]
    K[3,3] = -M[1,1] + M[2,2] - M[3,3]
    K[4,4] = -M[1,1] - M[2,2] + M[3,3]
    K[1,2] = M[2,3] - M[3,2]
    K[1,3] = M[3,1] - M[1,3]
    K[1,4] = M[1,2] - M[2,1]
    K[2,1] = M[2,3] - M[3,2]
    K[2,3] = M[1,2] + M[2,1]
    K[2,4] = M[3,1] + M[1,3]
    K[3,1] = M[3,1] - M[1,3]
    K[3,2] = M[1,2] + M[2,1]
    K[3,4] = M[2,3] + M[3,2]
    K[4,1] = M[1,2] - M[2,1]
    K[4,2] = M[3,1] + M[1,3]
    K[4,3] = M[2,3] + M[3,2]
    eigenval = broadcast(real,eigvals(K))
    rmsd = (abs(G_A + G_B - 2 * maximum(eigenval)) / size(A,1)) ^ 0.5
    
return rmsd
end

"""
    get_rmsd_qcp(t::PyObject,alignref::String="all",
    regselect::String="all",strucselect::String="all";mp::Bool=true)

Get RMSD between each frame of the trajectory of the selected structure,
using QCP algorithm. Return a n×n array, where n is the frame number.

# Arguments:
- `t::PyObject`: the trajectory.
- `alignref::String="all"`: structure used to align frames.
- `regselect::String="all"`: region of the system to calculate RMSD.
- `strucselect::String="all"`: structure in given region (all, backbone, 
sidechain, or nonh) to calculate RMSD.
- `mp::Bool=true`: if use multi-processing.

See also: [`get_rmsd_direct`](@ref)
"""
function get_rmsd_qcp(t::PyObject,alignref::String="all",regselect::String="all",
         strucselect::String="all";mp::Bool=true)

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
    # Get topology
    topology = t.topology
    # select
    regselection = topology.select("$regselect")
    t.restrict_atoms(regselection)
    # update topology
    topology = t.topology
    if strucselect == "nonh"
        strucselect = "not element H"
    end
    selection = topology.select("all and $strucselect")
    # calculate RMSD from coordinates
    # get coord
    coord = copy(t.xyz[:,(selection.+1),:])
    RMSD = zeros(t.n_frames,t.n_frames)
    frames = copy(t.n_frames)
    if mp
        Threads.@threads for i = 1:(frames-1)
            for j = (i+1):frames
                RMSD[i,j] = _rmsd(coord[i,:,:],coord[j,:,:])
            end
        end
    else 
        for i = 1:(frames-1)
            for j = (i+1):frames
                RMSD[i,j] = _rmsd(coord[i,:,:],coord[j,:,:])
            end
        end
    end
    RMSD = (RMSD + RMSD') * 10

return RMSD
end

"""
    get_rmsd_by_residues(t::PyObject, refframe::Int=0, alignref::String="all",
    regselect::String="all", strucselect::String="all")

Get RMSD of residues belong to the selected structure in each frame of the
trajectory, using direct RMSD algorithm. Return a n×m array, where n is the 
frame number, m is the residue number.

# Arguments:
- `t::PyObject`: the trajectory.
- `alignref::String="all"`: structure used to align frames.
- `refframe::Int=0`: the reference frame used to calculate RMSD.
- `regselect::String="all"`: region of the system to calculate RMSD.
- `strucselect::String="all"`: structure in given region (all, backbone, 
sidechain, or nonh) to calculate RMSD.

See also: [`get_rmsd_by_atoms`](@ref)
"""
function get_rmsd_by_residues(t::PyObject,refframe::Int=0,alignref::String="all",
		 regselect::String="all",strucselect::String="all")

	# strucselect could be one of "all", "backbone", "sidechain" and "nonh"
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
    regselection = topology.select("$regselect")
    t.restrict_atoms(regselection)
    # update topology
    topology = t.topology
    if strucselect == "nonh"
        strucselect = "not element H"
    end
    # calculate RMSD from coordinates
    RMSD = zeros(t.n_frames,t.n_residues)
    b = reshape(t.xyz[1,:,:],(1,t.n_atoms,3))
    c = repeat(b,t.n_frames,1,1)
    xyz = t.xyz .- c
    xyz = xyz .^ 2
    flag = 0
    for residue in topology.residues
        index = residue.index
        selection = topology.select("all and resid $index and $strucselect")
        flag = flag + 1
        #for j = 1:t.n_frames
        #    RMSD[j,flag] = sum(xyz[j,(selection.+1),:]) / size(selection,1)
        #end
        RMSD[:,flag] = sum(sum(xyz[:,(selection.+1),:], dims=2), dims=3) ./ size(selection,1)
    end
    RMSD = (RMSD .^ 0.5) * 10

return RMSD
end

"""
    get_rmsd_by_atoms(t::PyObject, refframe::Int=0, alignref::String="all",
    regselect::String="all", strucselect::String="all")

Get RMSD of atoms belong to the selected structure in each frame of the
trajectory, using direct RMSD algorithm. Return a n×m array, where n is the 
frame number, m is the atom number.

# Arguments:
- `t::PyObject`: the trajectory.
- `alignref::String="all"`: structure used to align frames.
- `refframe::Int=0`: the reference frame used to calculate RMSD.
- `regselect::String="all"`: region of the system to calculate RMSD.
- `strucselect::String="all"`: structure in given region (all, backbone, 
sidechain, or nonh) to calculate RMSD.

See also: [`get_rmsd_by_residues`](@ref)
"""
function get_rmsd_by_atoms(t::PyObject, refframe::Int=0, alignref::String="all",
         regselect::String="all", strucselect::String="all")

    # Get RMSD of each atom in trajectory to frame 0.

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
    regselection = topology.select("$regselect and $strucselect")
    t.restrict_atoms(regselection)
    # update topology
    topology = t.topology
    # calculate RMSD from coordinates
    RMSD = zeros(t.n_frames,t.n_atoms)
    b = reshape(t.xyz[1,:,:],(1,t.n_atoms,3))
    c = repeat(b,t.n_frames,1,1)
    xyz = t.xyz .- c
    xyz = xyz .^ 2
    for i = 1:t.n_atoms
        for j = 1:t.n_frames
            RMSD[j,i] = sum(xyz[j,i,:])
        end
    end
    RMSD = (RMSD .^ 0.5) * 10

return RMSD
end
