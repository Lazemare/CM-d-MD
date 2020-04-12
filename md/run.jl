""" 
    runcmdmd(configfile::String)

Run CM-d-MD simulation.
"""
function runcmdmd(configfile::String)

    system = inputparser(configfile)

    for i = 1:system.iter

        mkdir("ITER_$i")
        if i == 1
            cp("$(system.initcoord)","ITER_1/ITER_1.pdb")
            cp("$(system.topo)","ITER_1/ITER_1.psf")
            println("# ---------- ITER 1 ---------- #")
        else
            cp("ITER_$(i-1)/endpoint_all.pdb","ITER_$i/ITER_$i.pdb")
            cp("$(system.topo)","ITER_$i/ITER_$i.psf")
        end

        cd("ITER_$i")
        @label runNAMD
        write_conf_namd("../$(system.namdconfigurefile)","ITER_$i",
                        "ITER_$i.pdb","ITER_$i.psf",system.steps,system.dcdfreq)
        println("Runing NAMD job ITER_$i ... ")
        try
            TIME = @timed runnamd("ITER_$i.conf",system.nproc,system.gpu)
        catch 
            @error "NAMD terminaled abnormally! Reason: BAD STRUCTURE!"
            @warn "If above error shows up too offen, you should check your params." 
            println("changing factor ...")
            rm("old",recursive=true,force=true)
            mkdir("old")
            cp("NAMD.log","./old/NAMD.log")
            cp("ERROR.log","./old/ERROR.log")
            rm("ITER_$i.pdb")
            cp("../ITER_$(i-1)/ITER_$(i-1).pdb","./ITER_$i.pdb")
            @goto runNAMD 
        end
        try 
            open("ITER_$i.dcd","r") do f
            end
        catch 
            @error "NAMD terminaled abnormally! Reason: BAD STRUCTURE!"
            @warn "If above error shows up too offen, you should check your params." 
            println("changing factor ...")
            rm("old",recursive=true,force=true)
            mkdir("old")
            cp("NAMD.log","./old/NAMD.log")
            cp("ERROR.log","./old/ERROR.log")
            rm("ITER_$i.pdb")
            cp("../ITER_$(i-1)/ITER_$(i-1).pdb","./ITER_$i.pdb")
            @goto runNAMD 
        end
        # Dump first frame as initial coordinate.
        if i == 1
            initcoord_dcd = md.load("ITER_1.dcd",top="ITER_1.psf")
            initcoord_dcd = get(initcoord_dcd,1,nothing)
            initcoord_dcd.save_dcd("../initcoord.dcd")
        end
        # Load trajectory.
        """ The last frame of `traj` is the initial conformation. """
        traj_all = loadtraj("ITER_$i.dcd","ITER_$i.psf","all",system.freq)
        initframe_all = loadtraj("../initcoord.dcd","ITER_$i.psf","all",1)
        traj_all = md.join([traj_all, initframe_all])
        # make some tests.
        if (traj_all.n_frames < (system.steps/system.dcdfreq-1))
            @error "NAMD terminaled abnormally! Reason: SIMULATION BLOW UP!"
            @warn "If above error shows up too offen, you should check your params."
            println("changing factor ...")
            rm("old",recursive=true,force=true)
            mkdir("old")
            cp("NAMD.log","./old/NAMD.log")
            cp("ERROR.log","./old/ERROR.log")
            rm("ITER_$i.pdb")
            cp("../ITER_$(i-1)/ITER_$(i-1).pdb","./ITER_$i.pdb")
            @goto runNAMD 
        else 
            println("NAMD terminaled normally.")
            println("Job ITER_$i takes $(TIME[2]) seconds to complete.")
        end
        traj = loadtraj("ITER_$i.dcd","ITER_$i.psf",system.selec,system.freq)
        initframe = loadtraj("../initcoord.dcd","ITER_$i.psf",system.selec,1)
        traj = md.join([traj, initframe])
        println("Start to calculate Order Parameters and Commute Map ...")
        println("Calculating Order Parameters ...")
        TIME = @timed (
        if system.orderparam == "rmsd"
            if system.rmsdprecision == "residues"
                param = get_rmsd_by_residues(traj,system.refframe,system.alignref,
                           system.rmsdregionselection,system.rmsdstructureselection)
            elseif system.rmsdprecision == "atoms"
                param = get_rmsd_by_atoms(traj,system.refframe,system.alignref,
                        system.rmsdregionselection,system.rmsdstructureselection)
            end
        elseif system.orderparam == "distance"
            if system.disprecision == "atoms"
                param = get_distance_by_atoms(traj,system.disregionselection,
                        system.disstructureselection)[1]
            elseif system.disprecision == "residues"
                param = get_distance_by_residues(traj,system.disregionselection,
                        system.disstructureselection)[1]
            end
        elseif system.orderparam == "position"
            if system.posprecision == "atoms"
                param = get_position_by_atoms(traj,system.refframe,system.alignref,
                        system.posregionselection,system.posstructureselection)
            elseif system.posprecision == "residues"
                param = get_position_by_residues(traj,system.refframe,system.alignref,
                        system.posregionselection,system.posstructureselection)
            end
        end )
        println("Calculation of Order Parameters takes $(TIME[2]) seconds to complete.")
        println("Calculating Commute Map ...")
        println("Lag Time: $(system.lag)")
        println("Size of order parameters: $(size(param))")
        TIME = @timed rts, eigenvec, vector = cmap(param,system.lag,commute_map=system.commutemap)
        println("Calculation of Commute Map takes $(TIME[2]) seconds to complete.")
        println("There are $(size(vector,1)) points in the map.")
        @printf "%19s    %6s\n" "min" "max"
        @printf "%s %10.4f   %4.4f\n" " 1st coord" minimum(vector[1,:]) maximum(vector[1,:])
        @printf "%s %10.4f   %4.4f\n" " 2nd coord" minimum(vector[2,:]) maximum(vector[2,:])
        @printf "%s %10.4f   %4.4f\n" " 3rd coord" minimum(vector[3,:]) maximum(vector[3,:])
        @printf "%s %10.4f   %4.4f\n" " 4th coord" minimum(vector[4,:]) maximum(vector[4,:])
        @printf "%s %10.4f   %4.4f\n" " 5th coord" minimum(vector[5,:]) maximum(vector[5,:])
        @printf "%s %10.4f   %4.4f\n" " 6th coord" minimum(vector[6,:]) maximum(vector[6,:])
        @printf "%s %10.4f   %4.4f\n" " 7th coord" minimum(vector[7,:]) maximum(vector[7,:])
        @printf "%s %10.4f   %4.4f\n" " 8th coord" minimum(vector[8,:]) maximum(vector[8,:])
        @printf "%s %10.4f   %4.4f\n" " 5th coord" minimum(vector[9,:]) maximum(vector[9,:])
        @printf "%s %10.4f   %4.4f\n" "10th coord" minimum(vector[10,:]) maximum(vector[10,:])

        if i == 1

            traj = loadtraj("ITER_$i.dcd","ITER_$i.psf",system.selec,system.freq)
            initframe = loadtraj("../initcoord.dcd","ITER_$i.psf",system.selec,1)
            traj = md.join([traj, initframe])
            TIME = @timed basin, endpoint = pick_pmf2plus(vector[2:(size(vector,1)-1),:])
            println("pickPMF takes $(TIME[2]) seconds to complete.")
            println("frame $(system.freq*basin[1]) contains lowest free energy.")
            get(traj,basin[1],nothing).save_pdb("basin.pdb")
            println("file basin.pdb has been written!")
            println("frame $(system.freq*endpoint[1]) contains highest free energy.")
            println("Stretching ...")
            fa = system.stretchfactor
            println("Using stretch factor: $fa")
            TIME = @timed (if system.stretchprecision == "atoms"
                endpoint_,endpoint_all = stretch_by_atoms(traj,traj_all,system.selec,basin[1],
                                         endpoint[1],fa,system.alignref)
            elseif system.stretchprecision == "residues"
                endpoint_,endpoint_all = stretch_by_residues(traj,traj_all,system.selec,basin[1],
                                         endpoint[1],fa,system.alignref)
            end)
            println("Stretching takes $(TIME[2]) seconds to complete.")
            endpoint_.save_pdb("endpoint.pdb")
            println("file endpoint.pdb has been written!")
            endpoint_all.save_pdb("endpoint_all.pdb")
            println("file endpoint_all.pdb has been written!")
            open("../traj_info","a") do f
                @printf f "%s:%d  %d \n" "ITER" i Int(system.freq*endpoint[1])
            end
            println("# ---------- ITER 2 ---------- #")

        else

            traj = loadtraj("ITER_$i.dcd","ITER_$i.psf",system.selec,system.freq)
            initframe = loadtraj("../initcoord.dcd","ITER_$i.psf",system.selec,1)
            traj = md.join([traj, initframe])
            write_flag = 0
            # Check propagate directions.
            TIME = @timed basin, endpoint = pick_pmf2plus(vector[2:(size(vector,1)-1),:])
            println("pickPMF takes $(TIME[2]) seconds to complete.")
            dir1=sign(maximum(vector[:,1])+minimum(vector[:,1])-2*vector[size(vector,1),1])
            dir2=sign(maximum(vector[:,2])+minimum(vector[:,2])-2*vector[size(vector,1),2])
            # Make sure that we picked a point in the right direction.
            #dir_basin1=sign(maximum(vector[1:(size(vector,1)-1),1])+
            #                minimum(vector[1:(size(vector,1)-1),1])-2*vector[(basin[1]+1),1])
            #dir_basin2=sign(maximum(vector[1:(size(vector,1)-1),1])+
            #                minimum(vector[1:(size(vector,1)-1),1])-2*vector[(basin[1]+1),1])
            #if dir1 == dir_basin1
            #    restrict = 1
            #    for j=1:size(vector,1)
            #        if (((dir_basin1 == 1 && vector[j,1] > vector[(basin[1]+1),1]) || 
            #            (dir_basin1 == -1 && vector[j,1] < vector[(basin[1]+1),1]) &&
            #            ((dir1 == 1 && vector[j,1] > vector[size(vector,1),1]) ||
            #            (dir1 == -1 && vector[j,1] < vector[size(vector,1),1]))))
            #            restrict = [restrict i]
            #        end
            #    end
            #    if size(restrict,2) >= 2
            #        restrict = restrict[1:(size(restrict,2)-1)]
            #        basin, endpoint = pick_pmf2(vector[1:(size(vector,1)-1),:],1,2,restrict)
            #    else
            #        @warn "Bad dynamics! This CMdMD step maybe useless! Using basin as endpoint!"
            #        endpoint = basin
            #    end
            #if dir1 != dir_basin1
            #    println("dir1 != dir_basin1, using basin as endpoint.")
            #    if dir1 == 1
            #        order = sortperm(vector[(basin.+1),1],rev=true)
            #        endpoint = basin[order]
            #    else 
            #        order = sortperm(vector[(basin.+1),1])
            #        endpoint = basin[order]
            #    end
            #end
            println("frame $(system.freq*basin[1]) contains lowest free energy.")
            get(traj,basin[1],nothing).save_pdb("basin.pdb")
            println("file basin.pdb has been written!")
            println("frame $(system.freq*endpoint[1]) contains highest free energy.")
            println("Stretching ...")
            fa = system.stretchfactor
            println("Using stretch factor: $fa")
            TIME = @timed (if system.stretchprecision == "atoms"
                endpoint_,endpoint_all = stretch_by_atoms(traj,traj_all,system.selec,basin[1],
                                         endpoint[1],fa,system.alignref)
            elseif system.stretchprecision == "residues"
                endpoint_,endpoint_all = stretch_by_residues(traj,traj_all,system.selec,basin[1],
                                         endpoint[1],fa,system.alignref)
            end)
            println("Stretching takes $(TIME[2]) seconds to complete.")
            endpoint_.save_pdb("endpoint.pdb")
            println("file endpoint.pdb has been written!")
            endpoint_all.save_pdb("endpoint_all.pdb")
            println("file endpoint_all.pdb has been written!")
            # traj_basin_old = md.load("../ITER_$(i-1)/basin.pdb")
            # traj_basin = md.load("basin.pdb")
            # traj_tmp = md.join([traj_basin_old,traj_basin])
            # rmsd_basin = get_rmsd_qcp(traj_tmp,system.alignref,
            #             system.rmsdregionselection,system.rmsdstructureselection)
            # println("RMSD between last and current basin is: $(rmsd_basin[2])")
            # Hop
            ################## ! ! ! TODO ! ! ! ####################
            # rmsd_avg = (sum(RMSD[:,1].^2)./size(RMSD[:,1],1)) .^ 0.5
            ################## ! ! ! TODO ! ! ! ####################
            # if rmsd_basin[1] < (rmsd_avg / 5)
            #    println("Same basin! Use current basin and endpoint!")
            # elseif rmsd_basin[1] > (rmsd_avg / 5)
            #    println("Different basin! Will determine if the endpoints are the same!")
            #    traj_end_old = md.load("../ITER_$(i-1)/endpoint.pdb")
            #    traj_end = md.load("endpoint.pdb")
            #    get(traj_all,endpoint[1],nothing).save_pdb("endpoint_all.pdb")
            #    rmsd_end = md.rmsd(traj_basin,traj_basin_old)
            #    if rmsd_end[1] > (rmsd_avg / 5)
            #        println("Different endpoint! Use current basin and endpoint!")
            #    elseif rmsd_end[1] < (rmsd_avg / 5)
            #        println("Same endpoint! This CMdMD step maybe useless!")
            #    end
            #end
            
            open("../traj_info","a") do f
                if write_flag == 0
                    @printf f "%s:%d  %d \n" "ITER" i max(Int(system.freq*endpoint[1]),1)
                elseif write_flag == 1
                    @printf f "%s:%d  %d \n" "ITER" i max(Int(system.freq*endpoint_new[1]),1)  
                end
            end
            if i != system.iter
                println("# ---------- ITER $(i+1) ---------- #")
            end
        end
        cd("../")
    end
    println("# ---------- CMdMD terminaled normally. ---------- #")

return nothing
end
