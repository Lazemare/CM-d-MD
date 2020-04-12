"""
    pick_minmax(vec::Array{Float64}, index::Int=1)

Return the index of elements with max and min value of input array.
"""
function pick_minmax(vec::Array{Float64}, index::Int=1)

    # Determine the index of elements with max and min value of input array.
    order = sortperm(vec[:,index])
    minmax = [order[1] order[size(order,1)]]
    minmax = minmax .- 1

return minmax
end

"""
    pmf1(vec::Array{Float64}, index::Int=1; nbins::Any="auto")

Calculate 1-D pmf.

See also: [`pmf2`, `pmf2_gs`, `pmf2plus`](@ref)
"""
function pmf1(vec::Array{Float64},index::Int=1;nbins::Any="auto")

    if nbins == "auto"
        nbins = size(vec,1) .^ 0.5 |> floor |> Int
    elseif typeof(nbins) == String
        error("Acceptable value of nbins is one of \"auto\" or Int (must be bigger than 1)")
    end
    minbin = minimum(vec[:,index])
    maxbin = maximum(vec[:,index])
    bin = (maxbin - minbin) / nbins
    bins = Array{Float64,1}(UndefInitializer(),nbins)
    bins = bins .- bins
    for i = 1:nbins
        bins[i] = bins[i] + (i-1)*bin
    end
    bins = bins .+ minbin
    bins[nbins] = maxbin
    pmf = fill(1.0,nbins,1)
    for i = 1:size(vec,1)
        tmp = searchsortedfirst(bins, vec[i,index])
        pmf[tmp] = pmf[tmp] + 1
    end
    pmf = pmf ./ size(vec,1)
    pmf = broadcast(log,pmf) .* -1
    pmf = pmf .- minimum(pmf)

return pmf, nbins
end

"""
    pmf2(vec::Array{Float64}, index1::Int=1, index2::Int=2; nbins::Any="auto")

Calculate 2-D pmf.

See also: [`pmf1`, `pmf2_gs`, `pmf2plus`](@ref)
"""
function pmf2(vec::Array{Float64},index1::Int=1,index2::Int=2;nbins::Any="auto")

    if nbins == "auto"
        nbins = size(vec,1) ^ 0.5 |> floor |> Int
    elseif typeof(nbins) == String
        error("Acceptable value of nbins is one of \"auto\" or Int (must be bigger than 1)")
    end
    minbin1 = minimum(vec[:,index1])
    maxbin1 = maximum(vec[:,index1]) 
    minbin2 = minimum(vec[:,index2])
    maxbin2 = maximum(vec[:,index2])
    bin1 = (maxbin1 - minbin1) / nbins
    bin2 = (maxbin2 - minbin2) / nbins
    bins1 = zeros(nbins,1)[:,1]
    bins2 = zeros(nbins,1)[:,1]
    bins1 = bins1 .- bins1
    bins2 = bins2 .- bins2
    for i = 1:nbins
        bins1[i] = bins1[i] + (i-1)*bin1
        bins2[i] = bins2[i] + (i-1)*bin2
    end
    bins1 = bins1 .+ minbin1
    bins2 = bins2 .+ minbin2
    bins1[nbins] = maxbin1
    bins2[nbins] = maxbin2
    pmf2 = fill(1.0,nbins,nbins)
    for i = 1:size(vec,1)
        tmp1 = searchsortedfirst(bins1, vec[i,index1])
        tmp2 = searchsortedfirst(bins2, vec[i,index2])
        pmf2[tmp1,tmp2] = pmf2[tmp1,tmp2] + 1
    end
    pmf2 = pmf2 ./ size(vec,1)
    pmf2 = broadcast(log,pmf2) .* -1
    pmf2 = pmf2 .- minimum(pmf2)

return pmf2', nbins
end

"""
    pmf2(vec::Array{Float64}, index1::Int=1, index2::Int=2; nbins::Any="auto")

Calculate 2-D pmf, also return free energy of each point.
"""
function pmf2plus(vec::Array{Float64},index1::Int=1,index2::Int=2;nbins::Any="auto")

    if nbins == "auto"
        nbins = size(vec,1) ^ 0.5 |> floor |> Int
    elseif typeof(nbins) == String
        error("Acceptable value of nbins is one of \"auto\" or Int (must be bigger than 1)")
    end
    minbin1 = minimum(vec[:,index1])
    maxbin1 = maximum(vec[:,index1]) 
    minbin2 = minimum(vec[:,index2])
    maxbin2 = maximum(vec[:,index2])
    bin1 = (maxbin1 - minbin1) / nbins
    bin2 = (maxbin2 - minbin2) / nbins
    bins1 = zeros(nbins,1)[:,1]
    bins2 = zeros(nbins,1)[:,1]
    bins1 = bins1 .- bins1
    bins2 = bins2 .- bins2
    for i = 1:nbins
        bins1[i] = bins1[i] + (i-1)*bin1
        bins2[i] = bins2[i] + (i-1)*bin2
    end
    bins1 = bins1 .+ minbin1
    bins2 = bins2 .+ minbin2
    bins1[nbins] = maxbin1
    bins2[nbins] = maxbin2
    pmf2 = fill(1.0,nbins,nbins)
    for i = 1:size(vec,1)
        tmp1 = searchsortedfirst(bins1, vec[i,index1])
        tmp2 = searchsortedfirst(bins2, vec[i,index2])
        pmf2[tmp1,tmp2] = pmf2[tmp1,tmp2] + 1
    end
    pmf2 = pmf2 ./ size(vec,1)
    pmf2 = broadcast(log,pmf2) .* -1
    pmf2 = pmf2 .- minimum(pmf2)
    # Free energy of each point.
    pmf2_ = fill(0.0,size(vec,1),1)
    for i = 1:size(vec,1)
        tmp1 = searchsortedfirst(bins1, vec[i,index1])
        tmp2 = searchsortedfirst(bins2, vec[i,index2])
        pmf2_[i] = pmf2[tmp1,tmp2]
    end

return pmf2', pmf2_, nbins
end

"""
    pmf2(vec::Array{Float64}, index1::Int=1, index2::Int=2; nbins::Any="auto")

Calculate 2-D pmf with Gauss smoothing.

See also: [`pmf1`, `pmf2`, `pmf2plus`](@ref)
"""
function pmf2_gs(vec::Array{Float64},index1::Int=1,index2::Int=2;
                nbins::Any="auto",cutoff::Int=5,sigma::Number=1)

    if nbins == "auto"
        nbins = size(vec,1) ^ 0.5 |> floor |> Int
    elseif typeof(nbins) == String
        error("Acceptable value of nbins is one of \"auto\" or Int (must be bigger than 1)")
    end
    pi = MathConstants.pi
    e = MathConstants.e
    minbin1 = minimum(vec[:,index1])
    maxbin1 = maximum(vec[:,index1]) 
    minbin2 = minimum(vec[:,index2])
    maxbin2 = maximum(vec[:,index2])
    bin1 = (maxbin1 - minbin1) / nbins
    bin2 = (maxbin2 - minbin2) / nbins
    bins1 = zeros(nbins,1)[:,1]
    bins2 = zeros(nbins,1)[:,1]
    bins1 = bins1 .- bins1
    bins2 = bins2 .- bins2
    for i = 1:nbins
        bins1[i] = bins1[i] + (i-1)*bin1
        bins2[i] = bins2[i] + (i-1)*bin2
    end
    bins1 = bins1 .+ minbin1
    bins2 = bins2 .+ minbin2
    bins1[nbins] = maxbin1
    bins2[nbins] = maxbin2
    pmf2 = fill(1.0,nbins,nbins)
    # Perform PMF calculation with Gaussian smooth.
    for i = 1:size(vec,1)
        tmp1 = searchsortedfirst(bins1, vec[i,index1])
        tmp2 = searchsortedfirst(bins2, vec[i,index2])
        pmf2[tmp1,tmp2] = pmf2[tmp1,tmp2] + 1/(pi*2)^0.5/sigma
        for j = 1:cutoff
            if (tmp1-j) > 1
                pmf2[(tmp1-j),tmp2] += 1/(pi*2)^0.5*e^(-j^2/(2*sigma^2))/sigma
            end
            if (tmp2-j) > 1
                pmf2[tmp1,(tmp2-j)] += 1/(pi*2)^0.5*e^(-j^2/(2*sigma^2))/sigma
            end
            if (tmp1+j) < nbins
                pmf2[(tmp1+j),tmp2] += 1/(pi*2)^0.5*e^(-j^2/(2*sigma^2))/sigma
            end
            if (tmp2+j) < nbins
                pmf2[tmp1,(tmp2+j)] += 1/(pi*2)^0.5*e^(-j^2/(2*sigma^2))/sigma
            end
            for k = 1:cutoff
                if (tmp1-j) > 1 && (tmp2-k) > 1
                    pmf2[(tmp1-j),(tmp2-k)] += 1/(pi*2)^0.5*e^(-(j^2+k^2)/(2*sigma^2))/sigma
                end
                if (tmp1-j) > 1 && (tmp2+k) < nbins
                    pmf2[(tmp1-j),(tmp2+k)] += 1/(pi*2)^0.5*e^(-(j^2+k^2)/(2*sigma^2))/sigma
                end
                if (tmp1+j) < nbins && (tmp2-k) > 1
                    pmf2[(tmp1+j),(tmp2-k)] += 1/(pi*2)^0.5*e^(-(j^2+k^2)/(2*sigma^2))/sigma
                end
                if (tmp1+j) < nbins && (tmp2+k) < nbins
                    pmf2[(tmp1+j),(tmp2+k)] += 1/(pi*2)^0.5*e^(-(j^2+k^2)/(2*sigma^2))/sigma
                end
            end
        end
    end
    pmf2 = pmf2 ./ size(vec,1)
    pmf2 = broadcast(log,pmf2) .* -1
    pmf2 = pmf2 .- minimum(pmf2)

return pmf2', nbins
end

"""
    pick_pmf(vec::Array{Float64}, index::Int=1; nbins::Any="auto")

Pick the point with longest distance to the point with lowest PMF (1-D).

See also: [`pick_pmf2`, `pick_pmf2plus`](@ref)
"""
function pick_pmf(vec::Array{Float64},index::Int=1;nbins::Any="auto")

    # Since we assume equilibrium sampling, the pmf of rare conformation must 
    # be higher than the main conformation. Thus, we just find the point (with
    # maximum distance to the point with lowest pmf) with as the input of next step.
    # First we should find the position of lowest pmf.
    fe, n = pmf1(vec,1;nbins = nbins)
    order = sortperm(fe[:,1])
    minbin = minimum(vec[:,index])
    maxbin = maximum(vec[:,index])
    pos = minbin + (maxbin - minbin) / n * (order[1] - 1)
    # calculate distance to pos1 of each point.
    dis = broadcast(abs,(vec[:,index] .- pos))
    disorder = sortperm(dis[:,1],rev=true)
    # pick first and last 10 points.
    point = disorder[1:10] .- 1
    basin = disorder[(size(disorder,1)-9):size(disorder,1)] .- 1

return basin, point
end

"""
    pick_pmf2(vec::Array{Float64}, index1::Int=1, index2::Int=2,
    restrict::Union{Array{Int64,1},Nothing}=nothing; nbins::Any="auto")

Pick the point with longest distance to the point with lowest PMF (2-D).

# Arguments:
- `vec::Array{Float64}`: input data array.
- `index1::Int=1`: x-axis data series (column in vec).
- `index2::Int=2`: y-axis data series (column in vec).
- `restrict::Union{Array{Int64,1},Nothing}=nothing`: points that can be picked.
- `nbins::Any="auto"`: bin number used to calculate PMF.

See also: [`pick_pmf`, `pick_pmf2plus`](@ref)
"""
function pick_pmf2(vec::Array{Float64},index1::Int=1,index2::Int=2,
         restrict::Union{Array{Int64,1},Nothing}=nothing;nbins::Any="auto")

    # Not the same as pick_pmf, this function use pmf from first two coordinates.
    fe, n = pmf2(vec,index1,index2;nbins=nbins)
    fe = reshape(fe,(n^2,1))[:,1]
    order = sortperm(fe[:,1])
    order1 = (order[1]/n) |> floor |> Int
    order2 = rem(order[1],n)
    minbin1 = minimum(vec[:,index1])
    maxbin1 = maximum(vec[:,index1])
    minbin2 = minimum(vec[:,index2])
    maxbin2 = maximum(vec[:,index2])
    pos1 = minbin1 + (maxbin1 - minbin1) / n * (order1-1)
    pos2 = minbin2 + (maxbin2 - minbin2) / n * (order2-1)
    # calculate distance to (pos1 pos2) of each point.
    vec[:,index1] = vec[:,index1] .- pos1
    vec[:,index2] = vec[:,index2] .- pos2
    dis = (vec[:,index1] .^ 2 .+ vec[:,index2] .^ 2) .^ 0.5
    # pick first and last 10 points.
    if restrict == nothing
        disorder = sortperm(dis[:,1])
        basin = disorder[1:10] .- 1
        disorder = sortperm(dis[:,1],rev=true)
        point = disorder[1:10] .- 1
    else
        if size(restrict,1) >= 10
            disorder = sortperm(dis[:,1])
            basin = disorder[1:10] .- 1
            disorder = sortperm(dis[restrict,1],rev=true)
            point = disorder[1:10] .- 1
        else
            disorder = sortperm(dis[:,1])
            basin = disorder .- 1
            disorder = sortperm(dis[restrict,1],rev=true)
            point = disorder .- 1
        end
    end

return basin, point
end

"""
    pick_pmf2plus(vec::Array{Float64}, index1::Int=1, index2::Int=2,
    restrict::Union{Array{Int64,1},Nothing}=nothing; nbins::Any="auto")

Pick the point with longest distance to the point with lowest PMF (2-D).
Only points with PMF larger than the median PMF of all points will be picked.

# Arguments:
- `vec::Array{Float64}`: input data array.
- `index1::Int=1`: x-axis data series (column in vec).
- `index2::Int=2`: y-axis data series (column in vec).
- `restrict::Union{Array{Int64,1},Nothing}=nothing`: points that can be picked.
- `nbins::Any="auto"`: bin number used to calculate PMF.

See also: [`pick_pmf`, `pick_pmf2`](@ref)
"""
function pick_pmf2plus(vec::Array{Float64},index1::Int=1,index2::Int=2,
         restrict::Union{Array{Int64,1},Nothing}=nothing;nbins::Any="auto")

    # Not the same as pick_pmf, this function use pmf from first two coordinates.
    fe, n = pmf2(vec,index1,index2;nbins=nbins)
    fe = reshape(fe,(n^2,1))[:,1]
    order = sortperm(fe[:,1])
    order1 = (order[1]/n) |> floor |> Int
    order2 = rem(order[1],n)
    minbin1 = minimum(vec[:,index1])
    maxbin1 = maximum(vec[:,index1])
    minbin2 = minimum(vec[:,index2])
    maxbin2 = maximum(vec[:,index2])
    pos1 = minbin1 + (maxbin1 - minbin1) / n * (order1-1)
    pos2 = minbin2 + (maxbin2 - minbin2) / n * (order2-1)
    # calculate distance to (pos1 pos2) of each point.
    vec[:,index1] = vec[:,index1] .- pos1
    vec[:,index2] = vec[:,index2] .- pos2
    dis = (vec[:,index1] .^ 2 .+ vec[:,index2] .^ 2) .^ 0.5
    # We should also know energy of each point, and only pick 
    # points with energy higher than mean free energy.
    fe, fe_, n = pmf2plus(vec,index1,index2;nbins=nbins)
    fe_mean = middle(fe_)
    fe_sign = max.(sign.(fe_ .- fe_mean),0)
    # pick first and last 10 points.
    if restrict == nothing
        disorder = sortperm(dis[:,1])
        basin = disorder[1:10] .- 1
        # dorp low free energy points
        dis = dis .* fe_sign 
        disorder = sortperm(dis[:,1],rev=true)
        point = disorder[1:10] .- 1
    else
        if size(restrict,1) >= 10
            disorder = sortperm(dis[:,1])
            basin = disorder[1:10] .- 1
            # dorp low free energy points
            dis = dis .* fe_sign 
            disorder = sortperm(dis[restrict,1],rev=true)
            point = disorder[1:10] .- 1
        else
            disorder = sortperm(dis[:,1])
            basin = disorder .- 1
            # dorp low free energy points
            dis = dis .* fe_sign 
            disorder = sortperm(dis[restrict,1],rev=true)
            point = disorder .- 1
        end
    end

return basin, point
end

"""
    pick_maxpmf(vec::Array{Float64}, index::Int=1; nbins::Any="auto")

Pick points with largest PMF (1-D).

See also: [`pick_maxpmf2`](@ref)
"""
function pick_maxpmf(vec::Array{Float64},index::Int=1;nbins::Any="auto")

    # Find the points with maximum pmf.
    fe, n = pmf1(vec,1;nbins = nbins)
    order = sortperm(fe[:,1],rev=true)
    minbin = minimum(vec[:,index])
    maxbin = maximum(vec[:,index])
    pos = minbin + (maxbin - minbin) / n * (order[1] - 1)
    # calculate distance to pos1 of each point.
    dis = broadcast(abs,(vec[:,index] .- pos))
    disorder = sortperm(dis[:,1])
    # pick first and last 10 points.
    point = disorder[1:10] .- 1

return point
end

"""
    pick_maxpmf(vec::Array{Float64}, index1::Int=1, index2::Int=2; 
    nbins::Any="auto")

Pick points with largest PMF (2-D).

See also: [`pick_maxpmf`](@ref)
"""
function pick_maxpmf2(vec::Array{Float64},index1::Int=1,index2::Int=2;nbins::Any="auto")

    # Not the same as pick_pmf, this function use pmf from first two coordinates.
    fe, n = pmf2(vec,index1,index2;nbins=nbins)
    fe = reshape(fe,(n^2,1))[:,1]
    fe = fe .- maximum(fe)
    order = sortperm(fe[:,1],rev=true)
    order1 = 0; order2 = 0
    for i = 1:size(fe,1)
        if fe[order[i]] != 0
            order1 = (order[i]/n) |> floor |> Int
            order2 = rem(order[i],n)
            break
        end
    end
    minbin1 = minimum(vec[:,index1])
    maxbin1 = maximum(vec[:,index1])
    minbin2 = minimum(vec[:,index2])
    maxbin2 = maximum(vec[:,index2])
    pos1 = minbin1 + (maxbin1 - minbin1) / n * order1
    pos2 = minbin2 + (maxbin2 - minbin2) / n * order2
    # calculate distance to (pos1 pos2) of each point.
    vec[:,index1] = vec[:,index1] .- pos1
    vec[:,index2] = vec[:,index2] .- pos2
    dis = (vec[:,index1] .^ 2 .+ vec[:,index2] .^ 2) .^ 0.5
    disorder = sortperm(dis[:,1])
    # pick first and last 10 points.
    point = disorder[1:10] .- 1

return point
end

function pick_maxrmsd(traj::PyObject, reftraj::PyObject, index::Array{Int64,1},
         alignref::String="all",regselect::String="all",strucselect::String="all";
         mp::Bool=true)

    # pick one frame from traj, which contains max RMSD with the given 
    # reference conformation.
    refframe = index[1]
    traj = md.join([get(reftraj,refframe,nothing), traj])
    RMSD = get_rmsd_part_qcp(traj,alignref,regselect,strucselect)
    disorder = sortperm(RMSD[:,1],rev=true)
    point = disorder[1:10] .- 1

return point
end
