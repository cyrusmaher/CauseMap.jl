module CauseMap
using Base.LinAlg.BLAS 
using PyCall

include("plotting.jl")
include("coord_descent_tuning.jl")
include("calc_manifolds.jl")

export makesingleplot, makeoptimizationplots, optandcalcCCM, precalc_manif_dists, calcCCM, load_example_data

function load_example_data(dsname)
    if dsname == "ParaDidi"
        ds   = readdlm(Pkg.dir("CauseMap") * "/examples/vr.raw_fixed.txt", '\t', Float64, header=false)  # read in data
        para = ds[:, 2]
        didi = ds[:, 3]
        return para, didi
    else
        throw(ArgumentError("Did not recognize $dsname dataset"))
    end
end


function processpredvals_simple(predvals::Array{Float64,2}, targvals::Vector{Float64}, lib_size::Int)
    nanval = NaN

    if all(isnan(predvals))
        println("All predvals were nanval in process predvals of lib size of $lib_size")
        return nanval
    end
    rhos = fill(nanval, size(predvals,2))
    
    for xx in 1:size(predvals,2)
        touse = !isnan(predvals[:,xx]) & !isnan(targvals)
        
        if sum(touse) < 10
            rhos[xx] = nanval
        else
            rhos[xx] = cor(predvals[touse,xx], targvals[touse])
            if isnan(rhos[xx])
                warn("in process predvals, why is rho nan?")
            end
        end
    end
    toret = !isnan(rhos)
    if sum(toret) > 0
        return median(rhos[toret])
    else
        println("no non-nan values")
        return nanval
    end

end


function weightfunc(distances::Array{Float64}; kernelargs...)
    w = [exp(-d) for d in distances]
    return w / sum(w)
end


function getdist!(dist_top::Vector{Float64}, inds_touse::AbstractVector{Int}, num_neighbors::Int) 
    if dist_top[1] == 0
        for xx in 1:num_neighbors
            if xx == 1
                dist_top[xx] = 1
            else
                dist_top[xx] = 1e200
            end
        end
        return 1e-100
    else
        min_dist::Float64 = dist_top[1]
        for xx in num_neighbors:-1:1
            dist_top[xx] = dist_top[xx] / min_dist
        end
        return min_dist
    end
end


function getpredstartstop(nobs::Int, ll::Int, lib_size::Int, npred::Int, pred_start_min::Int)
    left::Int = iceil(npred / 2)
    right::Int = npred - left - 1
    midpoint::Int = iceil(ll + lib_size / 2)
    lp::Int = midpoint - left
    rp::Int = midpoint + right

    if lp < pred_start_min
        rp += (pred_start_min - lp)
        lp  = pred_start_min
    end

    if rp > nobs
        lp -= rp % nobs
        rp  = nobs
    end
    
    if (rp - lp + 1) != npred
        error("start and stop not equal to npred")
    end

    return lp, rp
end


### this function accounts for ~70% of algorithm run time
function calcdistslice!(source_dists::Array{Float64, 2}, dist_top::Vector{Float64}, slice_inds::AbstractVector{Int}, topred::Int, nn::Range1{Int})    
    if in(topred, slice_inds)
        slice_inds = convert(Vector{Int}, slice_inds)
        splice!(slice_inds, findfirst([yy == topred for yy in slice_inds])) # remove topred, but don't move any of the preceding elements in the array
    end

    # this line accounts for ~60% of algorithm run time
    inds_touse::Vector{Int} = slice_inds[sortperm(source_dists[slice_inds, topred])[nn]] 
    for xx in nn
        dist_top[xx] = source_dists[inds_touse[xx], topred]  # this is not returned because it is written in-place
    end
    return inds_touse   
end


function prepgetpred(source_manifold, lib_start, lib_size, tau_p, npred, num_neighbors; nboots=0)
    nobs::Int = size(source_manifold, 2)
    lib_stop::Int = nobs - lib_size - tau_p
    nlib::Int = lib_stop - lib_start + 1
    if nlib < 1
        println("No libs of this size in dataset. Lib_start is $lib_start, lib_stop is $lib_stop")
        return fill(NaN, 1, 1), fill(NaN, 1,1)
    end

        if nboots > 0
        min_distances = fill(NaN, npred, nboots)
        predvals = fill(NaN, npred, nboots)
        targvals = fill(NaN, npred)
    else
        min_distances = fill(NaN, npred, nlib)
        predvals = fill(NaN, npred, nlib)
        targvals = fill(NaN, npred)
    end
    lib_end::Int = lib_start + lib_size - 1
    dist_top = fill(NaN, num_neighbors)
    return nobs, lib_stop, nlib, predvals, targvals, min_distances, lib_end, dist_top
end


function getpredvals_boot(source_manifold::Array{Float64,2}, source_dists::Array{Float64,2}, 
                                             target_series::Vector{Float64}, lib_size::Int,
                                             lib_start::Int, num_neighbors::Int, 
                                             tau_p::Int, npred::Int, pred_start_min::Int, nboots::Int)
    # Use bootstrap samples to generate libraries
    
    nobs, lib_stop, nlib, predvals, targvals, min_distances, lib_end, dist_top = prepgetpred(source_manifold,
                                                                                                                                lib_start, lib_size, 
                                                                                                                                tau_p, npred, num_neighbors; 
                                                                                                                                nboots=nboots)
    # lib_stop is the last possible sliding window START index (stop of the start...confusing...)
    # lib_end is the last in the current window
    nn = 1:num_neighbors
    for (pp, topred) in enumerate(
                                                    pred_start_min:min(pred_start_min + npred, length(target_series))
                                                    )  
        targvals[pp] = target_series[topred]
        for xx in 1: nboots  
            slice_inds::Array{Int, 1} = [rand(lib_start:lib_stop) for xx in 1:lib_size]
            
            inds_touse = calcdistslice!(source_dists, dist_top, slice_inds, topred, nn)   
            min_distances[pp, xx] = getdist!(dist_top, inds_touse, num_neighbors)
            weights = weightfunc(dist_top)
            predvals[pp, xx] = dot(weights, target_series[inds_touse + tau_p]) 
        end   
    end

    rhos = processpredvals_simple(predvals, targvals, lib_size)
    return predvals, min_distances, rhos
end


function getpredvals_sw(source_manifold::Array{Float64,2}, source_dists::Array{Float64,2},
    target_series::Vector{Float64}, lib_size::Int,
    lib_start::Int, num_neighbors::Int, 
    tau_p::Int, npred::Int, pred_start_min::Int)
    # Use a sliding window for your library. This is the traditional form of CCM

    nobs, lib_stop, nlib, predvals, targvals, min_distances, lib_end, dist_top = prepgetpred(source_manifold,
                                                                                                                                lib_start, lib_size, 
                                                                                                                                tau_p, npred, num_neighbors)
    
    # lib_stop is the last possible sliding window START index (stop of the start..confusing...)
    # lib_end is the last in the current window
    nn = 1:num_neighbors
    
    start_count::Int = 1
    for ll in lib_start:lib_stop
        ll_end::Int     = ll + lib_size - 1
        slice_inds = ll:ll_end
        predstart, predstop = getpredstartstop(nobs, ll, lib_size, npred, pred_start_min)
        
        pred_count::Int=1
        for topred in predstart:predstop
            slice_inds = ll:ll_end
            inds_touse = calcdistslice!(source_dists, dist_top, slice_inds, topred, nn)
            
            min_distances[pred_count, start_count] = getdist!(dist_top, inds_touse, num_neighbors)
            
            weights = weightfunc(dist_top)
            predvals[pred_count, start_count] = dot(weights, target_series[inds_touse+tau_p])
            targvals[pred_count] = target_series[topred]
            
            pred_count += 1
        end
        start_count += 1
        lib_end += 1
    end
    rhos = processpredvals_simple(predvals, targvals, lib_size)
    return predvals, min_distances, rhos
end


function cross_mapping(source_manif_dict::Dict, source_dist_dict::Dict, 
                                        target_series::AbstractVector,
                                        nobs::Int,
                                        libsizemin::Int, libsizemax::Int, E::Int,
                                        tau_s::Int, tau_p::Int, npred::Int, 
                                        pred_start_min::Int, num_neighbors, args...; 
                                        lib_start::Int=0, nboots=0)
                                        

    res12 = fill(NaN, libsizemax - libsizemin + 1)
    count = 1

    nlibpluslibsize = nobs - tau_p - lib_start + 1

    for lib_size in libsizemin:libsizemax
        nlib = nlibpluslibsize - lib_size
        if nboots > 0
            predvals, min_distances, rhos = getpredvals_boot(source_manif_dict[tau_s][E], source_dist_dict[tau_s][E], 
                                target_series, lib_size, lib_start, num_neighbors, tau_p, npred, pred_start_min, nboots)
        else
            predvals, min_distances, rhos = getpredvals_sw(source_manif_dict[tau_s][E], source_dist_dict[tau_s][E], 
                                target_series, lib_size, lib_start, num_neighbors, tau_p, npred, pred_start_min)
        end
        res12[count] = rhos

        count += 1
    end
    if count == 1
        warn("Why was there no loop?")
    end
    return res12
end


function calclibstart(shadowmat_dict::Dict, E::Int, tau_s::Int)

    nanrows::Array{Bool, 1} = [in(NaN, shadowmat_dict[tau_s][E][:,xx]) for xx in 1:size(shadowmat_dict[tau_s][E], 2)]
    if in(true, nanrows)
        return (maximum(find(nanrows)) + 1)
    else
        return 1
    end
end

function calclibsizemax(nobs::Int, E::Int, tau_s::Int, tau_p::Int)  ## REMOVE
    return ifloor(nobs/tau_s) - E - tau_p
end

function check_lib_pred_defaults(libsizemin, libsizemax, lib_start, 
                                                        npred, pred_start,
                                                        nobs, E, b_offset, nlag, 
                                                        tau_p, tau_s, shadowmat_dict; 
                                                        quick=false, lmin_def=8, lmax_offset=2)
    if libsizemin == 0
        libsizemin = lmin_def
    end
    
    
    if libsizemax == 0
        libsizemax = calclibsizemax(nobs, E, tau_s, tau_p)
    end
    libsizemax = min(libsizemax, nobs - tau_p - lib_start)


    if quick 
        libsizemin = max(E + b_offset + 1, libsizemax - nlag)
    else
        libsizemin = max(libsizemin, E + b_offset + 1)
    end

    if lib_start == 0 
        lib_start = calclibstart(shadowmat_dict, E, tau_s)
    end 

    if npred == 0
        npred = nobs
    end
    
    if pred_start == 0
        pred_start = 1
    end
    pred_start = max(lib_start, pred_start)
    

    if pred_start + npred - 1 > nobs
        npred -= (pred_start + npred - 1) - nobs
    end

    return libsizemin, libsizemax, lib_start, npred, pred_start
end


function calcCCM(var1::AbstractVector, var2::AbstractVector,
                                shadowmat_dict::Dict, distmat_dict::Dict,
                                E::Int, tau_s::Int, tau_p::Int; 
                                libsizemin::Int=0, libsizemax::Int=0,
                                npred::Int=0, pred_start::Int=0,
                                lib_start::Int=0, b_offset::Int=1, 
                                quick=false, nlag=10, nboots=0)

    nobs = length(var1)
    num_neighbors = E + b_offset

    libsizemin, libsizemax, lib_start, npred, pred_start = check_lib_pred_defaults(libsizemin, libsizemax, lib_start, 
                                                                                                                            npred, pred_start,
                                                                                                                            nobs, E, b_offset, nlag, 
                                                                                                                            tau_p, tau_s, shadowmat_dict; quick=quick)

    ################## clean input parameters
    
    #########
    ##### Raise error or warning messages as needed.
    if libsizemin > libsizemax
        println("Libsizemin, libsizemax: ($libsizemin, $libsizemax)")
        println("E, tau_s, tau_p, nobs, lib_start : ($E, $tau_s, $tau_p, $nobs, $lib_start)")
        warn("why is libsizemin less than libsizemax? returning NaN")
        return 0:0, fill(NaN, 1)
    end

    ###################
    if nboots > 0
        warn("Bootstrapped libraries are still experimental")
    end
    res12 = cross_mapping(shadowmat_dict, distmat_dict, var2, nobs, libsizemin, libsizemax, E,
                            tau_s, tau_p, npred, pred_start, num_neighbors; lib_start = lib_start, nboots=nboots)

    return libsizemin:libsizemax, res12
end


function optandcalcCCM(vec1::AbstractVector, vec2::AbstractVector, 
                                        E_vals::AbstractVector,  tau_s_vals::AbstractVector, tau_p_vals::AbstractVector; 
                                        nreps=5, b_offset=1, nboots=0, 
                                        npred::Int=0, pred_start::Int=0,
                                        libsizemin::Int=0, libsizemax::Int=0)
    """
    vec1: Time series 1
    vec2: Time series 2
    libsizemin: Minimum library size
    libsizemax: Maximum library size
    E_vals: A vector of dimensions to try
    tau_s_vals: A vector of lag lengths to use for manifold reconstruction
    tau_p_vals: A vector of lag lengths to try for the causal effect 
    npred: Number of points to predict 
    pred_start: First point to predict
    ## kwargs
    nreps: Number of coordinate descent runs
    """

    shadowmat_dict_vec1, distmat_dict_vec1 = precalc_manif_dists(E_vals, tau_s_vals, vec1)
    shadowmat_dict_vec2, distmat_dict_vec2 = precalc_manif_dists(E_vals, tau_s_vals, vec2)

    println("CoordDescentOpt opt1")
    res12 = CoordDescentOpt(vec1, vec2, 
                                                shadowmat_dict_vec1, distmat_dict_vec1,  
                                                E_vals, tau_s_vals, tau_p_vals;
                                                libsizemin=libsizemin, libsizemax=libsizemax, 
                                                npred=npred, pred_start=pred_start,
                                                nreps=nreps, nboots=nboots)
    
    println("CoordDescentOpt opt2")
    res21 = CoordDescentOpt(vec2, vec1, 
                                                shadowmat_dict_vec2, distmat_dict_vec2,  
                                                E_vals, tau_s_vals, tau_p_vals;
                                                libsizemin=libsizemin, libsizemax=libsizemax, 
                                                npred=npred, pred_start=pred_start,
                                                nreps=nreps, nboots=nboots)
    
    libsizemin_12 = max(res12["E"] + b_offset + 1, libsizemin)
    libsizemin_21 = max(res21["E"] + b_offset + 1, libsizemin)
    
    println("starting calcCCM1")
    librange12, yval_12 = calcCCM(vec1, vec2, 
                                                        shadowmat_dict_vec1, distmat_dict_vec1,
                                                        res12["E"], res12["tau_s"], res12["tau_p"]; 
                                                        libsizemin=libsizemin_12, libsizemax=libsizemax,
                                                        npred=npred, pred_start=pred_start,
                                                        nboots=nboots)
    println("starting calcCCM2")
    librange21, yval_21 = calcCCM(vec2, vec1, 
                                                        shadowmat_dict_vec2, distmat_dict_vec2,  
                                                        res21["E"], res21["tau_s"], res21["tau_p"], ; 
                                                        libsizemin=libsizemin_21, libsizemax=libsizemax,
                                                        npred=npred, pred_start=npred,
                                                        nboots=nboots)

    return (librange12, yval_12), (librange21, yval_21)
end
########## end CCM functions
end # module
