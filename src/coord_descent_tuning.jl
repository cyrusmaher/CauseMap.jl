########### start optimization functions
function CoordDescentOpt(source_series::Vector{Float64}, target_series::Vector{Float64}, 
                                            shadowmat_dict::Dict, distmat_dict::Dict,
                                            E_vals::AbstractVector{Int}, tau_s_vals::AbstractVector{Int}, tau_p_vals::AbstractVector{Int}; 
                                            libsizemin::Int=0, libsizemax::Int=0,
                                            npred::Int=0, pred_start::Int=0,
                                            nreps = 5, nboots=0)
    
    bestres=Dict()
    count = 0
    for xx in 1:nreps
        res12, evalcount = _CoordDescentOpt(source_series, target_series, 
                                                                    shadowmat_dict, distmat_dict,
                                                                    E_vals, tau_s_vals, tau_p_vals; 
                                                                    libsizemin=libsizemin, libsizemax=libsizemax,
                                                                    npred=npred, pred_start=pred_start,
                                                                    nboots=nboots)
        count += evalcount
        if xx == 1
            bestres = res12
        elseif res12["rho"] > bestres["rho"]
            bestres = res12
        end
    end
    println("Finished $(count) CCM evaluations")
    return bestres
end

    
function _CoordDescentOpt(source_series::Vector{Float64}, target_series::Vector{Float64}, 
                                                shadowmat_dict::Dict, distmat_dict::Dict,
                                                E_vals::AbstractVector{Int}, tau_s_vals::AbstractVector{Int}, tau_p_vals::AbstractVector{Int}; 
                                                libsizemin::Int=0, libsizemax::Int=0,
                                                npred::Int=0, pred_start::Int=0,
                                                nboots=0)

    E_vals = convert(Vector{Int}, E_vals)
    tau_s_vals = convert(Vector{Int}, tau_s_vals)
    tau_p_vals = convert(Vector{Int}, tau_p_vals)

    toopt = ["E"; "tau_s"; "tau_p"]

    current_vals = {"E"=>E_vals[1], "tau_s"=> tau_s_vals[1], "tau_p"=>tau_p_vals[1]}
    librange, res12 = calcCCM(source_series, target_series, 
                                                shadowmat_dict, distmat_dict,  
                                                current_vals["E"], current_vals["tau_s"], current_vals["tau_p"]; 
                                                libsizemin=libsizemin, libsizemax=libsizemax,
                                                npred=npred, pred_start=pred_start,
                                                quick=true, nboots=nboots)
    
    rhoinit = getrho(res12)
    best_vals = merge(current_vals, ["rho"=>rhoinit])
    
    all_vals = ["E"=>E_vals, "tau_s"=>tau_s_vals, "tau_p"=>tau_p_vals]

    iternum = 1
    count = 0
    while true
        var_count = 0
        println("========Starting iteration $iternum===========")
        for var in shuffle(toopt)
            nochange_bool, evalcount = optvar(source_series, target_series, shadowmat_dict, distmat_dict, 
                                                                    all_vals, current_vals, best_vals, var, libsizemax, npred, pred_start; nboots=nboots)
            count += evalcount
            current_vals[var] = best_vals[var][1] # make sure optimization of next variable is done with best value of this one

            if nochange_bool
                var_count += 1
            end
        end

        if var_count == size(toopt, 1)  # check for convergence (all variables unchanged)
            break
        end
        iternum += 1
    end
    return best_vals, count
end


function optvar(source_series::AbstractVector, target_series::AbstractVector, shadowmat_dict::Dict, 
    distmat_dict::Dict, all_vals::Dict, current_vals::Dict, 
    best_vals::Dict, var::ASCIIString, libsizemax::Int, npred::Int, pred_start::Int; nlag::Int=10, b_offset=1, nboots=0)
    
    
    if length(all_vals[var]) < 2
        return true, 0
    end
    count::Int= 0
    looplist     = shuffle(setdiff(all_vals[var], current_vals[var]))
    val_count::Int    = 0
    update_count::Int = 0
    for val in looplist
        current_vals[var] = val # update variable of interest
        libsizemax = min(libsizemax, calclibsizemax(length(source_series), current_vals["E"], current_vals["tau_s"], current_vals["tau_p"]))
        libsizemin = max(current_vals["E"] + b_offset + 1, libsizemax-nlag)
        librange, res12 = calcCCM(source_series, target_series, 
                                                    shadowmat_dict, distmat_dict,
                                                    current_vals["E"], 
                                                    current_vals["tau_s"], current_vals["tau_p"];
                                                    libsizemin=libsizemin, libsizemax=libsizemax, 
                                                    npred=npred, pred_start=pred_start,
                                                    quick=true, nboots=nboots)
        rho = getrho(res12)
        if rho > best_vals["rho"] # update best value if you have an improvement
            best_vals[var] = val
            best_vals["rho"] = rho
            update_count += 1
        else
            val_count += 1
        end
        count += 1 
    end

    if val_count + update_count !=length(looplist)
        error("Missing iteration in _optvar!!") 
    end

    if val_count == size(all_vals[var],1)-1
        return true, count
    else
        return false, count
    end
end


function getrho(rhos::AbstractArray)
    rhos = rhos[!isnan(rhos)]   
    if length(rhos) > 0
        rho = median(rhos)
    else
        rho = NaN
    end
    return rho
end
### end functions for CoordDescentOpt optimization
