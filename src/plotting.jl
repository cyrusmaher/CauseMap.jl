function makesingleplot(vec1::AbstractVector, vec2::AbstractVector, 
                                        E::Int, tau_s::Int, tau_p::Int,  
                                        var1name::ASCIIString, var2name::ASCIIString; 
                                        lib_start=1, xmin=false, xmax=false, ymin=false, 
                                        ymax=false, nboots=0, plot=true, 
                                        libsizemin=0, libsizemax=0, 
                                        npred=0, pred_start=0)
    """
    vec1: Time series 1
    vec2: Time series 2
    libsizemin: Minimum library size
    ibsizemax: Maximum library size
    E: System dimensionality
    tau_s: Lag for manifold reconstruction 
    tau_p: Lag for causal effect
    npred: Number of points to predict
    pred_start: Start for prediction 
    var1name: Name of first variable
    var2name: Name of second variable
    lib_start: Start of library
    """
    @pyimport matplotlib.pyplot as plt
    shadowmat_dict_vec1, distmat_dict_vec1 = precalc_manif_dists(E:E, tau_s:tau_s, vec1)
    shadowmat_dict_vec2, distmat_dict_vec2 = precalc_manif_dists(E:E, tau_s:tau_s, vec2)

    librange12, yval_12 = calcCCM(vec1, vec2, 
                                                        shadowmat_dict_vec1, distmat_dict_vec1, 
                                                        E, tau_s, tau_p; libsizemin=libsizemin, libsizemax=libsizemax, 
                                                        npred=npred, pred_start=pred_start,
                                                        lib_start=lib_start, nboots=nboots
                                                    )

    librange21, yval_21 = calcCCM(vec2, vec1, 
                                                        shadowmat_dict_vec2, distmat_dict_vec2, 
                                                        E, tau_s, tau_p; libsizemin=libsizemin, libsizemax=libsizemax, 
                                                        npred=npred, pred_start=pred_start,
                                                        lib_start=lib_start, nboots=nboots
                                                    )

    plt.plot(librange12, yval_12, label = "$var2name influences $var1name?")
    plt.plot(librange21, yval_21, label = "$var1name influences $var2name?")
    plt.xlabel("L")
    plt.ylabel("\$\\rho_{ccm}\$")
    plt.legend(loc="lower right")
    ax = plt.gca()
    ax[:xaxis][:set_ticks_position]("bottom")
    ax[:yaxis][:set_ticks_position]("left")
    if xmin != false
        ax[:set_xlim]((xmin, xmax))
    end
    if ymin != false
        ax[:set_ylim]((ymin, ymax))
    end
    if plot
        plt.show()
    else
        plt.close()
    end
end


function makeoptimizationplots(vec1::AbstractVector, vec2::AbstractVector, E_vals::AbstractVector,  
                                                    tau_s_vals::AbstractVector, tau_p_vals::AbstractVector,
                                                    var1name::ASCIIString, var2name::ASCIIString; 
                                                    nreps=5, b_offset=1, ncols=28, left_E=false, left_tau_p=false, 
                                                    right_E=false, right_tau_p=false, lagunit=1, unit=false, 
                                                    imfont="medium", nboots=0, legend=true, 
                                                    lib_start=1, plot=true, libsizemin=0, libsizemax=0, 
                                                    pred_start=0, npred=0, show_tau_s=true)
        
    println("Calculating manifolds")
    shadowmat_dict_vec1, distmat_dict_vec1 = precalc_manif_dists(E_vals, tau_s_vals, vec1)
    shadowmat_dict_vec2, distmat_dict_vec2 = precalc_manif_dists(E_vals, tau_s_vals, vec2)

    println("\nCoordDescent opt1")
    res12 = CoordDescentOpt(vec1, vec2, 
                                                shadowmat_dict_vec1, distmat_dict_vec1, 
                                                E_vals, tau_s_vals, tau_p_vals; 
                                                libsizemin=libsizemin, libsizemax=libsizemax,
                                                npred=npred, pred_start=pred_start, 
                                                nreps=nreps, nboots=nboots)
    println("\nCoordDescent opt2")
    res21 =CoordDescentOpt(vec2, vec1, 
                                                shadowmat_dict_vec2, distmat_dict_vec2, 
                                                E_vals, tau_s_vals, tau_p_vals; 
                                                libsizemin=libsizemin, libsizemax=libsizemax,
                                                npred=npred, pred_start=pred_start, 
                                                nreps=nreps, nboots=nboots)
    
    libsizemin_12 = max(res12["E"] + b_offset + 1, 10)
    libsizemin_21 = max(res21["E"] + b_offset + 1, 10)
    
    println("\nstarting calcCCM1")
    librange12, yval_12 = calcCCM(vec1, vec2, 
                                                        shadowmat_dict_vec1, distmat_dict_vec1, 
                                                        res12["E"], res12["tau_s"], res12["tau_p"];
                                                        libsizemin=libsizemin_12, libsizemax=libsizemax,
                                                        npred=npred, pred_start=pred_start, 
                                                        lib_start=1, nboots=nboots)
    println("starting calcCCM2")
    librange21, yval_21 = calcCCM(vec2, vec1, 
                                                        shadowmat_dict_vec2, distmat_dict_vec2, 
                                                        res12["E"], res12["tau_s"], res12["tau_p"];
                                                        libsizemin=libsizemin_21, libsizemax=libsizemax,
                                                        npred=npred, pred_start=pred_start, 
                                                        lib_start=1, nboots=nboots)
    
    @pyimport matplotlib.pyplot as plt
    fig = plt.figure()
    label1 = "\n$(ucfirst(var2name)) influences $var1name?"
    label2 = "\n$(ucfirst(var1name)) influences $var2name?"
    
    rho1   = @sprintf("%.2f", res12["rho"])
    rho2   = @sprintf("%.2f", res21["rho"])
    
    if show_tau_s
        statsuff1 = "$(res12["tau_s"])"
        statsuff2 = "$(res21["tau_s"])"
    else
        statsuff1 = ""
        statsuff2 = ""
    end
    stats1 = "   E = $(res12["E"]), \$\\tau_{p}\$ = $(int(res12["tau_p"] * lagunit))" * statsuff1
    stats2 = "   E = $(res21["E"]), \$\\tau_{p}\$ = $(int(res21["tau_p"] * lagunit))" * statsuff2

    ax1 = plt.subplot2grid((2,ncols), (0,0), rowspan=2, colspan=ifloor(ncols/3))
    
    ax1[:plot](librange12, yval_12, label="$label1\n$stats1")
    ax1[:plot](librange21, yval_21, label="$label2\n$stats2")
    
    ax1[:yaxis][:set_ticks_position]("left")
    ax1[:xaxis][:set_ticks_position]("bottom")

    ax1[:set_ylabel]("\$\\rho_{ccm}\$", fontsize="x-large")
    ax1[:set_xlabel]("L", fontsize="large", labelpad=12)
    
    if legend == true
        ax1[:legend](loc=4, fontsize="x-small")
    end
    
    mat12, Es12, taus12 = get_E_taupcurves(vec1, vec2, 
                                                                        shadowmat_dict_vec1, distmat_dict_vec1, 
                                                                        libsizemax, res12["E"], res12["tau_s"], res12["tau_p"], 
                                                                        E_vals, tau_p_vals, npred, pred_start; 
                                                                        left_E=left_E, left_tau_p=left_tau_p, 
                                                                        right_E=right_E, right_tau_p=right_tau_p, nboots=nboots
                                                                    )
    
    mat21, Es21, taus21 = get_E_taupcurves(vec2, vec1, 
                                                                        shadowmat_dict_vec2, distmat_dict_vec2,
                                                                        libsizemax, res21["E"], res21["tau_s"], res21["tau_p"], 
                                                                        E_vals, tau_p_vals, npred, pred_start;
                                                                        left_E=left_E, left_tau_p=left_tau_p, 
                                                                        right_E=right_E, right_tau_p=right_tau_p, nboots=nboots
                                                                    )

    catmat = vcat(mat12, mat21)
    if minimum(catmat) < 0
        vmin = 0
        llab = "\$\\leq 0\$"
    else
        vmin = minimum(catmat)
        llab = @sprintf("%.2f", vmin)
    end

    vmax = 1  # maximum correlation
    
    if unit != false
        xlabel = "\$\\tau_p\$\ ($unit)"
    else
        xlabel =  "\$\\tau_p\$"
    end

    ax2 = plt.subplot2grid((2,ncols), (0,ifloor(ncols/3)), colspan=2*ifloor(ncols/3))
    imax12 = plotheatmap(ax2, mat12, vmin, vmax, label1, taus12, Es12, imfont, lagunit, xlabel)

    ax3 = plt.subplot2grid((2,ncols), (1,ifloor(ncols/3)), colspan=2*ifloor(ncols/3))
    imax21 = plotheatmap(ax3, mat21, vmin, vmax, label2, taus21, Es21, imfont, lagunit, xlabel)

    ax4 = plt.subplot2grid((2,ncols), (0,ncols-1), rowspan=2)  # axis for colorbar
    
    cbar = fig[:colorbar](imax21, cax=ax4, ticks=[vmin, vmax], format = "%.2f")
    cbar[:ax][:set_yticklabels]([llab, "\$1\$"], fontsize="large")
    cbar[:set_label]("\$max\\ \\rho_{ccm}\$", fontsize="large", labelpad=2)
    cbar[:ax][:get_yaxis]()[:labelpad]=0
    
    if plot
        plt.show()
    else
        plt.close()
    end
end

function get_totest(val::Int, numtotest::Int, minval::Int, maxval::Int)
    if maxval - minval - 1 < numtotest
        println("Maxval and minval for E_taup curves are inconsistent with number to test")
        println("Fixing maxval")
        maxval = minval + numtotest - 1  
    end
    
    left::Int  = val-iceil(numtotest/2) 
    right::Int = val+ifloor(numtotest/2)-1
    
    if left < minval
        right += (minval-left)
        left   = minval
    end

    if right >  maxval
        left -= right%maxval
        right = maxval
    end
    
    if right-left+1 != numtotest
        error("start and stop not equal to npred")
    elseif left < minval
        println("left: $left, right: $right, Max: $maxval, min: $minval.")
        error("Left value less than minval")
    end
    
    return left, right 
end


function get_E_taupcurves(source_series::AbstractVector, target_series::AbstractVector, shadowmat_dict::Dict,
                                            distmat_dict::Dict, libsizemax::Int, E::Int, tau_s::Int, tau_p::Int, 
                                            E_vals::AbstractVector, tau_p_vals::AbstractVector, npred::Int, pred_start::Int;  
                                            nlag=10, b_offset=1, num_E=6, Estart=2, num_tau_p=14, tau_p_start=0,
                                            left_E=false, left_tau_p=false, right_E=false, right_tau_p=false, nboots=0)

    num_E = min(maximum(E_vals) - Estart + 1, num_E)
    num_tau_p = min(maximum(tau_p_vals-tau_p_start+1), num_tau_p)

    if (typeof(left_E) == Bool) && (left_E == false)
        left_E, right_E = get_totest(E, num_E, Estart, maximum(E_vals))
    end
    
    if (typeof(left_tau_p) == Bool) && (left_tau_p == false)
        left_tau_p, right_tau_p = get_totest(tau_p, num_tau_p, tau_p_start, maximum(tau_p_vals))
    end

    res = fill(NaN, num_E, num_tau_p)

    for (count_E, Eval) in enumerate(left_E:right_E)
        println("\nCalculating CCM max for E: $Eval")
        for (count_tau_p, tau_p_val) in enumerate(left_tau_p:right_tau_p)
            # libsizemin will be ignored. libsizemax will be reduced if it is too large
            println("\ttau_p: $tau_p_val")
            
            librange, rhos = calcCCM(source_series, target_series, 
                                                        shadowmat_dict, distmat_dict,
                                                        Eval, tau_s, tau_p_val;
                                                        libsizemin=libsizemax - nlag, libsizemax=libsizemax,
                                                        npred=npred, pred_start=pred_start, 
                                                        quick=true, nboots=nboots) 
            rhos = rhos[!isnan(rhos)]
            
            if length(rhos) > 0
                rho = median(rhos)
            else
                rho = NaN
            end

            res[count_E, count_tau_p] = rho
        end
    end
    return res, left_E:right_E, left_tau_p:right_tau_p
end

function plotheatmap(ax, mat, vmin, vmax, label, taus, Es, imfont, lagunit, xlabel)
    imax = ax[:imshow](mat, vmin=vmin, vmax=vmax, origin="lower")
    ax[:set_title](label, fontsize="medium")
    
    ax[:xaxis][:set_ticks_position]("bottom")
    ax[:yaxis][:set_ticks_position]("right")

    xticks = 0:(length(taus)-1)
    ax[:set_xticks](xticks[1:2:end]) # take every other tick
    ax[:set_yticks](0:(length(Es)-1))

    ax[:set_xticklabels](
                                    format_ticks(
                                                        (taus * lagunit)[1:2:end]
                                                        ), fontsize=imfont)
    ax[:set_yticklabels](Es, fontsize=imfont)

    ax[:set_xlim]((0, length(taus)-1))
    ax[:set_ylim]((0, length(Es)-1))
    ax[:set_xlabel](xlabel, fontsize=imfont)
    ax[:set_ylabel]("E", fontsize=imfont, labelpad = 12)
    return imax
end

function format_ticks(xticks)
    
    try
        xticks = convert(Vector{Float64}, xticks) # float range to float vector (otherwise subtraction fails)
    catch
        return xticks # This may fail if you have an int range
    end

    xticks_int = fill(0, length(xticks))
    try
        xticks_int = convert(Vector{Int}, xticks)  # If values are too far from ints, this raises an InexactError
    catch
        return xticks
    end

    variance = var(xticks)
    if all(abs(xticks - xticks_int) .< (variance * 1e-12))  # possibly more stringent than InexactError if values are small
        return xticks_int
    else
        return xticks
    end  
end
