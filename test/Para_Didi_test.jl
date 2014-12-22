using CauseMap
using PyCall


function getvectors(vec1::AbstractVector, vec2::AbstractVector, libsizemin::Int64, libsizemax::Int64, E::Int64, 
    tau_s::Int64, tau_p::Int64, npred::Int64, pred_start::Int64, lib_start::Int64=1)

    shadowmat_dict_vec1, distmat_dict_vec1 = precalc_manif_dists(3:3, 1:1, vec1)
    shadowmat_dict_vec2, distmat_dict_vec2 = precalc_manif_dists(3:3, 1:1, vec2)

    librange12, yval_12 = calcCCM(vec1, vec2, 
                                                        shadowmat_dict_vec1, distmat_dict_vec1, 
                                                        E, tau_s, tau_p; 
                                                        lib_start=lib_start, libsizemin=libsizemin, libsizemax=libsizemax,
                                                        npred=npred, pred_start=pred_start, pred_rand=false)
    librange21, yval_21 = calcCCM(vec2, vec1, 
                                                        shadowmat_dict_vec2, distmat_dict_vec2, 
                                                        E, tau_s, tau_p; 
                                                        lib_start=lib_start, libsizemin=libsizemin, libsizemax=libsizemax,
                                                        npred=npred, pred_start=pred_start, pred_rand=false)
    return yval_12, yval_21
end


function validatevecs(yval_12::AbstractVector, yval_21::AbstractVector; margin=.01,
                                            fname="Para_xmap_Didi_values.txt", col12=2, col21=3)
    col12 = 2
    col21 = 3
    
    ds_validation = readdlm(fname, '\t', Float64, header=true)[1]
    sum1 = (yval_12 - ds_validation[:, col12])
    sum1 = dot(sum1, sum1)
    sum2 = (yval_21 - ds_validation[:, col21])
    sum2 = dot(sum2, sum2)
    fail1 =  (sum1 > (length(yval_12)*margin))
    fail2 = (sum2 > (length(yval_21)*margin))
    if fail1 | fail2
        if fail1
            println("Para didi xmap off by $sum1")
        end
        if fail2
            println("Didi para xmap off by $sum2")
        end
        error("Error greater than a margin of $margin")
    end
end

E = 3
tau_s = 1
tau_p = 0

libsizemin = 8 
libsizemax = 58
pred_start = 13
npred = 60

ds   = readdlm("vr.raw_fixed.txt", '\t', Float64, header=false)
para = ds[:, 2]
didi = ds[:, 3]

para_didi, didi_para = getvectors(para, didi, libsizemin, libsizemax, E, 
                                                    tau_s, tau_p, npred, pred_start)
println("Starting Para Didi test")
validatevecs(para_didi, didi_para)
println("Test passed!")

