using CauseMap

para, didi = load_example_data("ParaDidi")

E_vals     = 2:6  # range to test of system dimensionality
tau_s_vals = 1:1 # range for lag length for manifold reconstruction
tau_p_vals = 0:5  # range to test for time lag of causal effect

# makesingleplot(para, didi ,
#                                         4, 1, 0,  
#                                        "x", "y")

makeoptimizationplots(repeat(para, outer=[1]), 
                                    repeat(didi, outer=[1]),  
                                        2:2, 1:1, 1:1, 
                                        "Para.", "Didi."; plot=false
                                    )

maxTSlen = 6
TSlens = 1:maxTSlen
outarray = fill(NaN, (maxTSlen, 2))
for TSlen in 2:2
    tic()
    makeoptimizationplots(repeat(para, outer=[TSlen]), 
                                        repeat(didi, outer=[TSlen]),  
                                            E_vals, tau_s_vals, tau_p_vals, 
                                            "Para.", "Didi."; plot=true
                                        )
    time = toc()
    outarray[TSlen, :] = [TSlen * length(para), time]
end
writedlm("TSresults.txt", outarray, '\t') 