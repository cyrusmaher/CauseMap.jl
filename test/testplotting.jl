using CauseMap

nboots = 50  # use bootstrap libraries instead of sliding windows (reduces effect of secular trends)

E     = 4  
tau_s = 1
tau_p = 0


ds   = readdlm("vr.raw_fixed.txt", '\t', Float64, header=false)  # read in data
para = ds[:, 2]
didi = ds[:, 3]

println("Starting bootstrap test")
makesingleplot(para, didi, 
                            E, tau_s, tau_p,
                            "para", "didi"; nboots=nboots, plot=false)


E_vals     = 2:10  # range to test of system dimensionality
tau_s_vals = 1:1 # range for lag length for manifold reconstruction
tau_p_vals = 0:15  # range to test for time lag of causal effect

println("Starting optimization plot test")
makeoptimizationplots(para, didi,
                                        E_vals, tau_s_vals, tau_p_vals,
                                        "Para.", "Didi."; nreps=1, plot=false)