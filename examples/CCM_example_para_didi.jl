using CauseMap

# specify parameter ranges to test over
E_vals     = 2:10  # range to test of system dimensionality
tau_s_vals = 1:1 # range for lag length for manifold reconstruction
tau_p_vals = 0:15  # range to test for time lag of causal effect
#

# read in data
para, didi = load_example_data("ParaDidi")

# run analysis
makeoptimizationplots(para, didi,  
                                        E_vals, tau_s_vals, tau_p_vals, 
                                        "Para.", "Didi."; 
                                        nreps=10, left_E=2, left_tau_p=0,   # optional 
                                        right_E=7, right_tau_p=12, lagunit=.5, 
                                        unit="days", show_tau_s=false # optional
                                    )