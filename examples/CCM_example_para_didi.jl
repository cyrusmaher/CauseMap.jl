using CauseMap

E_vals     = 2:10  # range to test of system dimensionality
tau_s_vals = 1:1 # range for lag length for manifold reconstruction
tau_p_vals = 0:15  # range to test for time lag of causal effect

### You can specify these if you want to finer control
# libsizemin = 8   # minimum library size
# libsizemax = 58 # maximum library size
# pred_start = 13 # first point to predict
# npred      = 60  # number of time series points to predict (larger is better)

# read in data
para, didi = load_example_data("ParaDidi")

# run analysis
makeoptimizationplots(para, didi,  E_vals, tau_s_vals, tau_p_vals, "Para.", "Didi.")  # make diagnostic plot
                                                                # ; nreps=1, left_E=2, left_tau_p=0,   # tune output by adding args
                                                                # right_E=7, right_tau_p=12, lagunit=.5, unit="days", 
                                                                # libsizemin=libsizemin, libsizemax=libsizemax,
                                                                # npred=npred, pred_start=pred_start)