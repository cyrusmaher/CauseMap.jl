require("CauseMap.jl")

libsizemin = 8   # minimum library size
libsizemax = 58 # maximum library size
pred_start = 13 # first point to predict
npred      = 60  # number of time series points to predict (larger is better)

E_vals     = 2:10  # range to test of system dimensionality
tau_s_vals = 1:1 # range for lag length for manifold reconstruction
tau_p_vals = 0:15  # range to test for time lag of causal effect

ds   = readdlm("vr.raw_fixed.txt", '\t', Float64, has_header=false)  # read in data
para = ds[:, 2]
didi = ds[:, 3]

# run analysis
makeoptimizationplots(para, didi, libsizemin, libsizemax, E_vals, 
			tau_s_vals, tau_p_vals, npred, pred_start,
			"Para.", "Didi."; nreps=1, left_E=2, left_tau_p=0, 
                                    right_E=7, right_tau_p=12, lagunit=.5, unit="days")