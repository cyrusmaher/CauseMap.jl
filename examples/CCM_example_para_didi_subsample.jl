using CauseMap

function calcFullseries(vec1, vec2, E, tau_s, tau_p, npred, pred_start, lib_start; nboots=0)

    shadowmat_dict_vec1, distmat_dict_vec1 = precalc_manif_dists(3:3, 1:1, vec1)
    shadowmat_dict_vec2, distmat_dict_vec2 = precalc_manif_dists(3:3, 1:1, vec2)

    librange12, yval_12 = calcCCM(vec1, vec2, 
                                                        shadowmat_dict_vec1, distmat_dict_vec1, 
                                                        E, tau_s, tau_p; nboots=nboots,
                                                        libsizemin=libsizemin, libsizemax=libsizemax,
                                                        npred=npred, pred_start=pred_start)

    librange21, yval_21 = calcCCM(vec2, vec1, 
                                                        shadowmat_dict_vec2, distmat_dict_vec2, 
                                                        E, tau_s, tau_p; nboots=nboots,
                                                        libsizemin=libsizemin, libsizemax=libsizemax,
                                                        npred=npred, pred_start=pred_start)

    return librange12, librange21, yval_12, yval_21
end


function calcshortseries(vec1, vec2, E, tau_s, tau_p, npred, pred_start, lib_start, libsizemin, libsizemax, factor::Int64; nboots=0)
    npred = int(npred/factor)
    libsizemax = int(libsizemax/factor)
    xval = libsizemin:libsizemax
    vec1 = vec1[1:factor:end]
    vec2 = vec2[1:factor:end]

    shadowmat_dict_vec1, distmat_dict_vec1 = precalc_manif_dists(E_vals, tau_s_vals, vec1)
    shadowmat_dict_vec2, distmat_dict_vec2 = precalc_manif_dists(E_vals, tau_s_vals, vec2)

    lib12, res12 = calcCCM(vec1, vec2,
                                            shadowmat_dict_vec1, distmat_dict_vec1, 
                                            E, 1, 0; 
                                            nboots=nboots, lib_start=lib_start,
                                            libsizemin=libsizemin, libsizemax=libsizemax,
                                            npred=npred, pred_start=pred_start)

    lib21, res21 = calcCCM(vec2, vec1,
                                            shadowmat_dict_vec2, distmat_dict_vec2, 
                                            E, 1, 0; 
                                            nboots=nboots, lib_start=lib_start,
                                            libsizemin=libsizemin, libsizemax=libsizemax,
                                            npred=npred, pred_start=pred_start)
    return xval, res12, res21

end


E = 3
tau_s = 1
tau_p = 0
lib_start = 1
libsizemin = 5 
libsizemax = 58
pred_start = 13
npred = 60
E_vals = E:E
tau_s_vals = tau_s:tau_s
nboots = 0
nboots2 = 50

ds   = readdlm("vr.raw_fixed.txt", '\t', Float64, header=false)  # read in data
para = ds[:, 2]
didi = ds[:, 3]

librange12, librange21, yval_12, yval_21 = calcFullseries(para, didi, E, tau_s, tau_p, npred, pred_start, lib_start; nboots=nboots)

xval_f2, yval_12_f2, yval_21_f2 = calcshortseries(para, didi, E, tau_s, 
                                                                            tau_p, npred, pred_start, 
                                                                            lib_start, libsizemin, libsizemax, 2; nboots=nboots)

xval_f3, yval_12_f3, yval_21_f3 = calcshortseries(para, didi, E, tau_s, 
                                                                            tau_p, npred, pred_start, 
                                                                            lib_start, libsizemin, libsizemax, 3; nboots=nboots)
fac = 3
left_E = 2
right_E = 6
right_tau_p = int(12/fac)
E_vals = left_E:right_E
tau_p_vals = 0:right_tau_p


# makeoptimizationplots(para[1:fac:end], didi[1:fac:end], libsizemin, int(libsizemax/fac), E_vals, 
#             tau_s_vals, tau_p_vals, int(npred/fac), pred_start,
#             "Para.", "Didi."; nreps=1, left_E=3, left_tau_p=0, 
#                                     right_E=right_E, right_tau_p=right_tau_p, lagunit=.5 * fac, unit="days", nboots=nboots, legend=true)

using PyCall
@pyimport matplotlib.pyplot as plt
@pyimport seaborn as sns
sns.set_style("white")
sns.set_context("poster")

plt.plot(librange12, yval_12, color="k", alpha=.7, label="Full (60 points)")
plt.plot(librange21, yval_21, color="k", alpha=.5)
plt.plot(xval_f2, yval_12_f2, color="b", label="1/2 thinned (30 points)", alpha=.9)
plt.plot(xval_f2, yval_21_f2, color="b",  alpha=.5)
plt.plot(xval_f3, yval_12_f3, color="r", label="1/3 thinned (20 points)",  alpha=.9)
plt.plot(xval_f3, yval_21_f3, color="r",  alpha=.5)
plt.xlabel("Library size")
plt.ylabel("\$\\rho_{ccm}\$", fontsize=40)
plt.locator_params(axis="y", nbins=4)
plt.legend(loc=0)
plt.show()

