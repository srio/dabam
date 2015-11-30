
from matplotlib import pylab as plt
import dabam
import numpy

do_plots = 0

#create DABAM objects to load experimental and two simulated profiles

dm = [dabam.dabam(),dabam.dabam(),dabam.dabam()]

#
# load input from DABAM
#

dabam_entry = 12

#dm = dabam.dabam()
dm[0].set_input_entryNumber(dabam_entry)
dm[0].load()
if do_plots:
    dm[0].set_input_plot("psd_h")
    dm[0].plot()
    dm[0].set_input_plot("acov_h")
    dm[0].plot()

#
# define inputs for fractal and Gaussian simulated profiles
#

cl = 0.05 #dm.autocorrelation_heights() # 0.03
beta =  1.1 # 0.9 #(-dm.powerlaw["hgt_pendent"])

mirror_length = (dm[0].sy[-1]-dm[0].sy[0])
npoints = (dm[0].zprof.size)
step = mirror_length / (npoints - 1)
random_seed = 65451
error_type   = 1  # 0=HEIGHTS, 1=SLOPES
if error_type == 0:
    rms = dm[0].stdev_profile_heights()
else:
    rms = dm[0].stdev_profile_slopes()

from error_profile import simulate_profile_1D
from dabam import autocovariance_1D

#profile_type = [0]  # 0=GAUSSIAN, 1=FRACTAL
idabam = 0

for iprofile_type in range(2):
    idabam += 1

    x,y = simulate_profile_1D(step=step,
                                mirror_length=mirror_length,
                                profile_type=iprofile_type,
                                correlation_length=cl,
                                power_law_exponent_beta=beta,
                                random_seed=random_seed,
                                error_type=error_type,
                                rms=rms,
                                rms_heights = 1e-6,
                              )


    # from correlated_profiles import rsgeng1D
    # # N = npoints # 1000
    # # mirror_length = (dm.sy[-1]-dm.sy[0]) # 1.0
    # # height_rms = 3e-9
    #
    # #cl = 0.03
    # x, y = rsgeng1D(npoints,mirror_length,height_rms,cl)
    slopes = numpy.gradient(y,x[1]-x[0])

    if do_plots:
        f1 = plt.figure(1)
        plt.plot(x,y)
        plt.plot(dm[0].sy,dm[0].zprof)
        plt.xlabel("Y")
        plt.ylabel("heights profile Z")
        plt.title("ENTRY=%d"%(dabam_entry))
        plt.show()

    #
    #reload simulated profile in dabam
    #
    #dm[iprofile_type] = dabam.dabam()
    dm[idabam].set_input_setDetrending(-1)
    dm[idabam].load_external_profile(x,slopes,type='slopes')
    if do_plots:
        dm[idabam].set_input_plot("psd_h")
        dm[idabam].plot()
        dm[idabam].set_input_plot("acov_h")
        dm[idabam].plot()

    print("mirror length:       original=%g, simulated=%g"%(dm[0].sy[-1]-dm[0].sy[0],x[-1]-x[0]))
    print("mirror step:         original=%g, simulated=%g"%(dm[0].sy[1] -dm[0].sy[0],x[1]-x[0]))
    print("mirror points:       original=%g, simulated=%g"%(dm[0].zprof.size,y.size))

    print("mirror slope error:  original=%g, simulated=%g"%( dm[0].stdev_profile_slopes(),   dm[idabam].stdev_profile_slopes()))
    print("mirror height error: original=%g, simulated=%g"%( dm[0].stdev_profile_heights(),  dm[idabam].stdev_profile_heights()))
    print("mirror beta:         original=%g, simulated=%g"%(-dm[0].powerlaw["hgt_pendent"],- dm[idabam].powerlaw["hgt_pendent"]))
    print("correlation length   original=%g, simulated=%g"%( dm[0].autocorrelation_heights(),dm[idabam].autocorrelation_heights()))



#
#dump file
#
output_file = "simulated_heights.dat"
f = open(output_file,'w')
for i in range(npoints):
    f.write("%g   %g   %g   %g\n"%(dm[0].sy[i],dm[0].zprof[i],dm[1].zprof[i],dm[2].zprof[i]))
f.close()

output_file = "simulated_psd.dat"
f = open(output_file,'w')
for i in range(len(dm[0].f)):
    f.write("%g   %g   %g   %g\n"%(dm[0].f[i],dm[0].psdHeights[i],dm[1].psdHeights[i],dm[2].psdHeights[i]))
f.close()