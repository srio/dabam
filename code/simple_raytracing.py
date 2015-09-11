"""

simple_raytracing
   performs optics calculations using raytracing

   inputs: reads the heights profiles from tmpHeights.dat
           file produced by dabam.py with detrending, e.g.:
           python3 dabam.py 4

    output: some plots



"""

__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2015"


import numpy
from matplotlib import pylab as plt

def main():
    #
    # y axis is horizontal
    # z axis is vertical
    #

    #
    #define focal distances
    #
    p = 30.0
    q = 10.0
    theta_grazing = 3e-3

    #
    #compute mirror radius
    #
    R = 2 / numpy.sin(theta_grazing) / (1/p + 1/q)
    print("Mirror radius of curvature set to: %.3f m (p=%.3f m, q=%.3f m, theta=%.2f mrad))"%(R,p,q,theta_grazing*1e3))



    #
    #load height profile
    #
    input_file = "tmpHeights.dat"
    a = numpy.loadtxt(input_file)
    hy0 = a[:,0]
    hz0 = a[:,1]

    #
    #interpolate to increase the number of points ans statistics
    #
    do_interpolate = 0
    if do_interpolate:
        mirror_length = hy0.max() - hy0.min()
        hy = numpy.linspace(hy0.min(),hy0.max(),10000)
        hz = numpy.interp(hy,hy0,hz0)
    else:
        hy = hy0
        hz = hz0


    L = hy[-1]-hy[0]
    print("Mirror data from file: %s :"%input_file)
    print("    Mirror length is: %.3f m"%L)
    print("    Mirror aperture is: %.3f um"%(1e6*L*numpy.sin(theta_grazing)))
    N = hy.size
    print("    Mirror contains %d points"%N)


    #
    #compute slopes
    #
    sz = numpy.gradient(hz,(hy[1]-hy[0]))
    slope_errors_rms = sz.std()
    print("    Mirror slope error RMS is  %.3f urad = %.3f arcsec"%(slope_errors_rms*1e6,slope_errors_rms*180/numpy.pi*3600))

    #
    #project on optical axis
    #
    hy_projected = hy * numpy.sin(theta_grazing)

    #angle with respect to Y axis
    theta_incident = hy_projected / p


    #reflection law on the flat mirror
    theta_reflection = theta_incident

    #apply focusing
    theta_reflection =  theta_reflection - 2 * hy / R

    #apply slope error
    theta_reflection =  theta_reflection - 2 * sz


    #compute coordinates at the image position
    image_z = hy_projected + q * theta_reflection

    #
    #image histogram
    #
    image_histogram, bin_edges = numpy.histogram(image_z,bins=51)
    bin_centers = bin_edges[0:-1]
    bin_centers += (bin_edges[1] - bin_centers[0])

    # dump to file
    outFile = "tmpImage.dat"
    if outFile != "":
        dd = numpy.concatenate( (bin_centers.reshape(-1,1), image_histogram.reshape(-1,1)),axis=1)
        dd[:,0] *= -1e6 # in microns, inverted to agree with shadow
        dd[:,1] /= dd[:,1].max()
        numpy.savetxt(outFile,dd)
        print("File written to disk: %s"%outFile)

    #
    #CALCULATE fwhm
    #
    tt = numpy.where(image_histogram>=max(image_histogram)*0.5)
    if image_histogram[tt].size > 1:
        binSize = bin_edges[1]-bin_edges[0]
        fwhm = 1e6*binSize*(tt[0][-1]-tt[0][0])
        print('Image fwhm: %.3f um'%(fwhm))
    fwhm_theory = 2.35*2*slope_errors_rms*q
    print('Image 2*slope_error_rms*q: %.3f um'%(fwhm_theory*1e6))



    #
    #plots
    #
    do_plots = 0
    if do_plots:
        f1 = plt.figure(3)
        plt.plot(1e6*bin_centers,image_histogram)
        plt.xlim((-50,50))
        plt.title("image histogram FWHM = %.3f um, theory %.3f"%(fwhm,fwhm_theory*1e6))
        plt.xlabel("Y [um]")
        plt.show()

#
# main program
#
if __name__ == '__main__':
    main()