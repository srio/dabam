"""

dabam: (dataBase for metrology)
       python tools for processing remote files containing the results
       of metrology measurements on X-ray mirrors

       functions: 
             cdf (calculate cumulative distribution function)
             psd (calculate power spectral density)
             write_shadowSurface (writes file with a mesh for SHADOW)
 
       MODIFICATION HISTORY:
           20150828 srio@esrf.eu, written
           20131109 srio@esrf.eu, added command line arguments, access metadata

"""

__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2015"


import numpy
from matplotlib import pylab as plt

#lensF is in fact 2*F
def goFromTo(source, image, distance=1.0, lensF=None, slopeError=None, wavelength=1e-10):
    distance = numpy.array(distance)
    x1 = numpy.outer(source,numpy.ones(image.size))
    x2 = numpy.outer(numpy.ones(source.size),image)
    r = numpy.sqrt( numpy.power(x1-x2,2) + numpy.power(distance,2) )
    # add lens at the image plane
    if lensF != None:
      r = r - numpy.power(x1-x2,2)/lensF
    if slopeError != None:
      r = r + 2 * slopeError
    wavenumber = numpy.pi*2/wavelength
    return numpy.exp(1.j * wavenumber *  r)


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
    F = 1 / (1/p + 1/q)
    print("Mirror radius of curvature set to: %.3f m (p=%.3f m, q=%.3f m, theta=%.2f mrad))"%(R,p,q,theta_grazing*1e3))
    print("Mirror focal length set to: %.3f m "%(F))


    #
    #load height profile
    #
    input_file = "tmpProfile.dat"
    a = numpy.loadtxt(input_file)
    hy0 = a[:,0]
    hz0 = a[:,1]

    #
    #interpolate to increase the number of points ans statistics
    #
    do_interpolate = 0
    if do_interpolate:
        mirror_length = hy0.max() - hy0.min()
        hy = numpy.linspace(hy0.min(),hy0.max(),2000)
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





    # # dump to file
    # dd = numpy.concatenate( (bin_centers.reshape(-1,1), image_histogram.reshape(-1,1)),axis=1)
    # outFile = "tmpImage.dat"
    # dd[:,0] *= -1e6 # in microns, inverted to agree with shadow
    # dd[:,1] /= dd[:,1].max()
    # numpy.savetxt(outFile,dd)
    sourcepoints = 1000
    slitpoints = 1000
    detpoints = 1000

    wavelength   =   1e-10
    aperture_diameter   =   2 * hy_projected.max()

    airy_disk_theta = 1.22 * wavelength / aperture_diameter


    detector_size = 50 * airy_disk_theta * q
    fwhm_theory = 2 * 2.35 * slope_errors_rms * q


    print("aperture _diameter = %f um "%(aperture_diameter*1e6))
    print("detector_size = %f um"%(detector_size*1e6))
    print("FWHM theory (2 sigma q) = %f um"%(fwhm_theory*1e6))
    print("Airy disk is: %f urad = %f um"%(airy_disk_theta*1e6,airy_disk_theta*q*1e6))
    if airy_disk_theta*q >= detector_size:
        detector_size = 5 * airy_disk_theta * q
        print("detector_size NEW = %f um"%(detector_size*1e6))

    position1x = numpy.linspace(0,0,sourcepoints)
    position2x = numpy.linspace(-aperture_diameter/2,aperture_diameter/2,slitpoints)
    position3x = numpy.linspace(-detector_size/2,detector_size/2,detpoints)

    sz_projected_interpolated = numpy.interp(position2x, hy, sz * numpy.sin(theta_grazing) )
    # sz_projected_interpolated = None

    # fields12 = goFromTo(position1x,position2x,q, wavelength=wavelength, lensF=2*F)

    fields12 = goFromTo(position1x,position2x,p, lensF=2*F, slopeError=sz_projected_interpolated, wavelength=wavelength)
    fields23 = goFromTo(position2x,position3x,q, lensF=None,wavelength=wavelength)
    # from 1 to 3, matrix multiplication
    fields13 = numpy.dot(fields12,fields23)


    print ("Shape of fields12, fields23, fields13: ",fields12.shape,fields23.shape,fields13.shape)

    #prepare results
    fieldComplexAmplitude = numpy.dot(numpy.ones(sourcepoints),fields13)
    print ("Shape of Complex U: ",fieldComplexAmplitude.shape)
    print ("Shape of position1x: ",position1x.shape)
    fieldIntensity = numpy.power(numpy.abs(fieldComplexAmplitude),2)
    fieldPhase = numpy.arctan2(numpy.real(fieldComplexAmplitude), \
                               numpy.imag(fieldComplexAmplitude))


    #
    # write spec formatted file
    #
    # out_file = "simple_physicaloptics.spec"
    # f = open(out_file, 'w')
    # header="#F %s \n\n#S  1 fresnel-kirchhoff diffraction integral\n#N 3 \n#L X[m]  intensity  phase\n"%out_file
    #
    # f.write(header)
    #
    # for i in range(detpoints):
    #    out = numpy.array((position2x[i], fieldIntensity[i], fieldPhase[i]))
    #    f.write( ("%20.11e "*out.size+"\n") % tuple( out.tolist())  )
    #
    # f.close()
    # print ("File written to disk: %s"%out_file)

    fieldIntensity /= fieldIntensity.max()
    tmpAbscissas = position3x * 1e6
    dd=numpy.concatenate( (tmpAbscissas, fieldIntensity) ,axis=0).reshape(2,-1).transpose()

    outFile = "tmpPhysicalOptics.dat"
    numpy.savetxt(outFile,dd)
    print ("File "+outFile+" written to disk:\n")



    #
    #plots
    #
    do_plots = 0
    if do_plots:
        #
        #plots
        #
        from matplotlib import pylab as plt

        plt.figure(1)
        plt.plot(position3x*1e6,fieldIntensity)
        plt.title("Fresnel-Kirchhoff Diffraction")
        plt.xlabel("X [um]")
        plt.ylabel("Intensity [a.u.]")
        plt.show()


        # f1 = plt.figure(1)
        # plt.plot(hy*1e3,hz*1e6)
        # plt.title("heights")
        # plt.xlabel("Y [mm]")
        # plt.ylabel("Z [nm]")
        #
        # f2 = plt.figure(2)
        # plt.plot(hy*1e3,sz*1e6)
        # plt.title("slopes")
        # plt.xlabel("Y [mm]")
        # plt.ylabel("Z' [urad]")


        # f3 = plt.figure(3)
        # #plt.plot(1e6*bin_edges[0:-1],image_histogram)
        # plt.plot(1e6*bin_centers,image_histogram)
        # plt.xlim((-50,50))
        # plt.title("image histogram FWHM = %.3f um, theory %.3f"%(fwhm,fwhm_theory*1e6))
        # plt.xlabel("Y [um]")
        #
        # plt.show()

#
# main program
#
if __name__ == '__main__':
    main()