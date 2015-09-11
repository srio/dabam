"""

aspheric_fit
   tests for performing aspheric and elliptical fits on height and slopes profiles

   inputs: reads the heights profiles from tmpHeights.da and the slopes profiles
           from tmpSlopes.dat, files produced by dabam.py without detrending, e.g.:
           python3 dabam.py 4 -D -1

    output: some plots



"""

__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2015"


import numpy
from scipy.optimize import curve_fit, leastsq

def func_aspheric(x, radius, kappa):
    #https://en.wikipedia.org/wiki/Aspheric_lens
    tmp = numpy.power(x,2) / radius / ( 1.0 + numpy.sqrt(1.0 - (1.0 + kappa) * ( numpy.power(x / radius,2) )))
    return tmp

def func_ellipse_slopes(x, p, q, theta):
    #
    # returns y'(x), the slopes of an ellipse defined by p,q, and theta
    #


    a = (p + q) / 2
    b = numpy.sqrt( p * q) * numpy.sin(theta)
    c = numpy.sqrt(a*a - b*b)

    epsilon = c / a

    # (x0,y0) are the coordinates of the center of the mirror
    # x0 = (p*p - q*q) / 4 / c
    x0 = (p - q) / 2 / epsilon
    y0 = -b * numpy.sqrt(1 - ((x0/a)**2))

    # the versor normal to the surface at the mirror center is -grad(ellipse)
    xnor = -2 * x0 / a**2
    ynor = -2 * y0 / b**2
    modnor = numpy.sqrt(xnor**2 + ynor**2)
    xnor /= modnor
    ynor /= modnor
    # tangent  versor is perpendicular to normal versor
    xtan =  ynor
    ytan = -xnor


    A = 1/b**2
    B = 1/a**2
    C = A

    CCC = numpy.zeros(11)
    #CCC[1] = A
    CCC[2] = B*xtan**2 + C*ytan**2
    CCC[3] = B*xnor**2 + C*ynor**2
    #CCC[4] = .0
    CCC[5] = 2*(B*xnor*xtan+C*ynor*ytan)
    #CCC[6] = .0
    #CCC[7] = .0
    CCC[8] = .0
    CCC[9] = 2*(B*x0*xnor+C*y0*ynor)
    CCC[10]= .0

    # ellipse implicit eq is c2 x^2 + c3 y^2 + c5 x y + c8 x + c9 y + c10 = 0
    # AA y^2 + BB y + CC = 0
    AA = CCC[3]
    BB = CCC[5]*x + CCC[9]
    CC = CCC[2]*x*x + CCC[8]*x + CCC[10]
    DD = BB*BB-4*AA*CC
    #yell = (-BB - numpy.sqrt(DD) )/(2*AA)
    #yellp = numpy.gradient(yell,(x[1]-x[0]))

    #calculate derivatives (primes P)
    BBP = CCC[5]
    CCP = 2*CCC[2]*x+CCC[8]
    DDP = 2*BB*BBP -4*AA*CCP
    ells = (-1/2/AA) * (BBP + DDP/2/numpy.sqrt(DD))

    return ells


def func_ellipse(x, p, q, theta):
    #
    # returns y(x), the heights of an ellipse defined by p,q, and theta
    #

    a = (p + q) / 2
    b = numpy.sqrt( p * q) * numpy.sin(theta)
    c = numpy.sqrt(a*a - b*b)

    epsilon = c / a

    # (x0,y0) are the coordinates of the center of the mirror
    # x0 = (p*p - q*q) / 4 / c
    x0 = (p - q) / 2 / epsilon
    y0 = -b * numpy.sqrt(1 - ((x0/a)**2))

    # print(">>>> func_ellipse: a=%f, b=%f, c=%f"%(a,b,c))
    # print(">>>> func_ellipse: x0=%f, y0=%f"%(x0,y0))

    # the versor normal to the surface at the mirror center is -grad(ellipse)
    xnor = -2 * x0 / a**2
    ynor = -2 * y0 / b**2
    modnor = numpy.sqrt(xnor**2 + ynor**2)
    xnor /= modnor
    ynor /= modnor
    # tangent  versor is perpendicular to normal versor
    xtan =  ynor
    ytan = -xnor


    A = 1/b**2
    B = 1/a**2
    C = A

    CCC = numpy.zeros(11)
    #CCC[1] = A
    CCC[2] = B*xtan**2 + C*ytan**2
    CCC[3] = B*xnor**2 + C*ynor**2
    #CCC[4] = .0
    CCC[5] = 2*(B*xnor*xtan+C*ynor*ytan)
    #CCC[6] = .0
    #CCC[7] = .0
    CCC[8] = .0
    CCC[9] = 2*(B*x0*xnor+C*y0*ynor)
    CCC[10]= .0

    # ellipse implicit eq is c2 x^2 + c3 y^2 + c5 x y + c8 x + c9 y + c10 = 0
    # AA y^2 + BB y + CC = 0
    AA = CCC[3]
    BB = CCC[5]*x + CCC[9]
    CC = CCC[2]*x*x + CCC[8]*x + CCC[10]

    yell1 = (-BB + numpy.sqrt(BB*BB-4*AA*CC) )/(2*AA)
    yell2 = (-BB - numpy.sqrt(BB*BB-4*AA*CC) )/(2*AA)

    return yell2

def func_ellipse_slopes_amparo(x, p, q, theta):
    #
    # returns y'(x), the slopes of an ellipse defined by p,q, and theta using the fromula in REF
    #
    a = (p + q) / 2
    b = numpy.sqrt( p * q) * numpy.sin(theta)
    c = numpy.sqrt(a*a - b*b)

    epsilon = c / a

    # (x0,y0) are the coordinates of the center of the mirror
    # x0 = (p*p - q*q) / 4 / c
    x0 = (p - q) / 2 / epsilon

    y0 = (-b * numpy.sqrt(1 - ((x0/a)**2)))
    mu = numpy.arctan( - b * x0 / numpy.sqrt(a*a - x0*x0) )
    #from Amparo Excel
    # y0 = -0.99994
    # mu = 0.003894
    y0 = -b*numpy.sqrt(1-x0*x0/a/a)
    alpha = numpy.arcsin(p/2/c*numpy.sin(numpy.pi-2*theta))
    mu = alpha - theta

    print(">>>> func_ellipse_slopes_amparo: a=%f, b=%f, c=%f"%(a,b,c))
    print(">>>> func_ellipse_slopes_amparo: a^2=%f, b^2=%f, F^2=%f"%(a**2,b**2,c**2))
    print(">>>> func_ellipse_slopes_amparo: x0=%f, y0=%f, mu=%f deg, mu=%f rad"%(x0,y0,180/numpy.pi*mu,mu))

    den = -(x*numpy.cos(mu))**2 - 2*numpy.cos(mu)*x*x0 + a**2 -x0**2
    bk = x0 + x*numpy.cos(mu)
    bk = bk / numpy.sqrt(den)
    ells = (numpy.cos(mu))**2 * (b/a) * bk + numpy.cos(mu) * numpy.sin(-mu)

    return ells

def func_ellipse_amparo(x, p, q, theta):
    #
    # returns y(x), the heights of an ellipse defined by p,q, and theta using the formula in REF
    #

    a = (p + q) / 2
    b = numpy.sqrt( p * q) * numpy.sin(theta)
    c = numpy.sqrt(a*a - b*b)

    epsilon = c / a

    # (x0,y0) are the coordinates of the center of the mirror
    # x0 = (p*p - q*q) / 4 / c
    x0 = (p - q) / 2 / epsilon
    y0 = (-b * numpy.sqrt(1 - ((x0/a)**2)))
    mu = numpy.arctan( - b * x0 / numpy.sqrt(a*a - x0*x0) )
    #from Amparo Excel
    # y0 = -0.99994
    # mu = 0.003894
    y0 = -b*numpy.sqrt(1-x0*x0/a/a)
    alpha = numpy.arcsin(p/2/c*numpy.sin(numpy.pi-2*theta))
    mu = alpha - theta

    # print(">>>> func_ellipse_slopes_amparo: a=%f, b=%f, c=%f"%(a,b,c))
    # print(">>>> func_ellipse_slopes_amparo: a^2=%f, b^2=%f, F^2=%f"%(a**2,b**2,c**2))
    # print(">>>> func_ellipse_slopes_amparo: x0=%f, y0=%f, mu=%f deg, mu=%f rad"%(x0,y0,180/numpy.pi*mu,mu))

    brk0 = ( (x*numpy.cos(mu)+x0)/a )**2
    #print("??? brk0: ",brk0)
    brk1 = b*b * (1 - brk0 )
    brk = -numpy.sqrt(brk1) -y0
    #print("??? pnc0: ",pnc0)
    pnc = numpy.cos(mu) * ( y0 + b*numpy.sqrt( 1 - x0*x0/a/a) )
    yell2 = numpy.cos(mu) * brk + numpy.sin(-mu) * x * numpy.cos(mu) + pnc
    return yell2


def cdf(sy, sz, method = 1 ):
    """
     cdf: Calculates the profile from the slope by simple integration

      INPUTS:
           sy - 1D array of (equally-spaced) lengths.
           sz - 1D array of slopes.
      KEYWORDS
           method : 0 use simple sum as integration method
                    1 use trapezoidal rule (default)
      RESTRICTIONS:
          the abscissas step must be sorted, but may not be constant

           1D array with cdf

    """

    zprof = sz*0.0
    if method == 0:
        steps = sy[0:sz.size-1]
        steps = numpy.concatenate(([0],steps))
        steps[0] = steps[1]
        steps.shape = -1
        steps = sy - steps
        zprof = numpy.cumsum(sz*steps)
    else:
        for i in range(sz.size):
          zprof[i]= numpy.trapz(sz[0:i+1], x = sy[0:i+1])

    return zprof


def main():
    #
    # y axis is horizontal
    # z axis is vertical
    #

    #
    #load height profile
    #
    input_file = "tmpHeights.dat"
    a = numpy.loadtxt(input_file)
    y0 = a[:,0]
    z0 = a[:,1]
    z0 -= z0.min()

    #
    #load slopes profile
    #
    input_file = "tmpSlopes.dat"
    a = numpy.loadtxt(input_file)
    yp0 = a[:,0]
    zp0 = a[:,1]

    L = y0[-1]-y0[0]
    print("Mirror data from file: %s :"%input_file)
    print("    Mirror length is: %.3f m"%L)
    N = y0.size
    print("    Mirror contains %d points"%N)

    slope_error_rms = zp0.std()
    print("    Mirror slope error RMS is  %.3f urad = %.3f arcsec"%(slope_error_rms*1e6,slope_error_rms*180/numpy.pi*3600))


    #
    #aspheric fit
    #

    # fit_method = 0   fit aspherical profile, use curve_fit on heights (NOT WORKING)
    # fit_method = 1   fit aspherical profile, use leastsq on heights (NOT WORKING)
    # fit_method = 2   fit aspherical profile, use bounded leastsq on heights (NOT WORKING)
    # fit_method = 3   fit elliptical profile
    # fit_method = 4   fit elliptical slopes

    fit_method = 3

    if fit_method == 0: # use curve_fit
        popt, pcov = curve_fit(func_aspheric, y0, z0)
        print("popt: ",repr(popt))
        print("pcov: ",repr(pcov))
        z2 = func_aspheric(y0, popt[0], popt[1])


    if fit_method == 1: # use lestsq
        fitfunc = lambda p, x: func_aspheric(x, p[0], p[1])
        #z2 = fitfunc([radius*1.2,0.0],y0)
        errfunc = lambda p, x, y: fitfunc(p, x) - y
        #z2 = errfunc([radius*1.2,0.0],y0,z0)


        #print(errfunc([radius,0.0],y0,z0))
        p0 = [radius, 0.0] # initial guess

        #popt, success = leastsq(errfunc, p0[:], args=(y0, z0))
        #print("popt: ",repr(popt))
        #print("success: ",repr(success))

        popt, cov_x, infodic, mesg, ier = leastsq(errfunc, p0, args=(y0, z0), full_output=True)
        print("Standard Least Squares fitting results:")
        print("p0:", p0)
        print("cov_x:", cov_x)
        print("infodic['nfev']:", infodic['nfev'])
        print("infodic['fvec']:", infodic['fvec'])
        print("infodic['fjac']:", infodic['fjac'])
        print("infodic['ipvt']:", infodic['ipvt'])
        print("infodic['qtf']:", infodic['qtf'])
        print("mesg:", mesg)
        print("ier:", ier)
        print(">>>>> popt: ",repr(popt))
        print("")

        z2 = func_aspheric(y0, popt[0], popt[1])

    if fit_method == 2: # use lestsq
        fitfunc = lambda p, x: func_aspheric(x, p[0], p[1])
        #z2 = fitfunc([radius*1.2,0.0],y0)
        errfunc = lambda p, x, y: fitfunc(p, x) - y
        #z2 = errfunc([radius*1.2,0.0],y0,z0)

        #print(errfunc([radius,0.0],y0,z0))
        p0 = [radius, -0.9999999999] # initial guess

        # https://github.com/jjhelmus/leastsqbound-scipy
        from leastsqbound import leastsqbound
        bounds = [(radius*0.5, radius*1.5), (-1.0, -0.9)]
        popt, cov_x, infodic, mesg, ier = leastsqbound(errfunc, p0, args=(y0,z0),
                                                    bounds=bounds,full_output=True)

        # print out results

        print("Bounded Least Squares fitting with no bounds results:")
        print("p0:", p0)
        print("cov_x:", cov_x)
        print("infodic['nfev']:", infodic['nfev'])
        print("infodic['fvec']:", infodic['fvec'])
        print("infodic['fjac']:", infodic['fjac'])
        print("infodic['ipvt']:", infodic['ipvt'])
        print("infodic['qtf']:", infodic['qtf'])
        print("mesg:", mesg)
        print("ier:", ier)
        print(">>>>> popt: ",repr(popt))
        print("")

        z2 = func_aspheric(y0, popt[0], popt[1])
        #z2 = func_aspheric(y0, popt[0], -0.99999996025050053)


    if fit_method == 3: # ellipse fitting heights

        ibounded = 1 #=0 curve_fit (no guess), 1=leastsq, 2=bounded

        print("======== Fitting heights =======")

        if ibounded == 0:
            print("======== Curve fitting without guess =======")
            popt, cov_x = curve_fit(func_ellipse, y0, z0, maxfev=10000)
        else:
            #dabam-4 (Amparo)
            p0 = 98.00
            q0 = 0.0775
            theta0 = 3.9e-3
            # #
            # #dabam-6 (Soleil) p=499.14 mm q= 4500 mm teta = 34.99 mrad
            # p0 = 499.14e-3
            # q0 = 4500e-3
            # theta0 = 34.99e-3
            # #
            # #dabam-19 #design parameter of ellipse: entrance arm: 420000mm; exit arm: 900mm; angle of incidence 3.052mrad
            # p0 = 420.0
            # q0 = 0.9
            # theta0 =  3.052e-3
            # #
            # #dabam-20 #design parameter of ellipse: entrance arm 9000mm; exit arm: 350mm; angle of incidence: 2.5deg
            # p0 = 9.0
            # q0 = 0.35
            # theta0 =  2.5*numpy.pi/180
            # #TODO:
            # zp0 = -zp0
            # #
            # #dabam-21 #design parameter of ellipse: entrance arm: 7500mm; exit arm: 2500mm; angle of incidence 0.6deg
            # p0 = 7.5
            # q0 = 2.5
            # theta0 =  0.6*numpy.pi/180

            p_guess = [p0,q0,theta0]
            z1 = func_ellipse_amparo(yp0, p_guess[0], p_guess[1], p_guess[2])

            # ishift = 0
            # if ishift:
            #     print(">>>>>>>> slope zp0[0], zp1[0], diff: ",zp0[0],zp1[1],zp0[0]-zp1[1])
            #     print(">>>>>>>>>>>>>>>>>>> shift value: ",(zp0-zp1)[1])
            #     p_guess[3]  = (zp0-zp1)[1]
            #     zp1 = func_ellipse(yp0, p_guess[0], p_guess[1], p_guess[2], p_guess[3])

            print("p0,q0,theta: ",p0,q0,theta0)

            fitfunc_ell_heights = lambda p, x: func_ellipse_amparo(x, p[0], p[1], p[2])
            errfunc_ell_heights = lambda p, x, y: fitfunc_ell_heights(p, x) - y

            if ibounded == 1:
                print("======== Least Squares fitting =======")
                popt, cov_x, infodic, mesg, ier = leastsq(errfunc_ell_heights, p_guess, args=(y0, z0),
                                                          full_output=True)
            else:
                print("======== Bounded Least Squares fitting =======")
                # https://github.com/jjhelmus/leastsqbound-scipy
                from leastsqbound import leastsqbound
                bounds = [ (p0*0.998,p0*1.002) , (q0*0.8,q0*1.2), (theta0*0.98,theta0*1.02)]
                popt, cov_x, infodic, mesg, ier = leastsqbound(errfunc_ell_heights, p_guess, args=(y0,z0),
                                                            bounds=bounds,full_output=True)


            print("cov_x:", cov_x)
            # print("infodic['nfev']:", infodic['nfev'])
            # print("infodic['fvec']:", infodic['fvec'])
            # print("infodic['fjac']:", infodic['fjac'])
            # print("infodic['ipvt']:", infodic['ipvt'])
            # print("infodic['qtf']:", infodic['qtf'])
            # print("mesg:", mesg)
            # print("ier:", ier)

            print(">>>>> p_guess:", p_guess)
            #zp1 = func_ellipse_slopes(yp0, p_guess[0], p_guess[1], p_guess[2], 0.0 ) #p_guess[3]) #TODO!!
            dd=numpy.concatenate( (y0, z1) ,axis=0).reshape(2,-1).transpose()
            outFile = "tmp_z1.dat"
            numpy.savetxt(outFile,dd)
            print ("File "+outFile+" written to disk:\n")



        print(">>>>> popt (p,q,theta): ",popt)

        z2 = func_ellipse(y0, popt[0], popt[1], popt[2])
        #amparo z2 = func_ellipse(y0, 98.0, 81.826e-3, 3.865*1e-3)
        if ibounded != 0: print ('Height error RMS z0-z1:             %.3f nm'%(1e9*(z0-z1).std()))
        print ('height error RMS z0-z2:             %.3f nm'%(1e9*(z0-z2).std()))
        dd=numpy.concatenate( (y0, z2) ,axis=0).reshape(2,-1).transpose()
        outFile = "tmp_z2.dat"
        numpy.savetxt(outFile,dd)
        print ("File "+outFile+" written to disk:\n")

    if fit_method == 4: # ellipse, slopes fit

        ibounded = 1 #=0 curve_fit (no guess), 1=leastsq, 2=bounded

        print("======== Fitting slopes =======")

        if ibounded == 0:
            print("======== Curve fitting without guess =======")
            popt, cov_x = curve_fit(func_ellipse_slopes, yp0, zp0, maxfev=10000)
        else:
            ishift = 0

            #dabam-4 (Amparo)
            p0 = 98.00
            q0 = 0.0775
            theta0 = 3.9e-3

            # # # #
            # #dabam-6 (Soleil) p=499.14 mm q= 4500 mm teta = 34.99 mrad
            # p0 = 499.14e-3
            # q0 = 4500e-3
            # theta0 = 34.99e-3
            # ishift = 1
            # #
            # #
            # # #
            # #dabam-19 #design parameter of ellipse: entrance arm: 420000mm; exit arm: 900mm; angle of incidence 3.052mrad
            # p0 = 420.0            #WRONG INPUTS?
            # q0 = 0.9              #WRONG INPUTS?
            # theta0 =  3.052e-3    #WRONG INPUTS?
            # # yp0 = yp0 - yp0[int(yp0.size/2)]
            #
            # #
            # #dabam-20 #design parameter of ellipse: entrance arm 9000mm; exit arm: 350mm; angle of incidence: 2.5deg
            # p0 = 9.0
            # q0 = 0.35
            # theta0 =  2.5*numpy.pi/180
            # #TODO:
            # yp0 = yp0 - yp0[int(yp0.size/2)]
            # zp0 = -zp0
            #
            # #
            # #dabam-21 #design parameter of ellipse: entrance arm: 7500mm; exit arm: 2500mm; angle of incidence 0.6deg
            # p0 = 7.5
            # q0 = 2.5
            # theta0 =  0.6*numpy.pi/180
            # ishift = 1

            p_guess = [p0,q0,theta0]
            zp1 = func_ellipse_slopes(yp0, p_guess[0], p_guess[1], p_guess[2])

            if ishift:
                zp0 = zp0 + (zp1[0]-zp0[0])

            print("p0,q0,theta: ",p0,q0,theta0)

            fitfunc_ell_slopes = lambda p, x:   func_ellipse_slopes(x, p[0], p[1], p[2])
            errfunc_ell_slopes = lambda p, x, y: fitfunc_ell_slopes(p, x) - y

            if ibounded == 1:
                print("======== Least Squares fitting =======")
                popt, cov_x, infodic, mesg, ier = leastsq(errfunc_ell_slopes, p_guess, args=(yp0, zp0),
                                                          full_output=True)
            elif ibounded == 2:
                print("======== Bounded Least Squares fitting =======")
                # https://github.com/jjhelmus/leastsqbound-scipy
                from leastsqbound import leastsqbound
                bounds = [ (p0*0.998,p0*1.002) , (q0*0.8,q0*1.2), (theta0*0.98,theta0*1.02)]
                popt, cov_x, infodic, mesg, ier = leastsqbound(errfunc_ell_slopes, p_guess, args=(yp0,zp0),
                                                            bounds=bounds,full_output=True)


            print("cov_x:", cov_x)
            # print("infodic['nfev']:", infodic['nfev'])
            # print("infodic['fvec']:", infodic['fvec'])
            # print("infodic['fjac']:", infodic['fjac'])
            # print("infodic['ipvt']:", infodic['ipvt'])
            # print("infodic['qtf']:", infodic['qtf'])
            # print("mesg:", mesg)
            # print("ier:", ier)

            print(">>>>> p_guess:", p_guess)
            #zp1 = func_ellipse_slopes(yp0, p_guess[0], p_guess[1], p_guess[2])
            dd=numpy.concatenate( (yp0, zp1) ,axis=0).reshape(2,-1).transpose()
            outFile = "tmp_zp1.dat"
            numpy.savetxt(outFile,dd)
            print ("File "+outFile+" written to disk:\n")



        print(">>>>> popt (p,q,theta): ",popt)
        zp2 = func_ellipse_slopes(yp0, popt[0], popt[1], popt[2])
        #amparo  zp2 = func_ellipse_slopes(yp, 98.0, 82.0424e-3, 3.8754e-3, 0.0)

        if ibounded != 0: print ('Slope error RMS zp0-zp1:             %.3f urad'%(1e6*(zp0-zp1).std()))
        print ('Slope error RMS zp0-zp2:             %.3f urad'%(1e6*(zp0-zp2).std()))
        dd=numpy.concatenate( (yp0, zp2) ,axis=0).reshape(2,-1).transpose()
        outFile = "tmp_zp2.dat"
        numpy.savetxt(outFile,dd)
        print ("File "+outFile+" written to disk:\n")



    #
    #plots
    #
    do_plots = 1
    if do_plots:
        #
        #plots
        #
        from matplotlib import pylab as plt

        if fit_method == 4: #slopes
            f1 = plt.figure(1)
            ax = plt.gca()
            ax.get_xaxis().get_major_formatter().set_useOffset(False)
            ax.get_yaxis().get_major_formatter().set_useOffset(False)
            plt.plot(yp0*1e3,zp0*1e6)
            if ibounded != 0: plt.plot(yp0*1e3,zp1*1e6)
            plt.plot(yp0*1e3,zp2*1e6)
            plt.title("slopes data (blue), starting (green) and optimized (red) ellipse")
            plt.xlabel("Y [mm]")
            plt.ylabel("Zp [urad]")

            f2 = plt.figure(2)
            ax = plt.gca()
            ax.get_xaxis().get_major_formatter().set_useOffset(False)
            ax.get_yaxis().get_major_formatter().set_useOffset(False)
            plt.plot(yp0*1e3,(zp0-zp2)*1e6)
            plt.title("residual slopes")
            plt.xlabel("Y [mm]")
            plt.ylabel("Zp [urad]")
        elif fit_method == 3: #heights
            f1 = plt.figure(1)
            ax = plt.gca()
            ax.get_xaxis().get_major_formatter().set_useOffset(False)
            ax.get_yaxis().get_major_formatter().set_useOffset(False)
            plt.plot(y0*1e3,z0*1e6)
            if ibounded != 0: plt.plot(y0*1e3,z1*1e6)
            plt.plot(y0*1e3,z2*1e6)
            plt.title("height data (blue), starting (green) and optimized (red) ellipse")
            plt.xlabel("Y [mm]")
            plt.ylabel("Zp [um]")

            f2 = plt.figure(2)
            ax = plt.gca()
            ax.get_xaxis().get_major_formatter().set_useOffset(False)
            ax.get_yaxis().get_major_formatter().set_useOffset(False)
            plt.plot(y0*1e3,(z0-z2)*1e9)
            plt.title("residual heights")
            plt.xlabel("Y [mm]")
            plt.ylabel("Z [nm]")



        plt.show()

#
# main program
#
if __name__ == '__main__':
    main()