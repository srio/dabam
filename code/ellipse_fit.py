"""

ellipse_fit   tests for performing elliptical fits on height and slopes profiles

    inputs: reads the heights profiles from dabam database without detrending,

    output: some statistical figures, files and plots


    STATUS: under development!

    USAGE: check main program

"""

__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2015"


import numpy
from scipy.optimize import curve_fit, leastsq


def func_ellipse_slopes_dabam(x, p, q, theta, shift):
    #
    # returns y'(x), the slopes of an ellipse defined by p,q, and theta
    #

    a = (p + q) / 2
    b = numpy.sqrt( p * q) * numpy.sin(theta)
    c = numpy.sqrt(a*a - b*b)

    epsilon = c / a

    # (x0,y0) are the coordinates of the center of the mirrorc
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

    CCC[2] = B*xtan**2 + C*ytan**2
    CCC[3] = B*xnor**2 + C*ynor**2

    CCC[5] = 2*(B*xnor*xtan+C*ynor*ytan)

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

    return ells+shift


def func_ellipse_heights_dabam(x1, p, q, theta, hshift, vshift):
    #
    # returns y(x), the heights of an ellipse defined by p,q, and theta
    #

    # x = x1
    # vshift = 0.0
    x = x1 + hshift*1e-9

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

    yell1 = (-BB + numpy.sqrt(BB*BB-4*AA*CC) )/(2*AA)
    yell2 = (-BB - numpy.sqrt(BB*BB-4*AA*CC) )/(2*AA)

    return yell2+vshift*1e-9

def func_ellipse_slopes_amparo(x, p, q, theta, shift):
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

    # y0 = (-b * numpy.sqrt(1 - ((x0/a)**2)))
    # mu = numpy.arctan( - b * x0 / numpy.sqrt(a*a - x0*x0) )
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

    return ells+shift

def func_ellipse_heights_amparo(x1, p, q, theta, hshift, vshift):
    #
    # returns y(x), the heights of an ellipse defined by p,q, and theta using the formula in REF
    #

    x = x1 + hshift*1e-9

    a = (p + q) / 2
    b = numpy.sqrt( p * q) * numpy.sin(theta)
    c = numpy.sqrt(a*a - b*b)

    epsilon = c / a

    # (x0,y0) are the coordinates of the center of the mirror
    # x0 = (p*p - q*q) / 4 / c
    x0 = (p - q) / 2 / epsilon
    # y0 = (-b * numpy.sqrt(1 - ((x0/a)**2)))
    # mu = numpy.arctan( - b * x0 / numpy.sqrt(a*a - x0*x0) )
    #from Amparo Excel
    # y0 = -0.99994
    # mu = 0.003894
    y0 = -b*numpy.sqrt(1-x0*x0/a/a)
    alpha = numpy.arcsin(p/2/c*numpy.sin(numpy.pi-2*theta))
    mu = alpha - theta

    brk0 = ( (x*numpy.cos(mu)+x0)/a )**2
    brk1 = b*b * (1 - brk0 )
    brk = -numpy.sqrt(brk1) -y0
    pnc = numpy.cos(mu) * ( y0 + b*numpy.sqrt( 1 - x0*x0/a/a) )
    yell2 = numpy.cos(mu) * brk + numpy.sin(-mu) * x * numpy.cos(mu) + pnc
    return yell2+vshift*1e-9


def func_ellipse_heights_xianbo(x1, S1, S2, TH, hshift, vshift): # OFFSET, XC):
    #
    # Hi Manuel,
    # I have checked the mirror 20 and 21. I got quite good elliptical fit with my Igor fitting function. The figure error results are:
    # Mirror 20: sigma_s = 0.42 urad, sigma_h = 2.06 nm
    # Mirror 21: sigma_s = 0.50 urad, sigma_h = 6.62 nm
    #
    # The 1-D elliptical fitting function I used is:
    #
    # Variables: OFFSET, S1, S2 and TH, XC (x offset)
    #
    # y = OFFSET+ cos(TH-asin((S1*sin(2.0*TH))/sqrt(S1^2.0+S2^2.0+2*S1*S2*cos(2*TH))))*(sqrt(S1*S2*(1.0-(-2.0*S2*sqrt((S2+S1*cos(2*TH))^2.0/(S1^2.0+S2^2.0+2.0*S1*S2*cos(2.0*TH)))+sqrt(S1^2.0+S2^2.0+2.0*S1*S2*cos(2.0*TH)))^2.0/(S1+S2)^2.0)*sin(TH)^2.0)-sqrt(S1*S2*(1.0-(1.0/(S1+S2)^2.0)*(-2.0*S2*sqrt((S2+S1*cos(2.0*TH))^2.0/(S1^2.0+S2^2.0+2.0*S1*S2*cos(2.0*TH)))+sqrt(S1^2.0+S2^2.0+2.0*S1*S2*cos(2.0*TH))+2.0*(X-XC)*cos(TH-asin((S1*sin(2.0*TH))/sqrt(S1^2.0+S2^2.0+2.0*S1*S2*cos(2.0*TH)))))^2.0)*sin(TH)^2.0)+(X-XC)*sin(TH-asin((S1*sin(2.0*TH))/sqrt(S1^2.0+S2^2.0+2.0*S1*S2*cos(2*TH)))))

    X = x1 + hshift*1e-9
    XC = 0.0
    y =         numpy.cos(TH-numpy.arcsin((S1*numpy.sin(2.0*TH))/numpy.sqrt(S1**2.0+S2**2.0+2*S1*S2*numpy.cos(2*TH))))*\
                (numpy.sqrt(S1*S2*(1.0-(-2.0*S2*numpy.sqrt((S2+S1*numpy.cos(2*TH))**2.0/ \
                (S1**2.0+S2**2.0+2.0*S1*S2*numpy.cos(2.0*TH)))+numpy.sqrt(S1**2.0+S2**2.0+ \
                2.0*S1*S2*numpy.cos(2.0*TH)))**2.0/(S1+S2)**2.0)*numpy.sin(TH)**2.0)- \
                numpy.sqrt(S1*S2*(1.0-(1.0/(S1+S2)**2.0)*(-2.0*S2*numpy.sqrt((S2+S1*numpy.cos(2.0*TH))**2.0/ \
                (S1**2.0+S2**2.0+2.0*S1*S2*numpy.cos(2.0*TH)))+numpy.sqrt(S1**2.0+S2**2.0+ \
                2.0*S1*S2*numpy.cos(2.0*TH))+2.0*(X-XC)*numpy.cos(TH-numpy.arcsin((S1*numpy.sin(2.0*TH))/ \
                numpy.sqrt(S1**2.0+S2**2.0+2.0*S1*S2*numpy.cos(2.0*TH)))))**2.0)*numpy.sin(TH)**2.0)+ \
                (X-XC)*numpy.sin(TH-numpy.arcsin((S1*numpy.sin(2.0*TH))/numpy.sqrt(S1**2.0+S2**2.0+2.0*S1*S2*numpy.cos(2*TH)))))

    return y+vshift*1e-9


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


def ellipse_fit(entry_number = 21, fit_method=1, fit_function='dabam', ibounded=1,do_plots=1):
    """

    :param entry_number: dabam entry number
    :param fit_method: fit elliptical heights (0), fit elliptical slopes (1)
    :param fit_function: 'dabam' , 'amparo', 'xianbo'
    :param ibounded: 0 curve_fit (no guess), 1=leastsq, 2=bounded
    :param do_plot: 1 display plots
    :return:
    """

    #
    #get profiles
    #
    import dabam
    dm = dabam.dabam()
    dm.inputs["entryNumber"] = entry_number
    dm.inputs["setDetrending"] = -1
    dm.load()

    #
    #
    #
    p0 = dm.h["ELLIPSE_DESIGN_P"]
    q0 = dm.h["ELLIPSE_DESIGN_Q"]
    theta0 = dm.h["ELLIPSE_DESIGN_THETA"]

    # y axis is horizonta, z axis is vertical
    y0  = dm.sy
    z0  = dm.zprof
    zp0 = dm.sz

    if fit_method == 0: # ellipse fitting heights

        hshift = 0.0 #nm
        vshift = 0.0 # nm

        #preprocessors
        if entry_number == 4:
            z0 -= z0.min()
        elif entry_number == 6:
            imin = z0.argmin()
            z0 -= z0[imin]
            y0 -= y0[imin]
            y0 *= -1.0
        elif entry_number == 19:
            pass
        elif entry_number == 20:
            hshift = -15e3
        elif entry_number == 21:
            y0 -= y0[y0.size/2]


        if fit_function == "amparo":
            fitfunc_ell_heights = lambda p, x: func_ellipse_heights_amparo(x, p[0], p[1], p[2], p[3], p[4])
            print("Using function definition: AMPARO")
        elif fit_function == "xianbo":
            fitfunc_ell_heights = lambda p, x: func_ellipse_heights_xianbo(x, p[0], p[1], p[2], p[3], p[4])
            print("Using function definition: XIANBO")
        else:
            fitfunc_ell_heights = lambda p, x: func_ellipse_heights_dabam(x, p[0], p[1], p[2], p[3], p[4])
            print("Using function definition: DABAM")

        print("======== Fitting heights =======")
        errfunc_ell_heights = lambda p, x, y: fitfunc_ell_heights(p, x) - y

        p_guess = [p0,q0,theta0,hshift,vshift]

        if ibounded == 0:
            print("======== Curve fitting without guess =======")
            if fit_function == "amparo":
                popt, cov_x = curve_fit(func_ellipse_heights_amparo, y0, z0, maxfev=10000)
            elif fit_function == "xianbo":
                popt, cov_x = curve_fit(func_ellipse_heights_xianbo, y0, z0, maxfev=10000)
            else:
                popt, cov_x = curve_fit(func_ellipse_heights_dabam, y0, z0, maxfev=10000)
        elif ibounded == 1:
            print("======== Least Squares fitting =======")
            popt, cov_x, infodic, mesg, ier = leastsq(errfunc_ell_heights, p_guess,args=(y0, z0), full_output=True)
        else:
            print("======== Bounded Least Squares fitting =======")
            # https://github.com/jjhelmus/leastsqbound-scipy
            from leastsqbound import leastsqbound
            bounds = [ (p0*0.998,p0*1.002) , (q0*0.8,q0*1.2), (theta0*0.98,theta0*1.02), (-0.01,0.01), (-0.001,0.001)]
            popt, cov_x, infodic, mesg, ier = leastsqbound(errfunc_ell_heights, p_guess,
                                                           args=(y0,z0), full_output=True, bounds=bounds)
        print(">>>>> p_guess:", p_guess)
        print(">>>>> popt:", popt)

        z1 = fitfunc_ell_heights(p_guess, y0)
        z2 = fitfunc_ell_heights(popt, y0)

        zp1 = numpy.gradient(z1,y0[1]-y0[0])
        zp2 = numpy.gradient(z2,y0[1]-y0[0])


        rms_values = ( 1e9*(z0-z1).std(),1e9*(z0-z2).std(),
                       1e6*numpy.gradient(z0-z1,y0[1]-y0[0]).std(),1e6*numpy.gradient(z0-z2,y0[1]-y0[0]).std())

        print ('Height error RMS z0-z1:             %.3f nm'%(1e9*(z0-z1).std()))
        print ('Slope error RMS z0-z1:             %.3f urad'%(1e6*numpy.gradient(z0-z1,y0[1]-y0[0]).std()))

        print ('height error RMS z0-z2:             %.3f nm'%(1e9*(z0-z2).std()))
        print ('Slope error RMS z0-z2:             %.3f urad'%(1e6*numpy.gradient(z0-z2,y0[1]-y0[0]).std()))

        #dump file
        outFile = "fit.spec"
        f = open(outFile,'w')
        f.write("#F %s\n"%outFile)
        f.write("\n#S 1 heights profiles\n")
        f.write("#N 5\n")
        f.write("#L coordinate [mm]  heights_orig [um]  ellipse_guess [um]  ellipse_fit [um]  heights-fit [nm]\n")
        for i in range(y0.size):
            f.write("%f %f %f %f %f \n"%(y0[i]*1e3,z0[i]*1e6,z1[i]*1e6,z2[i]*1e6,(z0[i]-z2[i])*1e9))
        f.close()


    if fit_method == 1: # ellipse, slopes fit

        #preprocessors ellipse fitting slopes
        if entry_number == 4:
            pass
        elif entry_number == 6:
            pass
        elif entry_number == 19:
            pass
        elif entry_number == 20:
            pass
        elif entry_number == 21:
            pass

        if fit_function == "amparo":
            fitfunc_ell_slopes  =  lambda p, x: func_ellipse_slopes_amparo(x, p[0], p[1], p[2], p[3])
        elif fit_function == "dabam":
            fitfunc_ell_slopes  =  lambda p, x: func_ellipse_slopes_dabam(x, p[0], p[1], p[2], p[3])
        else:
            print("Not implemented function: ",fit_function)
            raise NotImplementedError

        print("======== Fitting slopes =======")

        errfunc_ell_slopes = lambda p, x, y: fitfunc_ell_slopes(p, x) - y
        p_guess = [p0,q0,theta0,0.0]

        zp1 = fitfunc_ell_slopes(p_guess, y0)

        if ibounded == 0:
            print("======== Curve fitting without guess =======")

            if fit_function == "amparo":
                popt, cov_x = curve_fit(func_ellipse_slopes_amparo, y0, zp0, maxfev=10000)
            elif fit_function == "dabam":
                popt, cov_x = curve_fit(func_ellipse_slopes_dabam, y0, zp0, maxfev=10000)
            else:
                print("Not implemented function: ",fit_function)
                raise NotImplementedError

        elif ibounded == 1:
                print("======== Least Squares fitting =======")
                popt, cov_x, infodic, mesg, ier = leastsq(errfunc_ell_slopes, p_guess, args=(y0, zp0),
                                                          full_output=True)
        else:
                print("======== Bounded Least Squares fitting =======")
                # https://github.com/jjhelmus/leastsqbound-scipy
                from leastsqbound import leastsqbound
                bounds = [ (p0*0.9,p0*1.1) , (q0*0.9,q0*1.1), (theta0*0.9,theta0*1.1), (-2.0,2.0)]
                popt, cov_x, infodic, mesg, ier = leastsqbound(errfunc_ell_slopes, p_guess, args=(y0,zp0),
                                                            bounds=bounds,full_output=True)

        print(">>>>> p_guess:", p_guess)
        print(">>>>> popt (p,q,theta): ",popt)

        zp2 = fitfunc_ell_slopes(popt, y0)

        z1 = cdf(y0,zp1)
        z2 = cdf(y0,zp2)

        rms_values = (1e9*(cdf(y0,zp0-zp1)).std(),1e9*(cdf(y0,zp0-zp2)).std(),
                      1e6*(zp0-zp1).std(),1e6*(zp0-zp2).std() )

        print ('Height error RMS z0-z1:            %.3f nm'%(1e9*(cdf(y0,zp0-zp1)).std()))
        print ('Slope error RMS zp0-zp1:             %.3f urad'%(1e6*(zp0-zp1).std()))

        print ('Height error RMS z0-z2:            %.3f nm'%(1e9*(cdf(y0,zp0-zp2)).std()))
        print ('Slope error RMS zp0-zp2:             %.3f urad'%(1e6*(zp0-zp2).std()))

        #dump file
        outFile = "fit.spec"
        f = open(outFile,'w')
        f.write("#F %s\n"%outFile)
        f.write("\n#S 1 slopes profiles\n")
        f.write("#N 5\n")
        f.write("#L coordinate [mm]  slopes_orig [urad]  ellipse_guess [urad]  ellipse_fit [urad]  slopes-fit [nrad]\n")
        for i in range(y0.size):
            f.write("%f %f %f %f %f \n"%(y0[i]*1e3,zp0[i]*1e6,zp1[i]*1e6,zp2[i]*1e6,(zp0[i]-zp2[i])*1e9))
        f.close()

        print(">>>>> popt (p,q,theta): ",popt)


        rms_values = ( 1e9*(z0-z1).std(),1e9*(z0-z2).std(),
           1e6*numpy.gradient(z0-z1,y0[1]-y0[0]).std(),1e6*numpy.gradient(z0-z2,y0[1]-y0[0]).std())

        if ibounded != 0: print ('Height error RMS z0-z1:             %.3f nm'%(1e9*(z0-z1).std()))
        print ('height error RMS z0-z2:             %.3f nm'%(1e9*(z0-z2).std()))
        dd=numpy.concatenate( (y0, z2) ,axis=0).reshape(2,-1).transpose()
        outFile = "tmp_z2.dat"
        numpy.savetxt(outFile,dd)
        print ("File "+outFile+" written to disk:\n")

    #
    #plots
    #
    if do_plots:
        #
        #plots
        #
        from matplotlib import pylab as plt


        if (fit_method == 0): #heights
            f1 = plt.figure(1)
            ax = plt.gca()
            ax.get_xaxis().get_major_formatter().set_useOffset(False)
            ax.get_yaxis().get_major_formatter().set_useOffset(False)
            plt.plot(y0*1e3,z0*1e6)
            plt.plot(y0*1e3,z1*1e6)
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
        elif fit_method == 1: #slopes
            f1 = plt.figure(1)
            ax = plt.gca()
            ax.get_xaxis().get_major_formatter().set_useOffset(False)
            ax.get_yaxis().get_major_formatter().set_useOffset(False)
            plt.plot(y0*1e3,zp0*1e6)
            if ibounded != 0: plt.plot(y0*1e3,zp1*1e6)
            plt.plot(y0*1e3,zp2*1e6)
            plt.title("slopes data (blue), starting (green) and optimized (red) ellipse")
            plt.xlabel("Y [mm]")
            plt.ylabel("Zp [urad]")

            f2 = plt.figure(2)
            ax = plt.gca()
            ax.get_xaxis().get_major_formatter().set_useOffset(False)
            ax.get_yaxis().get_major_formatter().set_useOffset(False)
            plt.plot(y0*1e3,(zp0-zp2)*1e6)
            plt.title("residual slopes")
            plt.xlabel("Y [mm]")
            plt.ylabel("Zp [urad]")


        plt.show()


    return rms_values
#
# main program
#
if __name__ == '__main__':

    loop = 1 # =0=single run, 1=loop

    fit_function = "dabam"   # options dabam amparo xianbo
    ibounded = 0             # 0=no guess, 1=with guess, 2=bounded guesses. ibounded=2 not working


    if loop:
        entry_number = [4,6,19,20,21]
        out = numpy.zeros( (6,2*len(entry_number))) + 22

        kk = -1
        for i,entry in enumerate(entry_number):
            for j in range(2):
                kk += 1
                out[0,kk] = entry
                out[1,kk] = j

                tmp = ellipse_fit(entry_number=entry,fit_method=j,ibounded=ibounded,do_plots=0)
                print(tmp)

                out[2,kk] = tmp[0]
                out[3,kk] = tmp[1]
                out[4,kk] = tmp[2]
                out[5,kk] = tmp[3]


        print("\n\nfit_function=%s, ibounded=%d"%(fit_function,ibounded))
        print("entry fit:hgt/slp  heights_error slopes_error")
        for i in range(2*len(entry_number)):
            print("%3d       %1d         %5.2f nm    %5.2f urad"%(out[0,i],out[1,i],out[3,i],out[5,i]))
    else:
        fit_method = 0           # 0=fit heights, 1=fit slopes
        entry = 6 # [4,6,19,20,21]
        tmp = ellipse_fit(entry_number=entry,fit_method=fit_method,fit_function=fit_function,ibounded=ibounded,do_plots=1)


