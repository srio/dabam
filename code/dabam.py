"""

dabam: (dataBase for metrology)
       python tools for processing remote files containing the results
       of metrology measurements on X-ray mirrors

       functions: 
             cdf (calculate cumulative distribution function)
             psd (calculate power spectral density)
             write_shadowSurface (writes file with a mesh for SHADOW)
             func_ellipse_slopes
             func_ellipse
             func_ellipse_slopes_amparo
             func_ellipse_amparo
             write_shadowSurface
             get_arguments
             get_metadata_and_data

 
       MODIFICATION HISTORY:
           20130902 srio@esrf.eu, written
           20131109 srio@esrf.eu, added command line arguments, access metadata

"""

__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2013-2015"


import numpy



#
#testing OO way
#
class dabam(object):
    def __init__(self,dictionnary={}):
        self.inputs = dictionnary
    def info(self):
        argsdict = self.inputs
        print ("-----------------------------------------------------")
        for i,j in argsdict.items():
            print ("dabam.info(): %s = %s" % (i,j))
        print ("-----------------------------------------------------")


#
# functions 
# 

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


def psd(x, y, onlyrange = None):
    """
     psd: Calculates the PSD (power spectral density) from a profile

      INPUTS:
           x - 1D array of (equally-spaced) lengths.
           y - 1D array of heights.
      OUTPUTS:
           f - 1D array of spatial frequencies, in units of 1/[x].
           s - 1D array of PSD values, in units of [y]^3.
      KEYWORD PARAMETERS:
           onlyrange - 2-element array specifying the min and max spatial
               frequencies to be considered. Default is from
               1/(length) to 1/(2*interval) (i.e., the Nyquist
               frequency), where length is the length of the scan,
               and interval is the spacing between points.
     
      PROCEDURE
            S=Length*ABS(ifft(Y*Window)^2
            Where Length is as described above, and Window is the value of
            the optional window function 
     
    """
    n_pts = x.size
    if (n_pts <= 1): 
        print ("psd: Error, must have at least 2 points.")
        return 0

    xx = x
    yy = y

    window=yy*0+1.
    length=xx.max()-xx.min()  # total scan length.
    # psd 
    s=length*numpy.absolute(numpy.fft.ifft(yy*window)**2) 

    s=s[1:(n_pts/2+1*numpy.mod(n_pts,2))]  # take an odd number of points.
    n_ps=s.size                       # number of psd points.
    interval=length/(n_pts-1)         # sampling interval.
    f_min=1./length                   # minimum spatial frequency.
    f_max=1./(2.*interval)            # maximum (Nyquist) spatial frequency.
    # spatial frequencies.
    f=numpy.arange(float(n_ps))/(n_ps-1)*(f_max-f_min)+f_min 

    if onlyrange != None:
        roi =  (f <= onlyrange[1]) * (f >= onlyrange[0]) 
        if roi.sum() > 0:
            roi = roi.nonzero()
            f = f[roi]
            s = s[roi]

    return s,f

def func_ellipse_slopes(x, p, q, theta, shift):


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

    # print(">>>> func_ellipse_slopes: a=%f, b=%f, c=%f"%(a,b,c))
    # print(">>>> func_ellipse_slopes: x0=%f, y0=%f"%(x0,y0))

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

    return ells+shift

def write_shadowSurface(s,xx,yy,outFile='presurface.dat'):
    """
      write_shadowSurface: writes a mesh in the SHADOW/presurface format

      SYNTAX: 
           out = write_shadowSurface(z,x,y,outFile=outFile)

      INPUTS:
           z - 2D array of heights
           x - 1D array of spatial coordinates along mirror width.
           y - 1D array of spatial coordinates along mirror length.
     
      OUTPUTS:
           out - 1=Success, 0=Failure
           outFile - output file in SHADOW format. If undefined, the
                     file is names "presurface.dat"
     
    """
    out = 1

    try:
       fs = open(outFile, 'w')
    except IOError:
       out = 0
       print ("Error: can\'t open file: "+outFile)
       return 
    else:
        # dimensions
        fs.write( "%d  %d \n"%(xx.size,yy.size))
        # y array
        for i in range(yy.size): 
            fs.write(' ' + repr(yy[i]) )
        fs.write("\n")
        # for each x element, the x value and the corresponding z(y) profile
        for i in range(xx.size): 
            tmps = ""
            for j in range(yy.size): 
                tmps = tmps + "  " + repr(s[j,i])
            fs.write(' ' + repr(xx[i]) + " " + tmps )
            fs.write("\n")
        fs.close()
        print ("File for SHADOW "+outFile+" written to disk.")

def get_arguments():
    import argparse # to manage input parameters from command-line argument

    #
    # define default aparameters taken from command arguments
    #
    parser = argparse.ArgumentParser(description="dabam.py: python program to access and evaluate DAta BAse for Metrology (DABAM) files. See http://ftp.esrf.eu/pub/scisoft/dabam/Readme.md")
    # main argument
    parser.add_argument('entryNumber', nargs='?', metavar='N', type=int, default=0,
        help='an integer indicating the DABAM entry number')
    # options
    parser.add_argument('-v', '--verbose', action='store_true',
        help='print some debugging messages. Default is No')

    parser.add_argument('-l', '--localFileRoot',
        help='Define the name of local DABAM file root '+
        '(<name>.dat for data, <name>.txt for metadata).'+
        ' If undefined use remote file')
    parser.add_argument('-r', '--rootFile', default='tmp',
        help='Define the root for output files. Default is "tmp".')

    parser.add_argument('-D', '--setDetrending', default=-2,
        help='Detrending: if >0 is the polynomial degree, -1=skip, -2=automatic, -3=ellipse. Default=-2')

    parser.add_argument('-B', '--calcBoundaries', action='store_true',
        help='if set, calculate specular-diffuse scattering boundary.')
    parser.add_argument('-H', '--histoCalc', action='store_true',
        help='if set, calculate histograms.')
    parser.add_argument('-b', '--histoB', default=1e-7,
        help='If histogram is calculated, this is the binsize in rads. Default is 1e-7')
    parser.add_argument('-s', '--shadowCalc', action='store_true',
        help='if set, write file with mesh for SHADOW')
    parser.add_argument('--shadowNy', default=199,
        help='For SHADOW file, the number of points along Y (length). Default=199')
    parser.add_argument('--shadowNx', default=11,
        help='For SHADOW file, the number of points along X (width). Default=11')
    parser.add_argument('--shadowWidth', default=6.0,
        help='For SHADOW file, the surface dimension along X (width) in cm. Default=6.0')

    parser.add_argument('-m','--multiply', default=1.0,
        help='Multiply input (slope) by this number (to play with RMS values). Default=1.0')

    parser.add_argument('-S', '--useHeightsOrSlopes', default=-1,
        help='using profile heights (0) or slopes (1). Overwrites FILE_FORMAT keyword. Default=1')

    parser.add_argument('-A', '--useAbscissasColumn', default=0,
        help='using abscissas column index. Default=0')

    parser.add_argument('-O', '--useOrdinatesColumn', default=1,
        help='using ordinates column index. Default=1')

    parser.add_argument('-P', '--plot', default="",
        help='plot: heights slopes psd_h psd_s cdf_h cdf_s')

    args = parser.parse_args()

    return args

def get_metadata_and_data(args):
    from io import StringIO
    #import urllib2  # for remote access

    try:
        # For Python 3.0 and later
        from urllib.request import urlopen
    except ImportError:
        # Fall back to Python 2's urllib2
        from urllib2 import urlopen


    import json # for decoding metadata


    #;
    #; inputs
    #;
    input_option = args.entryNumber

    if (args.localFileRoot == None):
        remoteAccess = 1  # 0=Local file, 1=Remote file
    else:
        remoteAccess = 0  # 0=Local file, 1=Remote file
        inFileRoot = args.localFileRoot

    if input_option == 0:
        remoteAccess = 0
        inFileRoot = args.localFileRoot
    else:
        inFileRoot = 'dabam-'+str(input_option)#;
    #;
    #; load file with slopes
    #;
    inFileDat = inFileRoot+'.dat'
    inFileTxt = inFileRoot+'.txt'
    if remoteAccess:
        myServer = 'http://ftp.esrf.eu/pub/scisoft/dabam/data/'
        # metadata
        myfileurl = myServer+inFileTxt
        myfileurl = myfileurl

        try:
            u = urlopen(myfileurl)
        except:
            print ("Error accessing remote file: "+myfileurl+" does not exist.")
            sys.exit()


        ur = u.read()
        ur1 = ur.decode(encoding='UTF-8')
        h = json.loads(ur1) # dictionnary with metadata

        # data
        myfileurl = myServer+inFileDat
        try:
            u = urlopen(myfileurl)
        except:
            print ("Error accessing remote file: "+myfileurl+" does not exist.")
            sys.exit()

        print()

        ur = u.read()
        #ur1 = ur.decode(encoding='UTF-8')
        skipLines = h['FILE_HEADER_LINES']

        import sys
        if sys.version_info[0] == 2:
            ur = StringIO( unicode(ur) )
        else:
            #print("+++++++++++++",ur.decode())
            #ur = StringIO( ur.decode(encoding='UTF-8') )
            ur = StringIO( ur.decode(encoding='ASCII') )

        a = numpy.loadtxt(ur, skiprows=skipLines )
        #print("+++++++++++++",a.shape)
        #a.shape = (-1,2)
        #print("+++++++++++++",a.shape)

    else:
        with open(inFileTxt, mode='r') as f1:
            h = json.load(f1)
        skipLines = h['FILE_HEADER_LINES']
        a = numpy.loadtxt(inFileDat, skiprows=skipLines) #, dtype="float64" )
        myServer = None

    return h,a,inFileTxt,inFileDat,myServer

def main():
    import json # for decoding metadata

    #
    # get arguments of dabam call
    #
    args = get_arguments()

    if (args.localFileRoot == None):
        if args.entryNumber == 0:
            raise Exception("Usage: python dabam.py <entry_number>")

    #
    #  list arguments
    #

    argsdict = vars(args)
    if (args.verbose == True):
        print ("-----------------------------------------------------")
        for i,j in argsdict.items():
            print ("%s = %s" % (i,j))
        print ("-----------------------------------------------------")


    dm = dabam(argsdict)

    # dm.info()

    #
    #retrieve data and metadata
    #
    h,a,inFileTxt,inFileDat,myServer = get_metadata_and_data(args)

    #;
    #; convert to SI units (m,rad)
    #;
    a[:,0] = a[:,0]*h['X1_FACTOR']
    a[:,1] = a[:,1]*h['Y1_FACTOR']
    ncols = a.shape[1]

    if int(h["FILE_FORMAT"]) <= 2:
        for i in range(2,ncols):    # X1 Y1 Y2 Y3...
            a[:,i] = a[:,i]*h['Y%d_FACTOR'%i]
    elif int(h["FILE_FORMAT"]) == 3: #X1 Y1 X2 Y2 etc
        ngroups = int(ncols / 2)
        icol = 1
        for i in range(2,ngroups):    # X1 Y1 Y2 Y3...
            icol += 1
            a[:,icol] = a[:,icol]*h['X%d_FACTOR'%i]
            icol += 1
            a[:,icol] = a[:,icol]*h['Y%d_FACTOR'%i]

    #
    #; apply multiplicative factor
    #
    if (args.multiply != 1.0):
        a[:,1] = a[:,1]  * float(args.multiply)

    #
    # select columns with abscissas and ordinates
    #
    col_abscissas = int(args.useAbscissasColumn)
    col_ordinates = int(args.useOrdinatesColumn)


    col_ordinates_title = 'unknown'
    if int(args.useHeightsOrSlopes) == -1:  #default, read from file
        if h['FILE_FORMAT'] == 1:  # slopes in Col2
            col_ordinates_title = 'slopes'
        if h['FILE_FORMAT'] == 2:  # heights in Col2
            col_ordinates_title = 'heights'
        if h['FILE_FORMAT'] == 3:  # slopes in Col2, file X1 Y1 X2 Y2
            col_ordinates_title = 'slopes'
    else:
        if int(args.useHeightsOrSlopes) == 0:
            col_ordinates_title = 'heights'
        if int(args.useHeightsOrSlopes) == 1:
            col_ordinates_title = 'slopes'

    print("Using: abscissas column index %d (mirror coordinates)"%(col_abscissas))
    print("       ordinates column index %d (profile %s)"%(col_ordinates,col_ordinates_title))


    #;
    #; Detrending:
    #; substract linear fit to the slopes (remove best circle from profile)
    #;
    if col_ordinates_title == 'slopes':
        sy = a[:,col_abscissas]
        sz = a[:,col_ordinates]
    elif col_ordinates_title == 'heights':
        sy = a[:,col_abscissas]
        #TODO we suppose that data are equally spaced. Think how to generalise
        sz = numpy.gradient(a[:,col_ordinates],(sy[1]-sy[0]))
    else:
        raise NotImplementedError

    sz1 = numpy.copy(sz)
    polDegree = int(args.setDetrending)

    # define detrending to apply: >0 polynomial prder, -1=None, -2=Default, -3=elliptical
    if polDegree == -2: # this is the default
        if (h['SURFACE_SHAPE']).lower() == "elliptical":
            polDegree = -3     # elliptical detrending
        else:
            polDegree = 1      # linear detrending


    if polDegree >= 0: # polinomial fit
        coeffs = numpy.polyfit(sy, sz, polDegree)
        pol = numpy.poly1d(coeffs)
        zfit = pol(sy)
        sz = sz - zfit

    if polDegree == -3: # ellipse
        try:
            from scipy.optimize import curve_fit
        except:
            raise ImportError("Cannot perform ellipse detrending: please install scipy")

        popt, cov_x = curve_fit(func_ellipse_slopes, sy, sz, maxfev=10000)
        zfit= func_ellipse_slopes(sy, popt[0], popt[1], popt[2], popt[3])
        sz = sz - zfit

    
    
    #;
    #; histogram (optional) (for fun, not needed)
    #;
    
    if (args.histoCalc == True):
        binsize = args.histoB # default is 1e-7 rads
        bins = numpy.ceil( (sz.max()-sz.min())/binsize )
        hz,hy = numpy.histogram(sz, bins = bins)
        
        #calculate positions of the center of the bins
        hyc = hy+0.5*(hy[1]-hy[0])
        hyc = hyc[0:-1]
        
        #Gaussian
        g = numpy.exp( -numpy.power(hyc-sz.mean(),2)/2/numpy.power(sz.std(),2) )
        g = g/g.sum()*hz.sum()
        
        # dump histogram to file (3 cols: center of bin, counts, Gaussian value)
        dd=numpy.concatenate( (hyc, hz, g ) ,axis=0).reshape(3,-1).transpose()
        outFile = args.rootFile+'Histo.dat'
        numpy.savetxt(outFile,dd)
        print ("File "+outFile+" written to disk. Cols: slope[rad], counts.")
    
    #;
    #; calculate heights by integrating the slope
    #;
    zprof = cdf(sy,sz)
    zprof1 = cdf(sy,sz1)

    #
    # Dump heights and slopes profiles to files
    #
    if args.rootFile != "":
        outFile = args.rootFile+'Heights.dat'
        dd=numpy.concatenate( (sy.reshape(-1,1), zprof.reshape(-1,1)),axis=1)
        numpy.savetxt(outFile,dd)
        print ("File "+outFile+" written to disk: Columns are:\n  coordinate(m),height(m).")

        outFile = args.rootFile+'Slopes.dat'
        dd=numpy.concatenate( (sy.reshape(-1,1), sz.reshape(-1,1)),axis=1)
        numpy.savetxt(outFile,dd)
        print ("File "+outFile+" written to disk: Columns are:\n  coordinate(m),slopes(rad).")
    
    #;
    #; calculate PSD on both profile and slope, and also then their cdf()
    #;
    psdHeights,f = psd(sy,zprof,onlyrange=None)
    psdSlopes,f = psd(sy,sz,onlyrange=None)
    cdfHeights = numpy.sqrt(cdf(f,psdHeights))
    cdfHeightsRMS = cdfHeights.max()
    cdfHeights = 1.0-cdfHeights/cdfHeightsRMS
    cdfSlopes = numpy.sqrt(cdf(f,psdSlopes)) 
    cdfSlopesRMS = cdfSlopes.max() 
    cdfSlopes = 1.0 - cdfSlopes/cdfSlopesRMS

    if args.rootFile != "":
        dd = numpy.concatenate( (f, psdHeights, psdSlopes, cdfHeights, cdfSlopes ) ,axis=0).reshape(5,-1).transpose()
        outFile = args.rootFile+'PSD.dat'
        numpy.savetxt(outFile,dd)
        print ("File "+outFile+" written to disk:  Columns are:\n"+
               "  freq (m^-1),psd_prof(m^3),psd_slope(rad^3),\n"+
               "              cdf(psd_prof),cdf(psd_slope).")
    
    
    #;
    #;  write file for SHADOW (optional)
    #;  replicate the (x,z) profile in a "s" mesh of npointsx * npointsy
    #; 
    
    if (args.shadowCalc == True):
        #inputs
        npointsy = int(args.shadowNy)
        npointsx = int(args.shadowNx)
        mirror_width = args.shadowWidth # in cm
        
        # units to cm
        y = sy * 100.0 # from m to cm
        z = zprof * 100.0 # from m to cm
        
        # set origin at the center of the mirror. TODO: allow any point for origin
        z = z - z.min()
        y = y - y[int(y.size/2)]


        # interpolate the profile (y,z) to have npointsy points (new profile yy,zz)
        mirror_length = y.max() - y.min()
        yy = numpy.linspace(-mirror_length/2.0,mirror_length/2.0,npointsy)
        zz = numpy.interp(yy,y,z)

        # dump to file interpolated profile (for fun)
        dd = numpy.concatenate( (yy.reshape(-1,1), zz.reshape(-1,1)),axis=1)
        outFile = args.rootFile + "ProfileInterpolated.dat"
        numpy.savetxt(outFile,dd)
        
        # fill the mesh arrays xx,yy,s with interpolated profile yy,zz
        xx=numpy.linspace(-mirror_width/2.0,mirror_width/2.0,npointsx)
        s = numpy.zeros( (npointsy,npointsx) )
        for i in range(npointsx):
            s[:,i]=zz
        
        # write Shadow file
        outFile = args.rootFile + "Shadow.dat"
        tmp = write_shadowSurface(s,xx,yy,outFile=outFile)
        
        # write header file
        outFile = args.rootFile + "Header.txt"
        with open(outFile, mode='w') as f1:
            json.dump(h, f1, indent=2)
        print ("File "+outFile+" containing metadata written to disk.")

    
    #; 
    #; separation between diffuse scattering and specular reflecting domains
    #; TODO: this is experimental, to refine the methodology
    #; 
    if (args.calcBoundaries == True): 
        wavelength  = 1.0e-10 # 0.8e-10 # m 
        theta = 3e-3 # 3.65e-3 # rad
        h_cut = wavelength/(4*numpy.pi*numpy.sin(theta)) #m 
        h_cut_normalized = h_cut/zprof.std()
        # note that numpy.interp requires increasing abscissas
        inew = cdfHeights.argsort()
        f_cut = numpy.interp( h_cut_normalized ,cdfHeights[inew],f[inew])



    #; 
    #; info
    #; 
    #
    print ('\n---------- profile results --------------------')
    if (args.localFileRoot == None):
        print ('Remote directory:\n   %s'%myServer)
    print ('Data File:     %s'%inFileDat)
    print ('Metadata File: %s'%inFileTxt)
    print ('Surface shape: %s '%(h['SURFACE_SHAPE']))
    print ('Facility:      %s'%(h['FACILITY']))
    print ('Scan length: %.3f mm'%(1e3*(sy[-1]-sy[0])))
    print ('Number of points: %d'%(len(sy)))

    print ('   ')
    if polDegree >= 0:
        if polDegree == 1:
            print ("Linear detrending: z'=%g x%+g"%(coeffs[0],coeffs[1]))
            print ('Radius of curvature: %.3F m'%(1.0/coeffs[-2]))
        else:
            print ('Polynomial detrending coefficients: '+repr(coeffs))
    elif polDegree == -1:
       print ('No detrending applied.')
    elif polDegree == -3:
       print ('Ellipse detrending applied.')

    print ('   ')
    print ('   ')
    print ('Slope error s_RMS:             %.3f urad'%(1e6*sz.std()))
    print ('         from PSD:             %.3f urad'%(1e6*cdfSlopesRMS))
    print ('         from USER (metadata): %s urad'%(h['CALC_SLOPE_RMS']))
    print ('Shape error h_RMS:              %.3f nm'%(1e9*zprof.std()))
    print ('         from PSD:              %.3f nm'%(1e9*cdfHeightsRMS))
    print ('         from USER (metadata):  %s nm'%(h['CALC_HEIGHT_RMS']))
    print ('PV of height profile (before detrend): %.3f nm'%(  1e9*(zprof1.max()-zprof1.min() )))
    print ('PV of height profile (after detrend):  %.3f nm'%(  1e9*(zprof.max()-zprof.min() )))
    print ('   ')
    if (args.calcBoundaries == True): 
        print ('Diffraction/reflection boundary zones: ')
        print ('   diffraction if h < h_cut, reflection if h > h_cut: ')
        print ('   at wavelength [A]: ')+repr(wavelength*1e10)
        print ('   at grazing angle [mrad]: ')+repr(theta*1e3)
        print ('   ')
        print ('   height cut h_cut [nm]: '+repr(  h_cut*1e9 ))
        print ('   h_cut/h_RMS: '+repr(  h_cut_normalized )) 
        print ('   frequency cut f_cut [m^-1]: '+repr(  f_cut ))
        print ('   length 1/f_cut [mm]: '+repr(  1.0/f_cut ))

    print ('%d & %d & %.2f (%.2f) & %.2f & %.2f (%.2f) & %.2f \\\\'%(args.entryNumber,
            int(1e3*(sy[-1]-sy[0])),
            1e6*sz.std(),   float(-1 if h['CALC_SLOPE_RMS'] is None else h['CALC_SLOPE_RMS'])  , 1e6*cdfSlopesRMS,
            1e9*zprof.std(),float(-1 if h['CALC_HEIGHT_RMS'] is None else h['CALC_HEIGHT_RMS']), 1e9*cdfHeightsRMS, ))
    print ('-------------------------------------------------\n')
    print (' ')

    if (args.plot != ""):
        try:
            from matplotlib import pylab as plt
        except:
            raise ImportError

        print("plotting: ",args.plot)
        if (args.plot == "heights" or args.plot == "all"):
            f1 = plt.figure(1)
            plt.plot(1e3*sy,1e6*zprof)
            plt.title("heights profile")
            plt.xlabel("Y [mm]")
            plt.ylabel("Z [um]")
        if (args.plot == "slopes" or args.plot == "all"):
            f2 = plt.figure(2)
            plt.plot(1e3*sy,1e6*sz)
            plt.title("slopes profile" or args.plot == "all")
            plt.xlabel("Y [mm]")
            plt.ylabel("Zp [urad]")
        if (args.plot == "psd_h" or args.plot == "all"):
            f3 = plt.figure(3)
            plt.loglog(f,psdHeights)
            plt.title("PSD of heights profile")
            plt.xlabel("f [m^-1]")
            plt.ylabel("PSD [m^3]")
        if (args.plot == "psd_s" or args.plot == "all"):
            f4 = plt.figure(4)
            plt.loglog(f,psdSlopes)
            plt.title("PSD of slopes profile")
            plt.xlabel("f [m^-1]")
            plt.ylabel("PSD [rad^3]")
        if (args.plot == "cdf_h" or args.plot == "all"):
            f5 = plt.figure(5)
            plt.semilogx(f,cdfHeights)
            plt.title("Lambda CDF(PDF) of heights profile")
            plt.xlabel("f [m^-1]")
            plt.ylabel("heights Lambda")
        if (args.plot == "cdf_s" or args.plot == "all"):
            f6 = plt.figure(6)
            plt.semilogx(f,cdfSlopes)
            plt.title("Lambda CDF(PDF) of slopes profile")
            plt.xlabel("f [m^-1]")
            plt.ylabel("slopes Lambda")
        plt.show()
#
# main program
#
if __name__ == '__main__':
    main()