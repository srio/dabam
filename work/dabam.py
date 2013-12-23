"""

dabam: (dataBase for metrology)
       python tools for processing remote files containing the results
       of metrology measurements on X-ray mirrors

       functions: 
             cdf (calculate cumulative distribution function)
             psd (calculate power spectral density)
             write_shadowSurface (writes file with a mesh for SHADOW)
 
       MODIFICATION HISTORY:
           20130902 srio@esrf.eu, written
           20131109 srio@esrf.eu, added command line arguments, access metadata
           20131223 srio@esrf.eu, added multi-column support

"""

__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2013"


import numpy

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
      OUTPUTS:
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
       fs = open(outFile, 'wb')
    except IOError:
       out = 0
       print ("Error: can\'t open file: "+outFile)
       return 
    else:
        # dimensions
        fs.write( repr(xx.size)+" "+repr(yy.size)+" \n" ) 
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


#
# main program
# 
if __name__ == '__main__':
    
    import argparse # to manage input parameters from command-line argument
    import StringIO, urllib2  # for remote access
    import json # for decoding metadata

    #
    # define default aparameters taken from command arguments
    #
    parser = argparse.ArgumentParser(description="dabam.py: python program to access and evaluate DAta BAse for Metrology (DABAM) files. See http://ftp.esrf.eu/pub/scisoft/DabamFiles/readme.txt")
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

    parser.add_argument('-D', '--polDegree', default=1, 
        help='degree of plynomial for detrending (<0 to skip). Default=1')

    parser.add_argument('-X', '--Xcol', default=1, 
        help='Index of abscissas column. Default=1')

    parser.add_argument('-Y', '--Ycol', default=-1, 
        help='Index of abscissas column. Default=-1 (last one)')

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


    args = parser.parse_args()
    argsdict = vars(args)


    #
    # list all keywords
    #
    #print 'entryNumber is: ',argsdict['entryNumber']," or: ",args.entryNumber
    if (args.verbose == True):
        print "-----------------------------------------------------"
        for i,j in argsdict.items():
            print "%s = %s" % (i,j)
        print "-----------------------------------------------------"

    #; 
    #; inputs
    #; 
    input_option = args.entryNumber  

    remoteAccess = 1  # 0=Local file, 1=Remote file

    if args.localFileRoot != None:
        remoteAccess = 0 
        inFileRoot = args.localFileRoot 

    #if input_option == 0:
    #    remoteAccess = 0 
    #    inFileRoot = args.localFileRoot 

    #;
    #; load file with slopes 
    #;
    if remoteAccess:  
        inFileRoot = 'dabam-'+str(input_option)
        inFileDat = inFileRoot+'.dat'
        inFileTxt = inFileRoot+'.txt'
        myServer = 'http://ftp.esrf.eu/pub/scisoft/DabamFiles/'
        print ("Accessing remote file: "+inFileDat)
        # metadata
        inFileTxtLong = myServer+inFileTxt
        try:
            u = urllib2.urlopen(inFileTxtLong)
        except:
            print ("Error accesing remote file: "+inFileTxtLong+" does not exist.")
            sys.exit()

        h = json.load(u) # dictionnary with metadata

        # data 
        inFileDatLong = myServer+inFileDat
        try:
            u = urllib2.urlopen(inFileDatLong)
        except:
            print ("Error accesing remote file: "+inFileDatLong+" does not exist.")
            sys.exit()

        s = StringIO.StringIO( u.read() )
        skipLines = h['FILE_HEADER_LINES']
        a = numpy.loadtxt(s, skiprows=skipLines )

    else:
        inFileTxtLong = inFileRoot+".txt"
        inFileDatLong = inFileRoot+".dat"
        print ("Accessing remote file: "+inFileTxtLong)
        with open(inFileTxtLong, mode='r') as f1: 
            h = json.load(f1)
        skipLines = h['FILE_HEADER_LINES']
        a = numpy.loadtxt(inFileDatLong, skiprows=skipLines) #, dtype="float64" )

    #; define the column-index with abscissas and ordinates.
    xcol = int(args.Xcol)-1
    if xcol < 0:
       xcol = a.shape[1] + xcol
    ycol = int(args.Ycol)
    if ycol < 0:
       ycol = a.shape[1] + ycol
    else:
       ycol = ycol -1

    #;
    #; convert to SI units (m,rad)
    #;
    a[:,xcol] = a[:,xcol]*h['X1_FACTOR']
    a[:,ycol] = a[:,ycol]*h['Y'+str(ycol)+'_FACTOR']
    #; apply multiplicative factor
    if (args.multiply != 1.0):
        a[:,ycol] = a[:,ycol]  * float(args.multiply)
        
    #;
    #; Detrending:
    #; substract linear fit to the slopes (remove best circle from profile)
    #;
    
    sy = a[:,xcol]
    sz = a[:,ycol]
    sz1 = numpy.copy(sz)  # keep a copy (undetrended)
    if int(args.polDegree) >= 0:
        coeffs = numpy.polyfit(sy, sz, args.polDegree)
        pol = numpy.poly1d(coeffs)
        zfit = pol(sy)
        sz = sz - zfit
        a[:,ycol] = sz
    
    
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
    #; calculate profile by integrating the slope
    #;
    zprof = cdf(sy,sz)
    zprof1 = cdf(sy,sz1)
    dd=numpy.concatenate( (sy.reshape(-1,1), zprof.reshape(-1,1)),axis=1)
    outFile = args.rootFile+'Profile.dat'
    numpy.savetxt(outFile,dd)
    print ("File "+outFile+" written to disk. Cols: coordinate (m), height (m).")
    
    #;
    #; calculate PSD on both profile and slope, and also then their cdf()
    #;
    psdProfile,f = psd(sy,zprof,onlyrange=None)
    psdSlopes,f = psd(sy,sz,onlyrange=None)
    cdfProfile = numpy.sqrt(cdf(f,psdProfile))
    cdfProfileRMS = cdfProfile.max()
    cdfProfile = 1.0-cdfProfile/cdfProfileRMS
    cdfSlopes = numpy.sqrt(cdf(f,psdSlopes)) 
    cdfSlopesRMS = cdfSlopes.max() 
    cdfSlopes = 1.0 - cdfSlopes/cdfSlopesRMS

    dd=numpy.concatenate( (f, psdProfile, psdSlopes, cdfProfile, cdfSlopes ) ,axis=0).reshape(5,-1).transpose()

    outFile = args.rootFile+'PSD.dat'
    numpy.savetxt(outFile,dd)
    print ("File "+outFile+" written to disk. Cols: freq (m^-1), psd_prof (m^3), psd_slope (rad^3), cdf(psd_prof), cdf(psd_slope).")
    
    
    #;
    #;  write file for SHADOW (optional)
    #;  replicate the (x,z) profile in a "s" mesh of npointsx * npointsy
    #; 
    
    if (args.shadowCalc == True):
        #inputs
        npointsy = args.shadowNy
        npointsx = args.shadowNx
        mirror_width = args.shadowWidth # in cm
        
        # units to cm
        y=sy*10.0 # from m to cm
        z=zprof*10.0 # from m to cm
        
        # set origing at the center of the mirror. TODO: allow any point for origin
        z = z - z.min()
        mirror_length = (y.max()-y.min())
        y = y - 0.5*mirror_length
        
        # interpolate the profile (y,z) to have npointsy points (new profile yy,zz)
        yy=numpy.linspace(-mirror_length/2.0,mirror_length/2.0,npointsy)
        zz=numpy.interp(yy,y,z)
        # dump to file interpolated profile (for fun)
        dd=numpy.concatenate( (yy.reshape(-1,1), zz.reshape(-1,1)),axis=1)
        outFile = args.rootFile + "ProfileInterpolated.dat"
        numpy.savetxt(outFile,dd)
        print ("File "+outFile+" used to prepare SHADOW file written to disk. Cols: coordinate (cm), height (cm).")
        
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
        inew = cdfProfile.argsort()
        f_cut = numpy.interp( h_cut_normalized ,cdfProfile[inew],f[inew])



    #; 
    #; info
    #; 
    #
    print ('\n---------- profile results --------------------')
    print ('Data File: %s, Cols: %d,%d'%(inFileDatLong,xcol+1,ycol+1))
    print ('Metadata File: %s'%inFileTxtLong)
    if int(args.polDegree) < 0:
       print ('No detrending applied')
    elif int(args.polDegree) == 1:
       print ('Linear fit coefficients: '+repr(coeffs))
       print ('Radius of curvature [m] : '+repr(1.0/coeffs[-2]))
    else:
       print ('Polynomial fit coefficients: '+repr(coeffs))

    print ('   ')
    print ('Slope error s_RMS [urad]:           '+repr(1e6*sz.std()))
    print ('                   from PSD [urad]: '+repr(1e6*cdfSlopesRMS))
    print ('Shape error h_RMS [nm]:             '+repr(1e9*zprof.std()))
    print ('            from PSD [nm]:          '+repr(1e9*cdfProfileRMS))
    print ('PV of height profile (before detrend) [nm]: '+repr(  1e9*(zprof1.max()-zprof1.min() )))
    print ('PV of height profile (after detrend) [nm]:  '+repr(  1e9*(zprof.max()-zprof.min() )))
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

    print ('-------------------------------------------------\n')
    print (' ')
