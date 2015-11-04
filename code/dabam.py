"""

dabam: (dataBase for metrology)
       python tools for processing remote files containing the results
       of metrology measurements on X-ray mirrors

       classes:
             dabam
       main functions:
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
           20151103 srio@esrf.eu, restructured to OO

"""

__author__ = "Manuel Sanchez del Rio"
__contact__ = "srio@esrf.eu"
__copyright = "ESRF, 2013-2015"


import numpy
import argparse # to manage input parameters from command-line argument
import sys
import json
import copy
from io import StringIO

try:
    # For Python 3.0 and later
    from urllib.request import urlopen
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen




class dabam(object):
    def __init__(self):
        self.description="dabam.py: python program to access and evaluate DAta BAse for Metrology (DABAM) files. See http://ftp.esrf.eu/pub/scisoft/dabam/Readme.md"

        self.server = "http://ftp.esrf.eu/pub/scisoft/dabam/data/"

        self.inputs = {
            'entryNumber':1,         # 'an integer indicating the DABAM entry number'
            'verbose':True, #False,         # 'print some debugging messages. Default is No'
            'localFileRoot':None,    # 'Define the name of local DABAM file root (<name>.dat for data, <name>.txt for metadata). Default=None, thus use remote file'
            'rootFile':"", # "tmp",        # 'Define the root for output files. Default is "tmp"'
            'setDetrending':-2,      # 'Detrending: if >0 is the polynomial degree, -1=skip, -2=automatic, -3=ellipse. Default=-2'
            'histoCalc':False,       # 'Calculate histograms. Default=No'
            'histoB':1e-7,           # 'If histogram is calculated, this is the binsize in rads. Default is 1e-7'
            'shadowCalc':False,      # 'Write file with mesh for SHADOW. Default=No'
            'shadowNy':199,          # 'For SHADOW file, the number of points along Y (length). Default=199'
            'shadowNx':11,           # 'For SHADOW file, the number of points along X (width). Default=11'
            'shadowWidth':6.0,       # 'For SHADOW file, the surface dimension along X (width) in cm. Default=6.0'
            'multiply':1.0,          # 'Multiply input profile (slope or height) by this number (to play with RMS values). Default=1.0'
            'useHeightsOrSlopes':-1, # 'Force calculations using profile heights (0) or slopes (1). Overwrites FILE_FORMAT keyword. Default=-1 (like FILE_FORMAT)'
            'useAbscissasColumn':0,  # 'Use abscissas column index. Default=0'
            'useOrdinatesColumn':1,  # 'Use ordinates column index. Default=1'
            'plot':None,#
            }
        #to load profiles:
        self.h           = None # metadata
        self.a           = None # raw datafile
        self.sy          = None # abscissa along the mirror
        self.sz1         = None # undetrended slope profile
        self.sz          = None # detrended slope profile
        self.zprof1      = None # undetrended heights profile
        self.zprof       = None # detrended heights profile
        self.coeffs      = None # information on detrending (polynomial coeffs)
        self.f           = None # frequency of Power Spectral Density
        self.psdHeights  = None # Power Spectral Density of Heights profile
        self.psdSlopes   = None # Power Spectral Density of slopes profile
        self.cdfHeights  = None # PDF integral of Heights profile
        self.cdfSlopes   = None # PDF integral of Slopes profile

    #
    #setters
    #

    #variables
    def set_input_entryNumber(self,value):
        self.inputs["entryNumber"] = value
    def set_input_verbose            (self,value):
        self.inputs["verbose"] = value
    def set_input_localFileRoot      (self,value):
        self.inputs["localFileRoot"] = value
    def set_input_rootFile           (self,value):
        self.inputs["rootFile"] = value
    def set_input_setDetrending      (self,value):
        self.inputs["setDetrending"] = value
    def set_input_histoCalc          (self,value):
        self.inputs["histoCalc"] = value
    def set_input_histoB             (self,value):
        self.inputs["histoB"] = value
    def set_input_shadowCalc         (self,value):
        self.inputs["shadowCalc"] = value
    def set_input_shadowNy           (self,value):
        self.inputs["shadowNy"] = value
    def set_input_shadowNx           (self,value):
        self.inputs["shadowNx"] = value
    def set_input_shadowWidth        (self,value):
        self.inputs["shadowWidth"] = value
    def set_input_multiply           (self,value):
        self.inputs["multiply"] = value
    def set_input_useHeightsOrSlopes (self,value):
        self.inputs["useHeightsOrSlopes"] = value
    def set_input_useAbscissasColumn (self,value):
        self.inputs["useAbscissasColumn"] = value
    def set_input_useOrdinatesColumn (self,value):
        self.inputs["useOrdinatesColumn"] = value
    def set_input_plot               (self,value):
        self.inputs["plot"] = value

    #others

    def set_from_command_line(self):
        #
        # define default aparameters taken from command arguments
        #
        parser = argparse.ArgumentParser(description=self.description)
        # main argument
        parser.add_argument('entryNumber', nargs='?', metavar='N', type=int, default=self.get_input_value('entryNumber'),
            help=self.get_input_value_help('entryNumber'))

        # options

        parser.add_argument('-'+self.get_input_value_short_name('verbose'),'--verbose', action='store_true', help=self.get_input_value_help('verbose'))

        parser.add_argument('-'+self.get_input_value_short_name('localFileRoot'), '--localFileRoot', help=self.get_input_value_help('localFileRoot'))

        parser.add_argument('-'+self.get_input_value('rootFile'), '--rootFile', default=self.get_input_value('rootFile'),
            help=self.get_input_value_help('rootFile'))

        parser.add_argument('-'+self.get_input_value_short_name('setDetrending'), '--setDetrending', default=self.get_input_value('setDetrending'),
            help=self.get_input_value_help('setDetrending'))


        parser.add_argument('-'+self.get_input_value_short_name('histoCalc'), '--histoCalc', action='store_true',
            help=self.get_input_value_help('histoCalc'))

        parser.add_argument('-'+self.get_input_value_short_name('histoB'), '--histoB', default=self.get_input_value('histoB'),
            help=self.get_input_value_help('histoB'))


        parser.add_argument('-'+self.get_input_value_short_name('shadowCalc'), '--shadowCalc', action='store_true',
            help=self.get_input_value_help('shadowCalc'))

        parser.add_argument('-'+self.get_input_value_short_name('shadowNy'), '--shadowNy', default=self.get_input_value('shadowNy'),
            help=self.get_input_value_help('shadowNy'))
        parser.add_argument('-'+self.get_input_value_short_name('shadowNx'), '--shadowNx', default=self.get_input_value('shadowNx'),
            help=self.get_input_value_help('shadowNx'))

        parser.add_argument('-'+self.get_input_value_short_name('shadowWidth'), '--shadowWidth', default=self.get_input_value('shadowWidth'),
            help=self.get_input_value_help('shadowWidth'))

        parser.add_argument('-'+self.get_input_value_short_name('multiply'), '--multiply', default=self.get_input_value('multiply'),
            help=self.get_input_value_help('multiply'))


        parser.add_argument('-'+self.get_input_value_short_name('useHeightsOrSlopes'), '--useHeightsOrSlopes', default=self.get_input_value('useHeightsOrSlopes'),
            help=self.get_input_value_help('useHeightsOrSlopes'))

        parser.add_argument('-'+self.get_input_value_short_name('useAbscissasColumn'), '--useAbscissasColumn', default=self.get_input_value('useAbscissasColumn'),
            help=self.get_input_value_help('useAbscissasColumn'))

        parser.add_argument('-'+self.get_input_value_short_name('useOrdinatesColumn'), '--useOrdinatesColumn', default=self.get_input_value('useOrdinatesColumn'),
            help=self.get_input_value_help('useOrdinatesColumn'))

        parser.add_argument('-'+self.get_input_value_short_name('plot'), '--plot', default=self.get_input_value('plot'),
            help=self.get_input_value_help('plot'))


        args = parser.parse_args()

        self.set_input_entryNumber(args.entryNumber)
        self.set_input_verbose(args.verbose)
        self.set_input_localFileRoot(args.localFileRoot)
        self.set_input_rootFile(args.rootFile)
        self.set_input_setDetrending(args.setDetrending)
        self.set_input_histoCalc(args.histoCalc)
        self.set_input_histoB(args.histoB)
        self.set_input_shadowCalc(args.shadowCalc)
        self.set_input_shadowNy(args.shadowNy)
        self.set_input_shadowNx(args.shadowNx)
        self.set_input_shadowWidth(args.shadowWidth)
        self.set_input_multiply(args.multiply)
        self.set_input_useHeightsOrSlopes(args.useHeightsOrSlopes)
        self.set_input_useAbscissasColumn(args.useAbscissasColumn)
        self.set_input_useOrdinatesColumn(args.useHeightsOrSlopes)
        self.set_input_plot(args.plot)


    def set_inputs_from_dictionary(self,dict):
        try:
            self.inputs["entryNumber"]        =  dict["entryNumber"]
            self.inputs["verbose"]            =  dict["verbose"]
            self.inputs["localFileRoot"]      =  dict["localFileRoot"]
            self.inputs["rootFile"]           =  dict["rootFile"]
            self.inputs["setDetrending"]      =  dict["setDetrending"]
            self.inputs["histoCalc"]          =  dict["histoCalc"]
            self.inputs["histoB"]             =  dict["histoB"]
            self.inputs["shadowCalc"]         =  dict["shadowCalc"]
            self.inputs["shadowNy"]           =  dict["shadowNy"]
            self.inputs["shadowNx"]           =  dict["shadowNx"]
            self.inputs["shadowWidth"]        =  dict["shadowWidth"]
            self.inputs["multiply"]           =  dict["multiply"]
            self.inputs["useHeightsOrSlopes"] =  dict["useHeightsOrSlopes"]
            self.inputs["useAbscissasColumn"] =  dict["useAbscissasColumn"]
            self.inputs["useOrdinatesColumn"] =  dict["useOrdinatesColumn"]
            self.inputs["plot"]               =  dict["plot"]
        except:
            raise Exception("Failed setting dabam input parameters from dictionary")

    #
    # tools
    #
    def remote_access(self):
        if (self.get_input_value("localFileRoot") == None):
            remoteAccess = 1  # 0=Local file, 1=Remote file
        else:
            remoteAccess = 0  # 0=Local file, 1=Remote file
        return remoteAccess

    #
    #getters
    #

    def get_input_value(self,key):
        try:
            return self.inputs[key]
        except:
            return None

    def get_inputs_as_dictionary(self):
        return copy.copy(self.inputs)

    def get_input_value_help(self,key):

        if key == 'entryNumber':        return 'an integer indicating the DABAM entry number'
        if key == 'verbose':            return 'print some debugging messages. Default is No'
        if key == 'localFileRoot':      return 'Define the name of local DABAM file root (<name>.dat for data, <name>.txt for metadata). Default=None, thus use remote file'
        if key == 'rootFile':           return 'Define the root for output files. Default is "tmp"'
        if key == 'setDetrending':      return 'Detrending: if >0 is the polynomial degree, -1=skip, -2=automatic, -3=ellipse. Default=-2'
        if key == 'histoCalc':          return 'Calculate histograms. Default=No'
        if key == 'histoB':             return 'If histogram is calculated, this is the binsize in rads. Default is 1e-7'
        if key == 'shadowCalc':         return 'Write file with mesh for SHADOW. Default=No'
        if key == 'shadowNy':           return 'For SHADOW file, the number of points along Y (length). Default=199'
        if key == 'shadowNx':           return 'For SHADOW file, the number of points along X (width). Default=11'
        if key == 'shadowWidth':        return 'For SHADOW file, the surface dimension along X (width) in cm. Default=6.0'
        if key == 'multiply':           return 'Multiply input profile (slope or height) by this number (to play with RMS values). Default=1.0'
        if key == 'useHeightsOrSlopes': return 'Force calculations using profile heights (0) or slopes (1). Overwrites FILE_FORMAT keyword. Default=-1 (like FILE_FORMAT)'
        if key == 'useAbscissasColumn': return 'Use abscissas column index. Default=0'
        if key == 'useOrdinatesColumn': return 'Use ordinates column index. Default=1'
        if key == 'plot':               return 'Plot what? heights slopes psd_h psd_s cdf_h cdf_s. Default=None'

        return ''


    def get_input_value_short_name(self,key):

        if key == 'entryNumber':         return 'N'
        if key == 'verbose':             return 'v'
        if key == 'localFileRoot':       return 'l'
        if key == 'rootFile':            return 'r'
        if key == 'setDetrending':       return 'D'
        if key == 'histoCalc':           return 'H'
        if key == 'histoB':              return 'b'
        if key == 'shadowCalc':          return 's'
        if key == 'shadowNy':            return 'y'
        if key == 'shadowNx':            return 'x'
        if key == 'shadowWidth':         return 'w'
        if key == 'multiply':            return 'm'
        if key == 'useHeightsOrSlopes':  return 'S'
        if key == 'useAbscissasColumn':  return 'A'
        if key == 'useOrdinatesColumn':  return 'O'
        if key == 'plot':                return 'P'

        return '?'


    #
    # file names
    #
    def file_metadata(self):
        return self._file_root()+'.txt'

    def file_data(self):
        return self._file_root()+'.dat'






    #
    # load profile and store data. This is the main action!!
    #

    def load(self):

        # load data and metadata

        self._load_file_metadata()
        self._load_file_data()

        # test consistency
        if (self.get_input_value("localFileRoot") == None):
            if self.get_input_value("entryNumber") == 0:
                raise Exception("Error: entry number cannot be zero for remote access.")

        #calculate detrended profiles
        self._load_profiles()

        #calculate psd
        self._calc_psd()


        #
        # write files
        #
        # write header file
        if self.get_input_value("rootFile") != "":
            outFile = self.get_input_value("rootFile") + "Header.txt"
            with open(outFile, mode='w') as f1:
                json.dump(self.h, f1, indent=2)
            if self.inputs["verbose"]:
                print ("File "+outFile+" containing metadata written to disk.")

        #
        # Dump heights and slopes profiles to files
        #
        if self.get_input_value("rootFile") != "":
            outFile = self.inputs["rootFile"]+'Heights.dat'
            dd=numpy.concatenate( (self.sy.reshape(-1,1), self.zprof.reshape(-1,1)),axis=1)
            numpy.savetxt(outFile,dd)
            if self.inputs["verbose"]:
                print ("File "+outFile+" written to disk: Columns are:\n  coordinate(m),height(m).")

            outFile = self.inputs["rootFile"]+'Slopes.dat'
            dd=numpy.concatenate( (self.sy.reshape(-1,1), self.sz.reshape(-1,1)),axis=1)
            numpy.savetxt(outFile,dd)
            if self.inputs["verbose"]:
                print ("File "+outFile+" written to disk: Columns are:\n  coordinate(m),slopes(rad).")


        #write psd file
        if self.get_input_value("rootFile") != "":
            dd = numpy.concatenate( (self.f, self.psdHeights, self.psdSlopes, self.cdfHeights, self.cdfSlopes ) ,axis=0).reshape(5,-1).transpose()
            outFile = self.get_input_value("rootFile")+'PSD.dat'
            numpy.savetxt(outFile,dd)
            if self.inputs["verbose"]:
                print ("File "+outFile+" written to disk:  Columns are:\n"+
                   "  freq (m^-1),psd_prof(m^3),psd_slope(rad^3),\n"+
                   "              cdf(psd_prof),cdf(psd_slope).")

        #calculate and write histograms
        if self.get_input_value("histoCalc"):
            hyc, hz, hg = self.calc_histograms()
            # dump histogram to file (3 cols: center of bin, counts, Gaussian value)
            if self.get_input_value("rootFile") != "":
                dd=numpy.concatenate( (hyc, hz, hg ) ,axis=0).reshape(3,-1).transpose()
                outFile = self.get_input_value("rootFile")+'Histo.dat'
                numpy.savetxt(outFile,dd)
                print ("File "+outFile+" written to disk. Cols: slope[rad], counts, Gaussian.")

        #shadow file
        if self.get_input_value("shadowCalc"):
            self.write_file_for_shadow()

        if self.inputs["verbose"]:
            self.info_profiles()


    #
    #calculations
    #

    def stdev_profile_heights(self):
        return self.zprof.std()

    def stdev_profile_slopes(self):
        return self.sz.std()

    def stdev_psd_heights(self):
        return ( numpy.sqrt(cdf(self.f,self.psdHeights)) ).max()

    def stdev_psd_slopes(self):
        return ( numpy.sqrt(cdf(self.f,self.psdHeights)) ).max()

    def stdev_user_heights(self):
        return self.h['CALC_HEIGHT_RMS']

    def stdev_user_slopes(self):
        return self.h['CALC_SLOPE_RMS']

    def stdev_summary(self):

        txt = ""
        txt += 'Slope error s_RMS:             %.3f urad\n'       %( 1e6*self.stdev_profile_slopes() )
        txt += '         from PSD:             %.3f urad\n'       %( 1e6*self.stdev_psd_slopes() )
        txt += '         from USER (metadata): %s urad\n'         %(self.stdev_user_slopes())
        txt += 'Shape error h_RMS:              %.3f nm\n'        %(1e9*self.stdev_profile_heights() )
        txt += '         from PSD:              %.3f nm\n'        %(1e9*self.stdev_psd_heights() )
        txt += '         from USER (metadata):  %s nm\n'          %(self.stdev_user_heights() )
        txt += 'PV of height profile (before detrend): %.3f nm\n' %(1e9*(self.zprof1.max() - self.zprof1.min() ))
        txt += 'PV of height profile (after detrend):  %.3f nm\n' %(1e9*(self.zprof.max() - self.zprof.min() ))

        return txt




    def calc_histograms(self):

        sy    = self.sy
        sz1    = self.sz1
        sz    = self.sz
        zprof1    = self.zprof1
        zprof     = self.zprof
        #;
        #; histogram
        #;
        binsize = self.get_input_value("histoB") # default is 1e-7 rads
        bins = numpy.ceil( (sz.max()-sz.min())/binsize )
        hz,hy = numpy.histogram(sz, bins = bins)

        #calculate positions of the center of the bins
        hyc = hy+0.5*(hy[1]-hy[0])
        hyc = hyc[0:-1]

        #Gaussian
        g = numpy.exp( -numpy.power(hyc-sz.mean(),2)/2/numpy.power(sz.std(),2) )
        g = g/g.sum()*hz.sum()


        return hyc, hz, g


    #
    # write things
    #
    def write_file_for_shadow(self):
        #;
        #;  write file for SHADOW (optional)
        #;  replicate the (x,z) profile in a "s" mesh of npointsx * npointsy
        #;

        #inputs
        npointsy = int(self.get_input_value("shadowNy"))
        npointsx = int(self.get_input_value("shadowNx"))
        mirror_width = self.get_input_value("shadowWidth") # in cm

        # units to cm
        y = self.sy * 100.0 # from m to cm
        z = self.zprof * 100.0 # from m to cm

        # set origin at the center of the mirror. TODO: allow any point for origin
        z = z - z.min()
        y = y - y[int(y.size/2)]


        # interpolate the profile (y,z) to have npointsy points (new profile yy,zz)
        mirror_length = y.max() - y.min()
        yy = numpy.linspace(-mirror_length/2.0,mirror_length/2.0,npointsy)
        zz = numpy.interp(yy,y,z)

        # dump to file interpolated profile (for fun)
        dd = numpy.concatenate( (yy.reshape(-1,1), zz.reshape(-1,1)),axis=1)
        outFile = self.get_input_value("rootFile") + "ProfileInterpolatedForShadow.dat"
        numpy.savetxt(outFile,dd)

        # fill the mesh arrays xx,yy,s with interpolated profile yy,zz
        xx=numpy.linspace(-mirror_width/2.0,mirror_width/2.0,npointsx)
        s = numpy.zeros( (npointsy,npointsx) )
        for i in range(npointsx):
            s[:,i]=zz

        # write Shadow file
        outFile = self.get_input_value("rootFile") + "Shadow.dat"
        tmp = write_shadowSurface(s,xx,yy,outFile=outFile)


    #
    # info
    #


    def info_profiles(self):
        # h = self._load_file_metadata()
        # sy,sz1,sz,zprof1,zprof = self._load_profiles()



        if int(self.get_input_value("setDetrending")) == -2: # this is the default
            if (self.h['SURFACE_SHAPE']).lower() == "elliptical":
                polDegree = -3     # elliptical detrending
            else:
                polDegree = 1      # linear detrending
        else:
            polDegree = int(self.get_input_value("setDetrending"))


        #;
        #; info
        #;
        #
        print ('\n---------- profile results --------------------')
        if (self.get_input_value("localFileRoot") == None):
            print ('Remote directory:\n   %s'%self.server)
        print ('Data File:     %s'%self.file_data())
        print ('Metadata File: %s'%self.file_metadata())
        print ('Surface shape: %s '%(self.h['SURFACE_SHAPE']))
        print ('Facility:      %s'%(self.h['FACILITY']))
        print ('Scan length: %.3f mm'%(1e3*(self.sy[-1]-self.sy[0])))
        print ('Number of points: %d'%(len(self.sy)))

        print ('   ')

        if polDegree >= 0:
            if polDegree == 1:
                print ("Linear detrending: z'=%g x%+g"%(self.coeffs[0],self.coeffs[1]))
                print ('Radius of curvature: %.3F m'%(1.0/self.coeffs[-2]))
            else:
                print ('Polynomial detrending coefficients: '+repr(self.coeffs))
        elif polDegree == -1:
           print ('No detrending applied.')
        elif polDegree == -3:
           print ('Ellipse detrending applied.')


        cdfHeightsRMS = (numpy.sqrt(cdf(self.f,self.psdHeights))).max()
        cdfSlopesRMS = (numpy.sqrt(cdf(self.f,self.psdSlopes))).max()

        # print ('   ')
        # print ('   ')
        # print ('Slope error s_RMS:             %.3f urad'%(1e6*self.sz.std()))
        # print ('         from PSD:             %.3f urad'%(1e6*cdfSlopesRMS))
        # print ('         from USER (metadata): %s urad'%(self.h['CALC_SLOPE_RMS']))
        # print ('Shape error h_RMS:              %.3f nm'%(1e9*self.zprof.std()))
        # print ('         from PSD:              %.3f nm'%(1e9*cdfHeightsRMS))
        # print ('         from USER (metadata):  %s nm'%(self.h['CALC_HEIGHT_RMS']))
        # print ('PV of height profile (before detrend): %.3f nm'%(  1e9*(self.zprof1.max() - self.zprof1.min() )))
        # print ('PV of height profile (after detrend):  %.3f nm'%(  1e9*(self.zprof.max() - self.zprof.min() )))
        # print ('   ')


        print (self.stdev_summary())


        # if (args.calcBoundaries == True):
        #     print ('Diffraction/reflection boundary zones: ')
        #     print ('   diffraction if h < h_cut, reflection if h > h_cut: ')
        #     print ('   at wavelength [A]: ')+repr(wavelength*1e10)
        #     print ('   at grazing angle [mrad]: ')+repr(theta*1e3)
        #     print ('   ')
        #     print ('   height cut h_cut [nm]: '+repr(  h_cut*1e9 ))
        #     print ('   h_cut/h_RMS: '+repr(  h_cut_normalized ))
        #     print ('   frequency cut f_cut [m^-1]: '+repr(  f_cut ))
        #     print ('   length 1/f_cut [mm]: '+repr(  1.0/f_cut ))

        #todo: remove
        print (self._latex_line())
        print ('-------------------------------------------------\n')
        print (' ')


    def plot(self):
        try:
            from matplotlib import pylab as plt
        except:
            raise ImportError

        what = self.get_input_value("plot")

        print("plotting: ",what)
        if (what == "heights" or what == "all"):
            f1 = plt.figure(1)
            plt.plot(1e3*self.sy,1e6*self.zprof)
            plt.title("heights profile")
            plt.xlabel("Y [mm]")
            plt.ylabel("Z [um]")
        if (what == "slopes" or what == "all"):
            f2 = plt.figure(2)
            plt.plot(1e3*self.sy,1e6*self.sz)
            plt.title("slopes profile" or what == "all")
            plt.xlabel("Y [mm]")
            plt.ylabel("Zp [urad]")
        if (what == "psd_h" or what == "all"):
            f3 = plt.figure(3)
            plt.loglog(self.f,self.psdHeights)
            plt.title("PSD of heights profile")
            plt.xlabel("f [m^-1]")
            plt.ylabel("PSD [m^3]")
        if (what == "psd_s" or what == "all"):
            f4 = plt.figure(4)
            plt.loglog(self.f,self.psdSlopes)
            plt.title("PSD of slopes profile")
            plt.xlabel("f [m^-1]")
            plt.ylabel("PSD [rad^3]")
        if (what == "cdf_h" or what == "all"):
            f5 = plt.figure(5)
            plt.semilogx(self.f,self.cdfHeights)
            plt.title("Lambda CDF(PDF) of heights profile")
            plt.xlabel("f [m^-1]")
            plt.ylabel("heights Lambda")
        if (what == "cdf_s" or what == "all"):
            f6 = plt.figure(6)
            plt.semilogx(self.f,self.cdfSlopes)
            plt.title("Lambda CDF(PDF) of slopes profile")
            plt.xlabel("f [m^-1]")
            plt.ylabel("slopes Lambda")
        plt.show()


    #
    # auxiliar methods for internal use
    #
    def _file_root(self):

        if self.remote_access():
            input_option = self.get_input_value("entryNumber")
            inFileRoot = 'dabam-'+str(input_option)
        else:
            inFileRoot = self.get_input_value("localFileRoot")

        return inFileRoot

    def _load_file_metadata(self):
        if self.remote_access():
            # metadata file
            myfileurl = self.server+self.file_metadata()
            u = urlopen(myfileurl)
            ur = u.read()
            ur1 = ur.decode(encoding='UTF-8')
            h = json.loads(ur1) # dictionnary with metadata
            self.h = h
        else:
            try:
                with open(self.file_metadata(), mode='r') as f1:
                    h = json.load(f1)
                self.h = h
            except:
                print ("_load_file_metadata: Error accessing local file: "+self.file_metadata())


    def _load_file_data(self):
        if self.remote_access():
            # data
            myfileurl = self.server+self.file_data()
            try:
                u = urlopen(myfileurl)
            except:
                print ("_load_file_data: Error accessing remote file: "+myfileurl+" does not exist.")
                return None

            ur = u.read()

            skipLines = self.h['FILE_HEADER_LINES']

            if sys.version_info[0] == 2:
                ur = StringIO(unicode(ur))
            else:
                ur = StringIO(ur.decode(encoding='ASCII'))

            a = numpy.loadtxt(ur, skiprows=skipLines )
            self.a = a
        else:
            try:
                with open(self.file_data(), mode='r') as f1:
                    h = json.load(f1)
                skipLines = h['FILE_HEADER_LINES']
                a = numpy.loadtxt(inFileDat, skiprows=skipLines) #, dtype="float64" )
                self.a = a
            except:
                print ("Error accessing local file: "+self.file_data())


    def _load_profiles(self):
        """
        Retrieve detrended profiles (slope and height): abscissa slope slope_detrended heights heights_detrended
        :return:
        """

        #;
        #; convert to SI units (m,rad)
        #;
        a = self.a.copy()

        a[:,0] = a[:,0]*self.h['X1_FACTOR']
        a[:,1] = a[:,1]*self.h['Y1_FACTOR']
        ncols = a.shape[1]

        if int(self.h["FILE_FORMAT"]) <= 2:
            for i in range(2,ncols):    # X1 Y1 Y2 Y3...
                a[:,i] = a[:,i]*h['Y%d_FACTOR'%i]
        elif int(self.h["FILE_FORMAT"]) == 3: #X1 Y1 X2 Y2 etc
            ngroups = int(ncols / 2)
            icol = 1
            for i in range(2,ngroups):    # X1 Y1 Y2 Y3...
                icol += 1
                a[:,icol] = a[:,icol]*self.h['X%d_FACTOR'%i]
                icol += 1
                a[:,icol] = a[:,icol]*self.h['Y%d_FACTOR'%i]

        #
        #; apply multiplicative factor
        #
        if (self.get_input_value("multiply") != 1.0):
            a[:,1] = a[:,1]  * float(get_value("multiply"))

        #
        # select columns with abscissas and ordinates
        #
        col_abscissas = int( self.get_input_value("useAbscissasColumn") )
        col_ordinates = int( self.get_input_value("useOrdinatesColumn") )


        col_ordinates_title = 'unknown'
        if int(self.get_input_value("useHeightsOrSlopes")) == -1:  #default, read from file
            if self.h['FILE_FORMAT'] == 1:  # slopes in Col2
                col_ordinates_title = 'slopes'
            if self.h['FILE_FORMAT'] == 2:  # heights in Col2
                col_ordinates_title = 'heights'
            if self.h['FILE_FORMAT'] == 3:  # slopes in Col2, file X1 Y1 X2 Y2
                col_ordinates_title = 'slopes'
        else:
            if int(self.get_input_value("useHeightsOrSlopes")) == 0:
                col_ordinates_title = 'heights'
            if int(self.get_input_value("useHeightsOrSlopes")) == 1:
                col_ordinates_title = 'slopes'

        if self.inputs["verbose"]:
            print("Using: abscissas column index %d (mirror coordinates)"%(col_abscissas))
            print("       ordinates column index %d (profile %s)"%(col_ordinates,col_ordinates_title))


        #;
        #; Detrending:
        #; substract linear fit to the slopes (remove best circle from profile)
        #;
        if col_ordinates_title == 'slopes':
            sy = a[:,col_abscissas]
            sz1 = a[:,col_ordinates]
        elif col_ordinates_title == 'heights':
            sy = a[:,col_abscissas]
            #TODO we suppose that data are equally spaced. Think how to generalise
            sz1 = numpy.gradient(a[:,col_ordinates],(sy[1]-sy[0]))
        else:
            raise NotImplementedError

        sz = numpy.copy(sz1)

        # define detrending to apply: >0 polynomial prder, -1=None, -2=Default, -3=elliptical
        if int(self.get_input_value("setDetrending")) == -2: # this is the default
            if (self.h['SURFACE_SHAPE']).lower() == "elliptical":
                polDegree = -3     # elliptical detrending
            else:
                polDegree = 1      # linear detrending
        else:
            polDegree = int(self.get_input_value("setDetrending"))


        if polDegree >= 0: # polinomial fit
            coeffs = numpy.polyfit(sy, sz1, polDegree)
            pol = numpy.poly1d(coeffs)
            zfit = pol(sy)
            sz = sz1 - zfit

        if polDegree == -3: # ellipse
            try:
                from scipy.optimize import curve_fit
            except:
                raise ImportError("Cannot perform ellipse detrending: please install scipy")

            popt, cov_x = curve_fit(func_ellipse_slopes, sy, sz1, maxfev=10000)
            zfit= func_ellipse_slopes(sy, popt[0], popt[1], popt[2], popt[3])
            sz = sz1 - zfit


        #;
        #; calculate heights by integrating the slope
        #;
        zprof = cdf(sy,sz)
        zprof1 = cdf(sy,sz1)

        self.sy = sy
        self.sz1 = sz1
        self.sz = sz
        self.zprof1 = zprof1
        self.zprof = zprof
        self.coeffs = coeffs


    def _calc_psd(self):
        sy    = self.sy
        sz1    = self.sz1
        sz    = self.sz
        zprof1    = self.zprof1
        zprof     = self.zprof

        #;
        #; calculate PSD on both profile and slope, and also then their cdf()
        #;
        psdHeights,f = psd(sy,zprof,onlyrange=None)
        psdSlopes,f = psd(sy,sz,onlyrange=None)
        cdfHeights = numpy.sqrt(cdf(f,psdHeights))
        cdfHeightsRMS = cdfHeights.max()
        cdfHeights = 1.0 - cdfHeights/cdfHeightsRMS
        cdfSlopes = numpy.sqrt(cdf(f,psdSlopes))
        cdfSlopesRMS = cdfSlopes.max()
        cdfSlopes = 1.0 - cdfSlopes/cdfSlopesRMS


        self.f = f
        self.psdHeights = psdHeights
        self.psdSlopes = psdSlopes
        self.cdfHeights = cdfHeights
        self.cdfSlopes = cdfSlopes
        return None


    def _latex_line(self):
        """
        to create a line with profile data latex-formatted for automatic compilation of tables in the paper
        :return:
        """
        return  ('%d & %d & %.2f (%.2f) & %.2f & %.2f (%.2f) & %.2f \\\\'%(self.get_input_value("entryNumber"),
        int(1e3*(self.sy[-1]-self.sy[0])),
        1e6*self.sz.std(),   float(-1 if self.h['CALC_SLOPE_RMS'] is None else self.h['CALC_SLOPE_RMS'])  , 1e6*self.stdev_psd_slopes(),
        1e9*self.zprof.std(),float(-1 if self.h['CALC_HEIGHT_RMS'] is None else self.h['CALC_HEIGHT_RMS']), 1e9*self.stdev_psd_heights(), ))


#
# main functions
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


#
# todo: reorganise this into class methods
#
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
        print ("File "+outFile+" written to disk (for SHADOW).")


#
# tests
#
def test_dabam_names():

    print("-------------------  test_dabam_names ------------------------------")
    dm = dabam()
    number_of_input_fields = len(dm.inputs)

    argsdict = dm.inputs
    names = []
    values = []
    for i,j in argsdict.items():
        names.append(i)
        values.append(j)

    #change values and reinsert in object
    values2 = copy.copy(values)
    for i in range(number_of_input_fields):
        if values[i] != None:
            values2[i] = 2*values[i]

    print ("-----------------------------------------------------")
    print ("--input_name value value2")
    for i in range(number_of_input_fields):
        print(i,names[i],values[i],values2[i])
        dm.inputs[names[i]] = values2[i]
    print ("-----------------------------------------------------")


    print ("-----------------------------------------------------")
    print ("--input_name input_name_short stored_value2, help")
    for i in range(number_of_input_fields):

        print(names[i],
            dm.get_input_value(names[i]),
            dm.get_input_value_short_name(names[i]),
            dm.inputs[names[i]],"\n",
            dm.get_input_value_help(names[i]),
              )
    print ("-----------------------------------------------------")


    #back to initial values
    dict2 = dm.get_inputs_as_dictionary()
    for i in range(number_of_input_fields):
        dict2[names[i]] = values[i]
        dm.inputs[names[i]] = values2[i]
    dm.set_inputs_from_dictionary(dict2)

    print ("--back to initial value")
    if (dm.inputs == dabam().inputs):
        print("Back to initial value: OK")
    else:
        raise Exception("Back to initial value: error returning to initial state")

def test_dabam_stdev_slopes():

    print("-------------------  test_dabam_slopes ------------------------------")
    stdev_profile_ok = [4.8651846142e-07,1.50962702525e-07]
    stdev_psd_ok = [2.47240356159e-08,2.44054273489e-09]

    for i in range(2):
        print(">> testing slopes stdev from profile number: ",i )
        dm = dabam()
        dm.set_input_verbose(False)
        dm.set_input_entryNumber(i+1)
        dm.load()
        stdev_profile = dm.stdev_profile_slopes()
        stdev_psd = dm.stdev_psd_slopes()
        print("stdev from profile ",stdev_profile,' OK: ',stdev_profile_ok[i])
        print("stdev from psd ",stdev_psd,' OK: ',stdev_psd_ok[i])
        assert abs(stdev_profile - stdev_profile_ok[i])<1e-10
        assert abs(stdev_psd - stdev_psd_ok[i])<1e-11

#
#
#
def main():

    # initialize
    dm = dabam()

    # get arguments of dabam command line
    dm.set_from_command_line()

    # access data
    dm.load()

    #dm.info_profiles()

    if dm.get_input_value("plot") != None:
        dm.plot()


    # test_dabam_names()
    #
    # test_dabam_stdev_slopes()
#
# main program
#
if __name__ == '__main__':
    main()