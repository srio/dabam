"""

dabam: (dataBase for metrology)
       python module for processing remote files containing the results
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
        self.description="dabam.py: python program to access and evaluate DAta BAse for Metrology (DABAM) files. See http://ftp.esrf.eu/pub/scisoft/dabam/README.md"

        self.server = "http://ftp.esrf.eu/pub/scisoft/dabam/data/"

        self.inputs = {
            'entryNumber':1,         # 'an integer indicating the DABAM entry number'
            'silent':False,          # 'Silent mode. Default is No'
            'localFileRoot':None,    # 'Define the name of local DABAM file root (<name>.dat for data, <name>.txt for metadata).'
            'outputFileRoot':"",     # 'Define the root for output files. Default is "", so no output files'
            'setDetrending':-2,      # 'Detrending: if >0 is the polynomial degree, -1=skip, -2=automatic, -3=ellipse. '
            'nbinS':101,             # 'number of bins of the slopes histogram in rads. '
            'nbinH':101,             # 'number of bins heights histogram in m. '
            'shadowCalc':False,      # 'Write file with mesh for SHADOW.'
            'shadowNy':-1,           # 'For SHADOW file, the number of points along Y (length). If negative, use the profile points. '
            'shadowNx':11,           # 'For SHADOW file, the number of points along X (width). '
            'shadowWidth':6.0,       # 'For SHADOW file, the surface dimension along X (width) in cm. '
            'multiply':1.0,          # 'Multiply input profile (slope or height) by this number (to play with StDev values). '
            'useHeightsOrSlopes':-1, # 'Force calculations using profile heights (0) or slopes (1). Overwrites FILE_FORMAT keyword. Default=-1 (like FILE_FORMAT)'
            'useAbscissasColumn':0,  # 'Use abscissas column index. '
            'useOrdinatesColumn':1,  # 'Use ordinates column index. '
            'plot':None,             # plot data
            'runTests':False,        # run tests cases
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
        self.histoSlopes = None # to store histogram

    #
    #setters
    #

    #variables
    def set_input_entryNumber(self,value):
        self.inputs["entryNumber"] = value
    def set_input_silent            (self,value):
        self.inputs["silent"] = value
    def set_input_localFileRoot      (self,value):
        self.inputs["localFileRoot"] = value
    def set_input_outputFileRoot     (self,value):
        self.inputs["outputFileRoot"] = value
    def set_input_setDetrending      (self,value):
        self.inputs["setDetrending"] = value
    def set_input_nbinS            (self,value):
        self.inputs["nbinS"] = value
    def set_input_nbinH            (self,value):
        self.inputs["nbinH"] = value
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
    def set_input_runTests           (self,value):
        self.inputs["runTests"] = value

    #others

    def set_from_command_line(self):
        #
        # define default aparameters taken from command arguments
        #
        parser = argparse.ArgumentParser(description=self.description)

        # main argument

        parser.add_argument('entryNumber', nargs='?', metavar='N', type=int, default=self.get_input_value('entryNumber'),
            help=self.get_input_value_help('entryNumber'))

        parser.add_argument('-'+self.get_input_value_short_name('runTests'), '--runTests', action='store_true',
            help=self.get_input_value_help('runTests'))

        # options (flags)

        parser.add_argument('-'+self.get_input_value_short_name('silent'),'--silent', action='store_true', help=self.get_input_value_help('silent'))

        #options (parameters)

        parser.add_argument('-'+self.get_input_value_short_name('localFileRoot'), '--localFileRoot', help=self.get_input_value_help('localFileRoot'))

        parser.add_argument('-'+self.get_input_value('outputFileRoot'), '--outputFileRoot', default=self.get_input_value('outputFileRoot'),
            help=self.get_input_value_help('outputFileRoot'))

        parser.add_argument('-'+self.get_input_value_short_name('setDetrending'), '--setDetrending', default=self.get_input_value('setDetrending'),
            help=self.get_input_value_help('setDetrending'))

        parser.add_argument('-'+self.get_input_value_short_name('nbinS'), '--nbinS', default=self.get_input_value('nbinS'),
            help=self.get_input_value_help('nbinS'))

        parser.add_argument('-'+self.get_input_value_short_name('nbinH'), '--nbinH', default=self.get_input_value('nbinH'),
            help=self.get_input_value_help('nbinH'))

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
        self.set_input_silent(args.silent)
        self.set_input_localFileRoot(args.localFileRoot)
        self.set_input_outputFileRoot(args.outputFileRoot)
        self.set_input_setDetrending(args.setDetrending)
        self.set_input_nbinS(args.nbinS)
        self.set_input_nbinH(args.nbinH)
        self.set_input_shadowCalc(args.shadowCalc)
        self.set_input_shadowNy(args.shadowNy)
        self.set_input_shadowNx(args.shadowNx)
        self.set_input_shadowWidth(args.shadowWidth)
        self.set_input_multiply(args.multiply)
        self.set_input_useHeightsOrSlopes(args.useHeightsOrSlopes)
        self.set_input_useAbscissasColumn(args.useAbscissasColumn)
        self.set_input_useOrdinatesColumn(args.useHeightsOrSlopes)
        self.set_input_plot(args.plot)
        self.set_input_runTests(args.runTests)


    def set_inputs_from_dictionary(self,dict):
        try:
            # self.inputs["entryNumber"]        =  dict["entryNumber"]
            # self.inputs["silent"]             =  dict["silent"]
            # self.inputs["localFileRoot"]      =  dict["localFileRoot"]
            # self.inputs["outputFileRoot"]     =  dict["outputFileRoot"]
            # self.inputs["setDetrending"]      =  dict["setDetrending"]
            # self.inputs["nbinS"]               =  dict["nbinS"]
            # self.inputs["nbinH"]               =  dict["nbinH"]
            # self.inputs["shadowCalc"]         =  dict["shadowCalc"]
            # self.inputs["shadowNy"]           =  dict["shadowNy"]
            # self.inputs["shadowNx"]           =  dict["shadowNx"]
            # self.inputs["shadowWidth"]        =  dict["shadowWidth"]
            # self.inputs["multiply"]           =  dict["multiply"]
            # self.inputs["useHeightsOrSlopes"] =  dict["useHeightsOrSlopes"]
            # self.inputs["useAbscissasColumn"] =  dict["useAbscissasColumn"]
            # self.inputs["useOrdinatesColumn"] =  dict["useOrdinatesColumn"]
            # self.inputs["plot"]               =  dict["plot"]
            # self.inputs["runTests"]           =  dict["runTests"]

            self.set_input_entryNumber        ( dict["entryNumber"]         )
            self.set_input_silent             ( dict["silent"]              )
            self.set_input_localFileRoot      ( dict["localFileRoot"]       )
            self.set_input_outputFileRoot     ( dict["outputFileRoot"]            )
            self.set_input_setDetrending      ( dict["setDetrending"]       )
            self.set_input_nbinS               ( dict["nbinS"]                )
            self.set_input_nbinH               ( dict["nbinH"]                )
            self.set_input_shadowCalc         ( dict["shadowCalc"]          )
            self.set_input_shadowNy           ( dict["shadowNy"]            )
            self.set_input_shadowNx           ( dict["shadowNx"]            )
            self.set_input_shadowWidth        ( dict["shadowWidth"]         )
            self.set_input_multiply           ( dict["multiply"]            )
            self.set_input_useHeightsOrSlopes ( dict["useHeightsOrSlopes"]  )
            self.set_input_useAbscissasColumn ( dict["useAbscissasColumn"]  )
            self.set_input_useOrdinatesColumn ( dict["useOrdinatesColumn"]  )
            self.set_input_plot               ( dict["plot"]                )
            self.set_input_runTests           ( dict["runTests"]            )
        except:
            raise Exception("Failed setting dabam input parameters from dictionary")

    #
    # tools
    #
    def is_remote_access(self):
        if (self.get_input_value("localFileRoot") == None):
            remoteAccess = 1  # 0=Local file, 1=Remote file
        else:
            remoteAccess = 0  # 0=Local file, 1=Remote file
        return remoteAccess

    def set_remote_access(self):
        self.set_input_localFileRoot(None)

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

        if key == 'entryNumber':        return 'An integer indicating the DABAM entry number or the remote profile files'
        if key == 'silent':             return 'Avoid printing information messages.'
        if key == 'localFileRoot':      return 'Define the name of local DABAM file root (<name>.dat for data, <name>.txt for metadata). If unset, use remore access'
        if key == 'outputFileRoot':     return 'Define the root for output files. Set to "" for no output.  Default is "'+self.get_input_value("outputFileRoot")+'"'
        if key == 'setDetrending':      return 'Detrending: if >0 is the polynomial degree, -1=skip, -2=automatic, -3=ellipse. Default=%d'%self.get_input_value("setDetrending")
        if key == 'nbinS':              return 'Number of bins for the slopes histogram in rads. Default is %d'%self.get_input_value("nbinS")
        if key == 'nbinH':              return 'Number of bins for the heights histogram in m. Default is %d'%self.get_input_value("nbinH")
        if key == 'shadowCalc':         return 'Write file with mesh for SHADOW. Default=No'
        if key == 'shadowNy':           return 'For SHADOW file, the number of points along Y (length). If negative, use the profile points. Default=%d'%self.get_input_value("shadowNy")
        if key == 'shadowNx':           return 'For SHADOW file, the number of points along X (width). Default=%d'%self.get_input_value("shadowNx")
        if key == 'shadowWidth':        return 'For SHADOW file, the surface dimension along X (width) in cm. Default=%4.2f'%self.get_input_value("shadowWidth")
        if key == 'multiply':           return 'Multiply input profile (slope or height) by this number (to play with StDev values). Default=%4.2f'%self.get_input_value("multiply")
        if key == 'useHeightsOrSlopes': return 'Force calculations using profile heights (0) or slopes (1). If -1, used metadata keyword FILE_FORMAT. Default=%d'%self.get_input_value("useHeightsOrSlopes")
        if key == 'useAbscissasColumn': return 'Use abscissas column index. Default=%d'%self.get_input_value("useAbscissasColumn")
        if key == 'useOrdinatesColumn': return 'Use ordinates column index. Default=%d'%self.get_input_value("useOrdinatesColumn")
        if key == 'plot':               return 'Plot: all heights slopes psd_h psd_s cdf_h cdf_s. histo_s histo_h. Default=%s'%repr(self.get_input_value("plot"))
        if key == 'runTests':           return 'Run test cases'
        return ''


    def get_input_value_short_name(self,key):

        if key == 'entryNumber':         return 'N'
        if key == 'silent':              return 's'
        if key == 'localFileRoot':       return 'l'
        if key == 'outputFileRoot':      return 'r'
        if key == 'setDetrending':       return 'D'
        if key == 'nbinS':               return 'b'
        if key == 'nbinH':               return 'e'
        if key == 'shadowCalc':          return 'S'
        if key == 'shadowNy':            return 'y'
        if key == 'shadowNx':            return 'x'
        if key == 'shadowWidth':         return 'w'
        if key == 'multiply':            return 'm'
        if key == 'useHeightsOrSlopes':  return 'Z'
        if key == 'useAbscissasColumn':  return 'A'
        if key == 'useOrdinatesColumn':  return 'O'
        if key == 'plot':                return 'P'
        if key == 'runTests':            return 'T'

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
            if self.get_input_value("entryNumber") <= 0:
                raise Exception("Error: entry number must be non-zero positive for remote access.")

        #calculate detrended profiles
        self._calc_detrended_profiles()

        #calculate psd
        self._calc_psd()

        #calculate histograms
        self._calc_histograms()

        #
        # write files
        #
        # write header file
        if self.get_input_value("outputFileRoot") != "":
            outFile = self.get_input_value("outputFileRoot") + "Header.txt"
            with open(outFile, mode='w') as f1:
                json.dump(self.h, f1, indent=2)
            if not(self.get_input_value("silent")):
                print ("File "+outFile+" containing metadata written to disk.")

        #
        # Dump heights and slopes profiles to files
        #
        if self.get_input_value("outputFileRoot") != "":
            outFile = self.get_input_value("outputFileRoot")+'Heights.dat'
            dd=numpy.concatenate( (self.sy.reshape(-1,1), self.zprof.reshape(-1,1)),axis=1)
            numpy.savetxt(outFile,dd,comments="#",header="F %s\nS 1  heights profile\nN 2\nL  coordinate[m]  height[m]"%(outFile))
            if not(self.get_input_value("silent")):
                print ("File "+outFile+" containing heights profile written to disk.")

            outFile = self.get_input_value("outputFileRoot")+'Slopes.dat'
            dd=numpy.concatenate( (self.sy.reshape(-1,1), self.sz.reshape(-1,1)),axis=1)
            numpy.savetxt(outFile,dd,comments="#",header="F %s\nS 1  slopes profile\nN 2\nL  coordinate[m]  slopes[rad]"%(outFile))
            if not(self.get_input_value("silent")):
                print ("File "+outFile+" written to disk.")


        #write psd file
        if self.get_input_value("outputFileRoot") != "":
            dd = numpy.concatenate( (self.f, self.psdHeights, self.psdSlopes, self.cdfHeights, self.cdfSlopes ) ,axis=0).reshape(5,-1).transpose()
            outFile = self.get_input_value("outputFileRoot")+'PSD.dat'
            header = "F %s\nS 1  power spectral density\nN 5\nL  freq[m^-1]  psd_heights[m^3]  psd_slopes[rad^3]  cdf(psd_h)  cdf(psd_s)"%(outFile)
            numpy.savetxt(outFile,dd,comments="#",header=header)
            if not(self.get_input_value("silent")):
                print ("File "+outFile+" written to disk.")


        # write slopes histogram
        if self.get_input_value("outputFileRoot") != "":
            dd=numpy.concatenate( (self.histoSlopes["x"],self.histoSlopes["y1"],self.histoSlopes["y2"] ) ,axis=0).reshape(3,-1).transpose()
            outFile = self.get_input_value("outputFileRoot")+'HistoSlopes.dat'
            numpy.savetxt(outFile,dd,header="F %s\nS  1  histograms of slopes\nN 3\nL  slope[rad] at bin center  counts  Gaussian with StDev = %g"%
                                            (outFile,self.stdev_profile_slopes()),comments='#')
            if not(self.get_input_value("silent")):
                print ("File "+outFile+" written to disk.")

        # heights histogrtam
        if self.get_input_value("outputFileRoot") != "":
            dd=numpy.concatenate( (self.histoHeights["x"],self.histoHeights["y1"],self.histoHeights["y2"] ) ,axis=0).reshape(3,-1).transpose()
            outFile = self.get_input_value("outputFileRoot")+'HistoHeights.dat'
            numpy.savetxt(outFile,dd,header="F %s\nS  1  histograms of heights\nN 3\nL  heights[m] at bin center  counts  Gaussian with StDev = %g"%
                                            (outFile,self.stdev_profile_heights()),comments='#')

        #shadow file
        if self.get_input_value("shadowCalc"):
            self.write_file_for_shadow()
            if not(self.get_input_value("silent")):
                outFile = self.get_input_value("outputFileRoot")+'Shadow.dat'
                print ("File "+outFile+" for SHADOW written to disk.")


        #info
        if not(self.get_input_value("silent")):
            print(self.info_profiles())

        if self.get_input_value("outputFileRoot") != "":
            outFile = self.get_input_value("outputFileRoot")+'Info.txt'
            f = open(outFile,'w')
            f.write(self.info_profiles())
            f.close()

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
        txt += 'Slope error s_StDev:           %.3f urad\n'       %( 1e6*self.stdev_profile_slopes() )
        txt += '         from PSD:             %.3f urad\n'       %( 1e6*self.stdev_psd_slopes() )
        txt += '         from USER (metadata): %s urad\n'         %(self.stdev_user_slopes())
        txt += 'Shape error h_StDev:           %.3f nm\n'        %(1e9*self.stdev_profile_heights() )
        txt += '         from PSD:             %.3f nm\n'        %(1e9*self.stdev_psd_heights() )
        txt += '         from USER (metadata): %s nm\n'          %(self.stdev_user_heights() )
        txt += 'PV of height profile (before detrend): %.3f nm\n' %(1e9*(self.zprof1.max() - self.zprof1.min() ))
        txt += 'PV of height profile (after detrend):  %.3f nm\n' %(1e9*(self.zprof.max() - self.zprof.min() ))

        return txt




    def _calc_histograms(self):
        """
        Calculates slopes and heights histograms and also the Gaussians with their StDev

        results are stored in:
        self.histoSlopes = {"x":hy_center, "y1":hz, "y2":g, "x_path":hy_path, "y1_path":hz_path, "y2_path":g_path}

        where:
          x is the abscissas (at bin center), y1 is the histogram, y2 is the Gaussian
          x_path is the abscissas with points at left and riggh edges of each bin, y1_path is the
        :return:
        """

        #
        # slopes histogram
        #

        # binsize = float(self.get_input_value("binS")) # default is 1e-7 rads
        # bins = numpy.ceil( (self.sz.max()-self.sz.min())/binsize )

        bins = int(self.get_input_value("nbinS"))
        hz,hy_left = numpy.histogram(self.sz, bins = bins)


        hy_center = hy_left[0:-1]+0.5*(hy_left[1]-hy_left[0]) #calculate positions of the center of the bins
        hy_right  = hy_left[0:-1]+1.0*(hy_left[1]-hy_left[0]) #calculate positions of the right edge of the bins

        hy_path = []
        hz_path = []
        for s,t,v in zip(hy_left,hy_right,hz):
            hy_path.append(s)
            hz_path.append(v)
            hy_path.append(t)
            hz_path.append(v)

        hy_path = numpy.array(hy_path)
        hz_path = numpy.array(hz_path)

        #Gaussian with StDev of data
        g = numpy.exp( -numpy.power(hy_center-self.sz.mean(),2)/2/numpy.power(self.stdev_profile_slopes(),2) )
        g = g/g.sum()*hz.sum()

        g_path = numpy.exp( -numpy.power(hy_path-self.sz.mean(),2)/2/numpy.power(self.stdev_profile_slopes(),2) )
        g_path = g_path/g_path.sum()*hz_path.sum()


        self.histoSlopes = {"x":hy_center, "y1":hz, "y2":g, "x_path":hy_path, "y1_path":hz_path, "y2_path":g_path}

        #
        # heights histogram
        #

        # binsize = float(self.get_input_value("binH"))
        # bins = numpy.ceil( (self.zprof.max()-self.zprof.min())/binsize )
        bins = int(self.get_input_value("nbinH"))
        hz,hy_left = numpy.histogram(self.zprof, bins = bins)

        hy_center = hy_left[0:-1]+0.5*(hy_left[1]-hy_left[0]) #calculate positions of the center of the bins
        hy_right  = hy_left[0:-1]+1.0*(hy_left[1]-hy_left[0]) #calculate positions of the right edge of the bins

        hy_path = []
        hz_path = []
        for s,t,v in zip(hy_left,hy_right,hz):
            hy_path.append(s)
            hz_path.append(v)
            hy_path.append(t)
            hz_path.append(v)

        hy_path = numpy.array(hy_path)
        hz_path = numpy.array(hz_path)

        #Gaussian with StDev of data
        g = numpy.exp( -numpy.power(hy_center-self.zprof.mean(),2)/2/numpy.power(self.stdev_profile_heights(),2) )
        g = g/g.sum()*hz.sum()

        g_path = numpy.exp( -numpy.power(hy_path-self.zprof.mean(),2)/2/numpy.power(self.stdev_profile_heights(),2) )
        g_path = g_path/g_path.sum()*hz_path.sum()

        self.histoHeights = {"x":hy_center, "y1":hz, "y2":g, "x_path":hy_path, "y1_path":hz_path, "y2_path":g_path}

    #
    # write things
    #
    def write_file_for_shadow(self):
        #
        #  write file for SHADOW (optional)
        #  replicate the (x,z) profile in a "s" mesh of npointsx * npointsy
        #

        #inputs
        npointsy = int(self.get_input_value("shadowNy"))
        npointsx = int(self.get_input_value("shadowNx"))
        mirror_width = float(self.get_input_value("shadowWidth")) # in cm

        # units to cm
        y = (self.sy).copy() * 100.0 # from m to cm
        z = (self.zprof).copy() * 100.0 # from m to cm

        # set origin at the center of the mirror. TODO: allow any point for origin
        z = z - z.min()
        y = y - y[int(y.size/2)]


        # interpolate the profile (y,z) to have npointsy points (new profile yy,zz)
        if npointsy > 0:
            mirror_length = y.max() - y.min()
            yy = numpy.linspace(-mirror_length/2.0,mirror_length/2.0,npointsy)
            zz = numpy.interp(yy,y,z)

            # dump to file interpolated profile (for fun)
            if self.get_input_value("outputFileRoot") != "":
                dd = numpy.concatenate( (yy.reshape(-1,1), zz.reshape(-1,1)),axis=1)
                outFile = self.get_input_value("outputFileRoot") + "ProfileInterpolatedForShadow.dat"
                numpy.savetxt(outFile,dd)
                if not(self.get_input_value("silent")):
                    print("File %s with interpolated heights profile for SHADOW written to disk."%outFile)
        else:
            yy = y
            zz = z
            npointsy = yy.size

        # fill the mesh arrays xx,yy,s with interpolated profile yy,zz
        xx=numpy.linspace(-mirror_width/2.0,mirror_width/2.0,npointsx)
        s = numpy.zeros( (npointsy,npointsx) )
        for i in range(npointsx):
            s[:,i]=zz

        # write Shadow file
        outFile = self.get_input_value("outputFileRoot") + "Shadow.dat"
        tmp = write_shadowSurface(s,xx,yy,outFile=outFile)


    #
    # info
    #
    def info_profiles(self):

        txt = ""

        if int(self.get_input_value("setDetrending")) == -2: # this is the default
            if (self.h['SURFACE_SHAPE']).lower() == "elliptical":
                polDegree = -3     # elliptical detrending
            else:
                polDegree = 1      # linear detrending
        else:
            polDegree = int(self.get_input_value("setDetrending\n"))



        #;
        #; info
        #;
        #
        txt += '\n---------- profile results -------------------------\n'
        if (self.get_input_value("localFileRoot") == None):
            txt += 'Remote directory:\n   %s\n'%self.server
        txt += 'Data File:     %s\n'%self.file_data()
        txt += 'Metadata File: %s\n'%self.file_metadata()
        txt += 'Surface shape: %s\n'%(self.h['SURFACE_SHAPE'])
        txt += 'Facility:      %s\n'%(self.h['FACILITY'])
        txt += 'Scan length: %.3f mm\n'%(1e3*(self.sy[-1]-self.sy[0]))
        txt += 'Number of points: %d\n'%(len(self.sy))

        txt += '   '

        if polDegree >= 0:
            if polDegree == 1:
                txt += "Linear detrending: z'=%g x%+g"%(self.coeffs[0],self.coeffs[1])+"\n"
                txt += 'Radius of curvature: %.3F m'%(1.0/self.coeffs[-2])+"\n"
            else:
                txt += 'Polynomial detrending coefficients: '+repr(self.coeffs)+"\n"
        elif polDegree == -1:
           txt += 'No detrending applied.\n'
        elif polDegree == -3:
           txt += 'Ellipse detrending applied.\n'

        txt += self.stdev_summary()

        txt += '----------------------------------------------------\n'
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


        return txt


    def plot(self):
        try:
            from matplotlib import pylab as plt
        except:
            print("Cannot make plots. Please install matplotlib.")
            return None

        what = self.get_input_value("plot")

        if what == "all":
            what = ["heights","slopes","psd_h","psd_s","cdf_h","cdf_s","histo_s","histo_h"]
        else:
            what = what.split(" ")

        for i,iwhat in enumerate(what):
            print("plotting: ",iwhat)
            if (iwhat == "heights" ):
                f1 = plt.figure(1)
                plt.plot(1e3*self.sy,1e6*self.zprof)
                plt.title("heights profile")
                plt.xlabel("Y [mm]")
                plt.ylabel("Z [um]")
            elif (iwhat == "slopes"):
                f2 = plt.figure(2)
                plt.plot(1e3*self.sy,1e6*self.sz)
                plt.title("slopes profile")
                plt.xlabel("Y [mm]")
                plt.ylabel("Zp [urad]")
            elif (iwhat == "psd_h"):
                f3 = plt.figure(3)
                plt.loglog(self.f,self.psdHeights)
                plt.title("PSD of heights profile")
                plt.xlabel("f [m^-1]")
                plt.ylabel("PSD [m^3]")
            elif (iwhat == "psd_s"):
                f4 = plt.figure(4)
                plt.loglog(self.f,self.psdSlopes)
                plt.title("PSD of slopes profile")
                plt.xlabel("f [m^-1]")
                plt.ylabel("PSD [rad^3]")
            elif (iwhat == "cdf_h"):
                f5 = plt.figure(5)
                plt.semilogx(self.f,self.cdfHeights)
                plt.title("Lambda CDF(PDF) of heights profile")
                plt.xlabel("f [m^-1]")
                plt.ylabel("heights Lambda")
            elif (iwhat == "cdf_s"):
                f6 = plt.figure(6)
                plt.semilogx(self.f,self.cdfSlopes)
                plt.title("Lambda CDF(PDF) of slopes profile")
                plt.xlabel("f [m^-1]")
                plt.ylabel("slopes Lambda")
            elif (iwhat == "histo_s" ):
                f7 = plt.figure(7)
                plt.plot(1e6*self.histoSlopes["x_path"],self.histoSlopes["y1_path"])
                plt.plot(1e6*self.histoSlopes["x_path"],self.histoSlopes["y2_path"])
                plt.title("slopes histogram and Gaussian with StDev: %10.3f urad"%(1e6*self.stdev_profile_slopes()))
                plt.xlabel("Z' [urad]")
                plt.ylabel("counts")
            elif (iwhat == "histo_h" ):
                f8 = plt.figure(8)
                plt.plot(1e9*self.histoHeights["x_path"],self.histoHeights["y1_path"])
                plt.plot(1e9*self.histoHeights["x_path"],self.histoHeights["y2_path"])
                plt.title("heights histogram and Gaussian with StDev: %10.3f nm"%(1e9*self.stdev_profile_heights()))
                plt.xlabel("Z [nm]")
                plt.ylabel("counts")
            else:
                print("Plotting options are: heights slopes psd_h psd_s cdf_h cdf_s")
                return None
        plt.show()


    #
    # auxiliar methods for internal use
    #
    def _file_root(self):

        if self.is_remote_access():
            input_option = self.get_input_value("entryNumber")
            inFileRoot = "dabam-"+"%03d"%(input_option)
        else:
            inFileRoot = self.get_input_value("localFileRoot")

        return inFileRoot

    def _load_file_metadata(self):
        if self.is_remote_access():
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
        if self.is_remote_access():
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


    def _calc_detrended_profiles(self):
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
                a[:,i] = a[:,i]*self.h['Y%d_FACTOR'%i]
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
                col_ordinates = 1
                col_ordinates_title = 'slopes'
            if self.h['FILE_FORMAT'] == 2:  # heights in Col2
                col_ordinates = 2
                col_ordinates_title = 'heights'
            if self.h['FILE_FORMAT'] == 3:  # slopes in Col2, file X1 Y1 X2 Y2
                col_ordinates = 1
                col_ordinates_title = 'slopes'
            if self.h['FILE_FORMAT'] == 4:  # heights in Col2, file X1 Y1 X2 Y2
                col_ordinates = 2
                col_ordinates_title = 'heights'
        else:
            if int(self.get_input_value("useHeightsOrSlopes")) == 0:
                col_ordinates_title = 'heights'
            if int(self.get_input_value("useHeightsOrSlopes")) == 1:
                col_ordinates_title = 'slopes'

        if not(self.get_input_value("silent")):
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
        else:
            coeffs = None

        if polDegree == -3: # ellipse
            coeffs = None
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
        cdfHeightsStDev = cdfHeights.max()
        cdfHeights = 1.0 - cdfHeights/cdfHeightsStDev
        cdfSlopes = numpy.sqrt(cdf(f,psdSlopes))
        cdfSlopesStDev = cdfSlopes.max()
        cdfSlopes = 1.0 - cdfSlopes/cdfSlopesStDev


        self.f = f
        self.psdHeights = psdHeights
        self.psdSlopes = psdSlopes
        self.cdfHeights = cdfHeights
        self.cdfSlopes = cdfSlopes


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

    b = numpy.sqrt( numpy.abs(p * q)) * numpy.sin(theta)
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
        #print ("File "+outFile+" written to disk (for SHADOW).")


#
# tests
#
def test_dabam_names():
    """
    Tests that the I/O methods work well for the list of input values
    :return:
    """

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
    """
    Tests the slope error value for the nmax first profiles (from remote server)
    :return:
    """

    print("-------------------  test_dabam_slopes ------------------------------")
    stdev_profile_ok = [4.8651846141972904e-07, 1.5096270252538352e-07, 1.7394444580303415e-07, 1.3428739185534941e-07, 8.4197811681221573e-07, 1.0097219914863226e-06, 5.74153915948042e-07, 5.7147678897188605e-07, 4.3527688789008779e-07, 2.3241765005153794e-07]
    stdev_psd_ok = [2.4724035615909884e-08, 2.4405427348911064e-09, 8.122795547837512e-10, 1.4019943864619925e-10, 1.6826933750035566e-08, 2.0711769898782215e-09, 3.2138903485739521e-10, 2.5457246098948428e-09, 1.9318084893422374e-09, 3.3646805118574403e-09]

    nmax = 10

    tmp_profile = []
    tmp_psd = []
    for i in range(nmax):
        print(">> testing slopes stdev from profile number: ",i )
        dm = dabam()
        dm.set_input_silent(True)
        dm.set_input_entryNumber(i+1)
        dm.load()
        stdev_profile = dm.stdev_profile_slopes()
        stdev_psd = dm.stdev_psd_slopes()
        tmp_profile.append(stdev_profile)
        tmp_psd.append(stdev_psd)

    print("stdev from profile:          ",repr(tmp_profile))
    print("stdev OK (stored) profile:   ",repr(stdev_profile_ok))
    print("stdev from psd:              ",repr(tmp_psd))
    print("stdev OK (stored)psd:        ",repr(stdev_psd_ok))

    for i in range(nmax):
        assert abs(tmp_profile[i] - stdev_profile_ok[i])<1e-10
        assert abs(tmp_psd[i] - stdev_psd_ok[i])<1e-11


#
#
#
def main():

    # initialize
    dm = dabam()

    dm.set_input_outputFileRoot("tmp") # write files by default
    dm.set_from_command_line()   # get arguments of dabam command line

    if dm.get_input_value("runTests"): # if runTests selected
        dm.set_input_outputFileRoot("")      # avoid output files
        test_dabam_names()
        test_dabam_stdev_slopes()
    else:

        dm.load()        # access data
        #todo: remove
        print(dm._latex_line())




        if dm.get_input_value("plot") != None:
            dm.plot()

#
# main program
#
if __name__ == '__main__':
    main()