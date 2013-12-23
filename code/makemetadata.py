#
# simple script to create a template of DABAM metadata file
#
# 20131113  srio@esrf.eu 
#
import json
#import numpy
from collections import OrderedDict

# set the file index
for idx in range(4): 
    #a = {}
    a = OrderedDict()
    #a['FILE_ID'] = 2  # unique entry number in the database
    a['FILE_FORMAT'] = 1  #1: ESRF-slope: X1(mm),Y1(urad)
    a['FILE_HEADER_LINES'] = 4
    #a['FILE_HEADER_TEXT'] = "\n 1360mmx160mmx50mm \n Scan Length = 1200 mm \n measurement step = 1mm \n x(mm) Slope (urad)" # can be read automatically
    a['X1_FACTOR'] = 1e-3 #multiplicative factor from user units to SI length units [m]
    a['Y1_FACTOR'] = 1e-6 #multiplicative factor from user units to SI angle units [rad]
    a['YEAR_FABRICATION'] = None  # not available
    a['SURFACE_SHAPE'] = "plane"  
    a['FUNCTION'] = None
    a['LENGTH'] = 1360.0e-3
    a['WIDTH'] = 160.0e-3
    a['THICK'] = 50.0e-3
    a['LENGTH_OPTICAL'] = 1200.0e-3
    a['SUBSTRATE'] = "silicon" 
    a['COATING'] = None
    a['FACILITY'] = "ESRF"
    a['INSTRUMENT'] = None
    a['POLISHING'] = None
    a['ENVIRONMENT'] = None
    a['SCAN_DATE'] = None
    a['PLOT_TITLE_X1'] = "x(mm)"
    a['PLOT_TITLE_Y1'] = "Slope (urad)"
    a['CALC_HEIGHT_RMS'] = None
    a['CALC_HEIGHT_RMS_FACTOR'] = None  # for the detrended profile
    a['CALC_SLOPE_RMS'] = None
    a['CALC_SLOPE_RMS_FACTOR'] = None  # for the detrended profile
    a['USER_EXAMPLE'] = "This is an example of user keyword"
    #a['TEST_NUMPY'] = numpy.array([1,2,3]).tolist()
    
    if (idx == 1):
        a['LENGTH'] = 400.0e-3
        a['WIDTH'] = 70.0e-3
        a['THICK'] = 80.0e-3
        a['LENGTH_OPTICAL'] = 360.0e-3
        a['SURFACE_SHAPE'] = "plane"  
        
    if (idx == 2):
        a['LENGTH'] = 145.0e-3
        a['WIDTH'] = 45.0e-3
        a['THICK'] = 45.0e-3
        a['LENGTH_OPTICAL'] = 118.0e-3
        a['SURFACE_SHAPE'] = "spherical"  
        
    if (idx == 3):
        a['LENGTH'] = 40.0e-3
        a['WIDTH'] = 20.0e-3
        a['THICK'] = 20.0e-3
        a['LENGTH_OPTICAL'] = 32.0e-3
        a['SURFACE_SHAPE'] = "elliptical"  
        
    
    #
    # dump data into file
    #
    #filename = "dabam-"+str(a['FILE_ID'])+".txt"
    filename = "dabam-"+str(1+idx)+".txt"
    
    with open(filename, mode='w') as f: 
        json.dump(a, f, indent=2)
    
    print "File written to disk: "+filename
    
    #
    # for test purposes, read file into a new dictionnary
    #
    with open(filename, mode='r') as f: 
        b = json.load(f)
    
    #
    # list all non-empty keywords
    #
    print "-----------------------------------------------------"
    for i,j in a.items():
        if (j != None):
            print "%s = %s" % (i,j)
    print "-----------------------------------------------------"
    
    #keys = b.keys()
    #values = b.values()
    #indx = 0 
    #print "-----------------------------------------------------"
    #for i in keys:
    #    print keys[indx]," = ",values[indx]
    #    indx = indx + 1
    #print "-----------------------------------------------------"
    
    
