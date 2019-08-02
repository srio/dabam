dabam
=====


** NOTE THAT THIS IS A WORKING REPOSITORY **
** data directory contains tha master repository of DABAM data **
** The official version of dabam.py is maintained at:
https://github.com/oasys-kit/SR-xraylib/blob/master/srxraylib/metrology/dabam.py **


DABAM: an open-source database of x-ray mirrors metrology


DABAM: an open-source database of x-ray mirrors metrology to be used with 
SHADOW or any other ray-tracing and wave propagation codes that can 
simulate the effect of the surface errors on the performance of a beamline.
        
This is the repository containing example files that can be used 
as templates for the database.
It also contains dabam.py, a python code that can be downloaded and run 
locally to retrieve the files and perform operations with them. 

We are seeking for collaborations with synchrotron facilities and metrology 
laboratories who can join the initiative by sending a few mirror metrology 
files in order to build a first useful database.

Please send datafiles, comments, ideas, contributions, etc. to 
Manuel Sanchez del Rio srio@esrf.eu

For discussion on this and other related topics, it is suggested to subscribe 
to the xrayoptics list serve mailing list:
https://lists.bnl.gov/mailman/listinfo/xrayoptics-l
                                                         

Important files in this directory: 
----------------------------------

doc/MetrologyDBProposal [.docx .pdf] Draft documentation of the project.

data/dabam-<N>.dat  datafile fir the Nth entry.

data/dabam-<N>.txt  metadata file for the Nth entry.

code/dabam.py  python program to access and manipulate the datafiles.  

code/makemetadata.py an auxiliary (not needed) python script to help creating 
    metadata files. 


