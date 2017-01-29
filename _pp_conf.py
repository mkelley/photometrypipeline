"""
configuration file for photometry pipeline

2016-01-27, michael.mommert@nau.edu
"""

import os
import sys
import logging
import warnings

# import pipeline-specific modules
from toolbox import *


def setup_diagnostics():
    """This function sets the current data and diagnostics root
    directories. The setup can and should be re-run in the current working
    directory"""

    # set up data and diagnostics directories
    dataroot = os.getcwd()+'/'
    diagroot = dataroot+'.diagnostics/'

    # create diagnostics directory
    if not os.path.exists(diagroot):
        os.mkdir(diagroot) 

    # define diagnostic website filenames
    index_filename = 'diagnostics.html'
    reg_filename   = '.diagnostics/' + 'registration.html'
    cal_filename   = '.diagnostics/' + 'calibration.html'
    res_filename   = '.diagnostics/' + 'results.html'


    return (dataroot, diagroot, index_filename,
            reg_filename, cal_filename, res_filename)



### suppress runtime and astropy warnings
warnings.simplefilter(action = "ignore", category = RuntimeWarning)
from astropy import wcs
from astropy.io import fits
warnings.filterwarnings('ignore', category=wcs.FITSFixedWarning)
warnings.filterwarnings('ignore', category=fits.column.VerifyWarning)
warnings.filterwarnings('ignore', category=fits.card.VerifyWarning)


### read photometry pipeline root path from environment variable
rootpath = os.environ.get('PHOTPIPEDIR')
if rootpath == None:
    print('ERROR: PHOTPIPEDIR variable has not been set')
    sys.exit(0)


# read telescope setup parameters (obsparam)
exec(compile(open(rootpath+'/setup/telescopes.py').read(), rootpath+'/setup/telescopes.py', 'exec'))


#### set diagnostics root directory and logging file

# master_diagroot is the directory where the pipeline has been started from
# i.e., potentially, the root directory of all underlying data
# this definition is outside setup_diagnostics(), as it is called only once
diagnostics_summary     = os.getcwd()+'/summary.html'
use_diagnostics_summary = False  # set this to True, if you want each 
                                 # pp_run process to report into a 
                                 # summary html catalog


# setting up directory paths and logging file
dataroot, diagroot, index_filename, \
    reg_filename, cal_filename, res_filename = setup_diagnostics()


### logging setup
log_formatline = '%(filename)s: %(message)s [%(asctime)s]'
log_level      = logging.DEBUG
log_datefmt    = '%m/%d/%Y %H:%M:%S'
log_filename = diagroot+'LOG'

# start pp_process_idx counter (if using 'pp_run all')
pp_process_idx = 0


### potential FITS header keywords for looking up the instrument
# any unique header keyword works as a potential identifier
instrument_keys = ['INSTRUME', 'LCAMMOD']


### available catalogs

# translate PP catalog identifier to Vizier identifier
# http://vizier.u-strasbg.fr/viz-bin/vizHelp?cats/U.htx
allcatalogs = {'URAT-1'  : 'urat1',
               '2MASS'   : '2mass-psc',
               'SDSS-R9' : 'sdss9',
               'APASS9'  : 'apass9',
               'CMC15'   : 'cmc15',
               'PPMXL'   : 'ppmxl',
               'USNO-B1' : 'usno-b1'}

# catalog magnitude systems
allcatalogs_magsys = {'URAT-1'  : 'Vega',
                      '2MASS'   : 'Vega',
                      'SDSS-R9' : 'AB',
                      'APASS9'  : 'Vega',
                      'CMC15'   : 'Vega',
                      'PPMXL'   : 'Vega',
                      'USNO-B1' : 'Vega'}

# catalog brightness fields for sorting
allcatalogs_mag = {'URAT-1'  : 'Vmag',
                   '2MASS'   : 'Jmag',
                   'SDSS-R9' : 'gmag',
                   'APASS9'  : 'Vmag',
                   'CMC15'   : 'r_mag',
                   'PPMXL'   : 'r1mag',
                   'USNO-B1' : 'r1mag'}




##### pipeline preferences 
# (if you don't know what you're doing, better don't mess around here)


# minimum number of reference sources in astrometric reference catalog
# in one of the frames (middle one of sequence); try other catalog
# if number of sources less than this number
min_sources_astrometric_catalog = 20

# how often to run SCAMP using one single catalog?
n_registration_repetitions = 2

# minimum number of reference sources in photometric reference catalog
min_sources_photometric_catalog = 3


# SCAMP contrast criteria for a good fit
scamp_as_contrast_limit = 2.5
scamp_xy_contrast_limit = 2.5

# positional uncertainty (arcsec) for target identification and
# cross-matching used in pp_photometry
pos_epsilon = 0.5

# flux threshold and margin for finding optimum aperture radius in
# pp_photometry
fluxlimit_aprad = 0.7
fluxmargin_aprad = 0.05

# minimum number of stars (integer number) or fraction (float) to use in
# photometric calibration
minstars = 0.5


#### support 

# path to local variable star database file
# if you 
vsx_database_file = os.environ.get('PPVARSTARSDB')
