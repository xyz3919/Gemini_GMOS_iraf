#! /usr/bin/env python

#import sys
import os
import numpy as np
from datetime import date
import copy
from pyraf import iraf
from pyraf.iraf import gemini, gemtools, gmos, onedspec, gnirs
import fileSelect as fs

# universal parameters

# VODKA
raw_path = '../raw'
obs = 'N' # N or S
observatory = 'Gemini-North'
extinction = 'gmos$calib/mkoextinct.dat'
ccdbin = '2 1'
roi = 'CentSp'
#apermask = '0.75arcsec'
apermask = '0.5arcsec'
sky_range = '350:450,575:625,750:850'
sky_range_no_offset = '350:450,575:625'
sky_range_offset = '575:625,750:850'
text_offset = '0 0\n0 -185.8\n0 -185.8\n0 0'

# DECam
"""
raw_path = '../raw'
obs = 'N' # N or S
ccdbin = '2 2'
#roi = 'Full'
roi = 'CentSp'
apermask = '0.75arcsec'
#sky_range = '980:1020,1070:1110,1160:1200'
sky_range = '190:230,300:340'
#text_offset = '0 0\n0 -92.9\n0 -92.9\n0 0'
"""

dbFile = '%s/obsLog.sqlite3' % raw_path

lambda1 = "4000" # lower bound
lambda2 = "10400" # upper bound
dx = "3.93" # A/pix

roi_std = 'CentSp'
#sky_range_std =  '190:230,300:340'
sky_range_std = '420:470,590:640'

centwave_std = 710
std_dir_name = 'onedstds$spec50cal/'
std_dir_name = 'gmos$calib/'
senfunc_std_name = 'eg131'


def get_month_range(central_date,date_range=(30,30)):
    from datetime import timedelta,date,datetime
    minus = timedelta(days=date_range[0])
    plus = timedelta(days=date_range[1])
    begin = date.isoformat(datetime.strptime(central_date,"%Y-%m-%d")-minus)
    end = date.isoformat(datetime.strptime(central_date,"%Y-%m-%d")+plus)
    return "%s:%s" %(begin,end)

def write_text_file(filename,text):

    f = open(filename, "w")
    f.write(text)
    f.close()


def move_fits(name,wavelengths,std=False):

    if not os.path.exists(name):
        os.makedirs(name)
    os.system("mv *gs*.fits %s/" % name)
    if not std: os.system("mv *%s*.fits %s/" % (name[0:5],name))

    if std: os.system("mv MCbiasstd.fits %s/" % name)
    else: os.system("mv MCbias.fits %s/" % name)

    for wave in wavelengths:
        os.system("mv MCflat%s.fits %s/" % (wave,name) )


def create_bias_files(exp_date,std=False):

    ###################
    # Creat Bias File #
    ###################
    if std: 
        suffix = "std"
        date_range = (7,7)
        roi_used = roi_std
    else: 
        suffix = ""
        date_range = (7,7)
        roi_used = roi


    # Select bias exposures within  months of the target observations:
    bias_time_range = get_month_range(exp_date,date_range=date_range)
    #for wave in wavelengths:
    #    qd[wave].update({'DateObs':bias_time_range})
    flag = {'use_me':1,'Instrument':'GMOS-%s' % obs,'RoI':roi_used,'CcdBin':ccdbin,\
            'DateObs':bias_time_range}

    print (" --Creating Bias MasterCal--")

    # Set the task parameters.
    gemtools.gemextn.unlearn()    # Disarm a bug in gbias
    gmos.gbias.unlearn()
    biasFlags = {'logfile':'biasLog.txt','rawpath':raw_path,\
                 'fl_vardq':'yes','verbose':'no'}

    # The following SQL generates the list of full-frame files to process.
    SQL = fs.createQuery('bias', flag)
    biasFiles = fs.fileListQuery(dbFile, SQL, flag)
    if len(biasFiles) > 1:
        print("Running gbias with %i bias files" % len(biasFiles))
        gmos.gbias(','.join(str(x) for x in biasFiles), 'MCbias'+suffix,**biasFlags)
    else: print("No Bias file found")

    # Clean up
    iraf.imdel('g%s20*.fits' % obs)

def create_flat_files(qd,wavelengths,std=False):

    ###################
    # Creat Flat file #
    ###################

    print (" -- Creating GCAL Spectral Flat-Field MasterCals --")

    if std: suffix = "std"
    else: suffix = ""

    gmos.gsflat.unlearn()
    # Normalize the spectral flats per CCD.
    # The response fitting should be done interactively.

    flatFlags = {
        'fl_over':'yes','fl_trim':'yes','fl_bias':'yes','fl_dark':'no',
        'fl_fixpix':'no','fl_oversize':'no','fl_vardq':'yes','fl_fulldq':'yes',
        'rawpath':raw_path,'fl_inter':'no','fl_detec':'yes',
        'function':'spline3','order':'13,11,28',
        'logfile':'gsflatLog.txt','verbose':'no'
        }
    for wave in wavelengths:
        flatFiles = fs.fileListQuery(dbFile, fs.createQuery('gcalFlat',\
                                     qd[wave]), qd[wave])
        if len(flatFiles) > 0:
            gmos.gsflat (','.join(str(x) for x in flatFiles), 'MCflat'+wave,
                    bias='MCbias'+suffix, **flatFlags)
        else: print("No Flat file found")

    # Clean up
    iraf.imdel('g%s20*.fits' % obs)


def basic_processing(qd,wavelengths,std=False):

    ####################
    # Basic Processing #
    ####################

    # apply bias and flat to arc and sci exposures

    print ("=== Processing Science Files ===")
    print (" -- Performing Basic Processing --")

    if std: suffix = "std"
    else: suffix = ""

    # Set task parameters.
    gmos.gsreduce.unlearn()
    sciFlags = {
        'fl_over':'yes','fl_trim':'yes','fl_bias':'yes','fl_gscrrej':'yes',
        'fl_dark':'no','fl_flat':'yes','fl_gmosaic':'yes','fl_fixpix':'yes',
        'fl_gsappwave':'yes','fl_oversize':'no','fl_crspec':'yes',
        'fl_vardq':'yes','fl_fulldq':'yes','rawpath':raw_path,
        'fl_inter':'no','logfile':'gsreduceLog.txt','verbose':'no'
    }
    arcFlags = copy.deepcopy(sciFlags)
    arcFlags.update({'fl_flat':'no','fl_vardq':'no','fl_fulldq':'no',
                     'fl_gscrrej':'yes','fl_crspec':'yes'})
    stdFlags = copy.deepcopy(sciFlags)

    # Perform basic reductions on all exposures for science targets.
    print ("  - Arc exposures -")
    for wave in wavelengths:
        arcFiles = fs.fileListQuery(dbFile, fs.createQuery('arc', qd[wave]),\
                                    qd[wave])
        if len(arcFiles) > 0:
            gmos.gsreduce (','.join(str(x) for x in arcFiles),
                           bias='MCbias'+suffix, **arcFlags)
    if std: img_type = 'std'
    else: img_type = 'sciSpec'

    print ("  - %s exposures -" % img_type)
    for wave in wavelengths:
        sciFiles = fs.fileListQuery(dbFile, fs.createQuery(img_type,
                                    qd[wave]), qd[wave])
        if len(sciFiles) > 0:
            gmos.gsreduce (','.join(str(x) for x in sciFiles),
                           bias='MCbias'+suffix, flatim='MCflat'+wave,**sciFlags)

    # Clean up
    iraf.imdel('g%s20*.fits' % obs)

def load_wavelength_cal(qd,wavelengths):

    #######################
    # Load Wavelength Cal #
    #######################

    # Set task parameters
    gmos.gswavelength.unlearn()
    waveFlags = {
        'coordlist':'gmos$data/CuAr_GMOS.dat','fwidth':3,'nsum':50,
        'function':'chebyshev','order':4,
        'fl_inter':'yes','logfile':'gswaveLog.txt','verbose':'no'
        }
    # Must select specific wavecals to match science exposures.

    prefix = 'gs'
    for wave in wavelengths:
        arcFiles = fs.fileListQuery(dbFile, fs.createQuery('arc', qd[wave]),\
                                    qd[wave])
        for arc in arcFiles:
            if not os.path.exists("database/id%s_001" % (prefix+arc)):
                gmos.gswavelength (prefix+arc, **waveFlags)
            else: print("arc id files exist for %s" % (prefix+arc))

def wavelength_calibration(qd,wavelengths,std=False):

    ##########################
    # wavelength calibration #
    ##########################

    prefix = 'gs'
    # Set task parameters.
    gmos.gstransform.unlearn()
    transFlags = {'fl_vardq':'yes','interptype':'spline3','fl_flux':'yes',\
            'logfile':'gstransLog.txt','lambda1':lambda1,'lambda2':lambda2,'dx':dx}

    if std: img_type = 'std'
    else: img_type = 'sciSpec'

    for wave in wavelengths:

        sciFiles = fs.fileListQuery(dbFile, fs.createQuery(img_type, qd[wave]),qd[wave])
        arcFiles = fs.fileListQuery(dbFile, fs.createQuery('arc', qd[wave]),qd[wave])
        for arc in arcFiles:
            gmos.gstransform (prefix+arc, wavtraname=prefix+arc, **transFlags)
        for sci in sciFiles:
            gmos.gstransform (prefix+sci, wavtraname=prefix+arcFiles[0], **transFlags)


    # Clean up
    iraf.imdel('%s%s20*.fits' % (prefix,obs))

def sky_subtraction(qd,wavelengths,std=False):

    ###################
    # sky subtraction #
    ###################

    prefix = 'tgs'
    gmos.gsskysub.unlearn()
    skyFlags = {'fl_oversize':'no','fl_vardq':'yes','fl_inter':'yes','verbose':'yes',\
            'logfile':'gsskysubLog.txt'
               }
    if std: 
        img_type = 'std'
        sky_range_used = sky_range_std
    else: 
        img_type = 'sciSpec'
        sky_range_used = sky_range

    for wave in wavelengths:

        # Fix up the target name for the output file
        sciFiles = fs.fileListQuery(dbFile, fs.createQuery(img_type, qd[wave]),qd[wave])
        arcFiles = fs.fileListQuery(dbFile, fs.createQuery('arc', qd[wave]),qd[wave])
        for sci in sciFiles:
            gmos.gsskysub (prefix+sci, long_sample=sky_range_used, **skyFlags)

def combine_img(qd,wavelengths,name,std=False):

    ##################
    # combine images #
    ##################

    prefix = 'stgs'
    gemtools.gemcombine.unlearn()
    offset_filename = 'offset.txt'
    combFlags = {'reject':'ccdclip','fl_vardq':'yes','verbose':'yes','offset':offset_filename}
    #combFlags = {'reject':'ccdclip','fl_vardq':'yes','verbose':'yes'}

    if std: img_type = 'std'
    else: img_type = 'sciSpec'

    write_text_file(offset_filename,text_offset)
    sciFiles_all = []
    for wave in wavelengths:
        sciFiles = fs.fileListQuery(dbFile, fs.createQuery(img_type, qd[wave]),qd[wave])
        for sciFile in sciFiles:
            sciFiles_all.append(prefix+sciFile)
    text_all_sci = ",".join(sciFiles_all)
    gemtools.gemcombine(text_all_sci, output=name[0:5], **combFlags)


def fit_sensitivity_function(qd):

    ##########################
    # sensitivity calibration #
    ###########################

    # extract 1d spectrum
    prefix = 'stgs'
    gmos.gsextract.unlearn()
    extrFlags = {
        'apwidth':3.,'fl_inter':'yes','find':'yes',
        'trace':'yes','tfunction':'chebyshev','torder':'6','tnsum':20,
        'background':'fit','bfunction':'chebyshev','border':2,
        'fl_vardq':'no','logfile':'gsextrLog.txt'
        }
    stdFiles = fs.fileListQuery(dbFile, fs.createQuery('std', qd['std']), qd['std'])
    for std in stdFiles:
        gmos.gsextract (prefix+std, **extrFlags)

    # fit sensitivity function
    prefix = 'estgs'
    gmos.gsstandard.unlearn()
    sensFlags = {
        'fl_inter':'yes','starname':senfunc_std_name,'caldir':std_dir_name,
        'observatory':observatory,'extinction':extinction,
        'function':'chebyshev','order':9,'verbose':'no','logfile':'gsstdLog.txt'
        }
    stdFiles = fs.fileListQuery(dbFile, fs.createQuery('std', qd['std']), qd['std'])
    for std in stdFiles:
        gmos.gsstandard (prefix+std, sfile='std.txt', sfunction='sens',**sensFlags)

    # avoid small mag fit at the edge -> causing 'Floating point overflow'
    #iraf.imreplace('sens',-30,lower='INDEF',upper=-30)


def flux_calibration(qd,wavelengths,name,std=False):

    ####################
    # Flux calibration #
    ####################

    prefix = 'stgs'
    gmos.gscalibrate.unlearn()
    calibFlags = {
        'extinction':extinction,'fl_ext':'yes','fl_scale':'no',
        'sfunction':'sens','fl_vardq':'yes','logfile':'gscalibrateLog.txt',
        'observatory':observatory}
    """
    # Science target
    for wave in wavelengths:
        sciFiles = fs.fileListQuery(dbFile, fs.createQuery('sciSpec', qd[wave]),qd[wave])
        for sci in sciFiles:
            gmos.gscalibrate(prefix+sci, **calibFlags)
    """
    if std :
        stdFiles = fs.fileListQuery(dbFile, fs.createQuery('std', qd['std']), qd['std'])
        gmos.gscalibrate(prefix+stdFiles[0], **calibFlags)
    else:
        gmos.gscalibrate(name[0:5], **calibFlags)
        extrFlags = {
        'apwidth':3.,'fl_inter':'yes','find':'yes',
        'trace':'yes','tfunction':'chebyshev','torder':'6','tnsum':20,
        'background':'fit','bfunction':'chebyshev','border':2,
        'fl_vardq':'no','logfile':'gsextrLog.txt'
        }
        gmos.gsextract ("c"+name[0:5], **extrFlags)




def reduce_std(name,exp_date):

    print ("### Begin Processing GMOS/Longslit Images ###")
    print ("=== Creating MasterCals ===")

    sci_time_range = get_month_range(exp_date,date_range=(1,1))

    qd = {"std":{'use_me':1,
           'Instrument':'GMOS-%s' % obs,'CcdBin':ccdbin,'RoI':roi_std,
           'Disperser':'R150+_%','CentWave':centwave_std,'AperMask':apermask,
           'Object':name+'%','DateObs':sci_time_range}
         }

    wavelengths = ["std"]

    # create bias files
    #create_bias_files(exp_date,std=True)

    # create flat files
    #create_flat_files(qd,wavelengths,std=True)

    # basic processing
    #basic_processing(qd,wavelengths,std=True)

    # load wavelength cal
    #load_wavelength_cal(qd,wavelengths)

    # wavelength calibration
    #wavelength_calibration(qd,wavelengths,std=True)

    # sky subtraction
    #sky_subtraction(qd,wavelengths,std=True)

    # fit sensitivity function
    fit_sensitivity_function(qd)

    # flux calibration
    flux_calibration(qd,wavelengths,name='std',std=True)

    # move fits files to sub dir
    move_fits(name,wavelengths,std=True)



def gmos_ls_proc(name,exp_date):

    print ("### Begin Processing GMOS/Longslit Images ###")
    print ("###")
    print ("=== Creating MasterCals ===")

    sci_time_range = get_month_range(exp_date,date_range=(1,1))

    qd = {'700':{'use_me':1,
           'Instrument':'GMOS-%s' % obs,'CcdBin':ccdbin,'RoI':roi,
           'Disperser':'R150+_%','CentWave':700.0,'AperMask':apermask,
           'Object':name+'%','DateObs':sci_time_range}
         }
    qd['720'] = copy.deepcopy(qd['700'])
    qd['720'].update({'CentWave':720.0})

    wavelengths = ['700','720']

    # create bias files
    #create_bias_files(exp_date)

    # create flat files
    #create_flat_files(qd,wavelengths)

    # basic processing
    #basic_processing(qd,wavelengths)

    # load wavelength cal
    #load_wavelength_cal(qd,wavelengths)

    # wavelength calibration
    #wavelength_calibration(qd,wavelengths)

    # sky subtraction
    #sky_subtraction(qd,wavelengths)

    # image combination and cosmic ray rejection
    #combine_img(qd,wavelengths,name)

    # flux calibration
    flux_calibration(qd,wavelengths,name)

    # move fits files to sub dir
    move_fits(name,wavelengths)

if __name__ == '__main__':

    #name = 'J1852+4833'
    #exp_date = '2020-07-02'
    #name = 'J1857+7048'
    #exp_date = '2020-07-29'
    #name = 'J1732-1335'
    #exp_date = '2020-08-13'
    #name = 'J1804+3230'
    #exp_date = '2020-08-14'
    #name = 'J214701.62+252315.3'
    #exp_date = '2020-07-02'
    #name = 'J220729.16+252203.3'
    #exp_date = '2020-07-30'
    #name = 'J225147.83+001640.6'
    #exp_date = '2020-07-31'
    #name = 'J093207.08+072251.9'
    #exp_date = '2020-02-04'
    #name = 'J081830.47+060138.1'
    #exp_date = '2020-02-03'
    name = 'J0749+2255'
    exp_date = '2022-04-20'
    std_name = 'EG131'
    std_exp_date = '2022-04-20'

    #name = 'J003715.89+205825.6'
    #exp_date = '2021-01-05'
    #name = 'J011812.03-010442.6'
    #exp_date = '2021-01-07'
    #name = 'J101142.54-030213.3'
    #exp_date = '2021-01-05'

    gmos_ls_proc(name,exp_date)
    #reduce_std(std_name,std_exp_date)
        

