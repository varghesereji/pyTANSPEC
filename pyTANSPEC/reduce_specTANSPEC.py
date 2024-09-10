#!/usr/bin/env python
#!/usr/anaconda3/envs/my_env_py3/bin python
#This script is to semi-automate basic data reduction of TANSPEC cross-dispersed spectroscopic data.
#Highlight ---> No 'IRAF' task has been used. Solely in Python.
#---------------------------------------------supriyo+indiajoe
import os
import os.path
import glob
import argparse
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import sys
import numpy as np
import warnings
import re
import shlex
import readline
import shutil
import subprocess
import time
import pkgutil
from astropy.stats import mad_std
from astropy.stats import biweight_location
from scipy import signal, ndimage

import fnmatch
import time
import datetime
from astropy.time import Time
from astropy.io import ascii
import astropy.table as table 
import astroscrappy
import configparser
from .libs import imarith
import WavelengthCalibrationTool.reidentify as reident
from WavelengthCalibrationTool import recalibrate
from WavelengthCalibrationTool.recalibrate import scale_interval_m1top1
import SpectrumExtractor.spectrum_extractor as specextractor

import warnings
warnings.filterwarnings("ignore")

def fxn():
    warnings.warn("deprecated", DeprecationWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()

try:
    from scipy.ndimage import filters
except ImportError:
    print('Scipy module missing.')
    print('You will need to install Scipy if and only if you need to do median filtered continuum division')
  
try:
    from functools32 import partial
    import ConfigParser
except ImportError:
    from functools import partial
    import configparser as  ConfigParser
    
       
def WriteSpecToFitsFile(SpecFlux, Wavel, fitsheader, VarianceF=None, Name1=None, Name2=None):
    # Creating a new file to store data (flux and wavelength)
    hdul1 = fits.PrimaryHDU(SpecFlux, header=fitsheader)
    if Name2 is None:
        hdul2 = fits.ImageHDU(Wavel)
    else:
        hdul2 = fits.ImageHDU(Wavel, name=Name2)
    if VarianceF is not None:
        hdul3 = fits.ImageHDU(VarianceF, name=Name1)
        hdulist = fits.HDUList([hdul1, hdul3, hdul2])
    else:
        hdulist = fits.HDUList([hdul1, hdul2])
    return hdulist  


def SpecMake(InputFiles, method = None, ScaleF = None):
    """ First scaled the spectra based on background windows and spectra extraction window
    then Combined if there are multiple wavelength calibrated spectra in the InputFiles.

    Parameters
    ----------
    InputFiles: list if number of file (N) is > 1
        List of filename of the input fits to combine
    method : {'sum','mean','median'}, optional 
        Combining method to use.
        sum: Sum of the spectra
        mean: Mean of the spectra
        median: Median of spectra

    ScaleF: scaling the spectra based on the background windows and spectra extraction window
    -------
    Outputhdulist: hdulist
        Returns output fits (if N > 1 then combined) file in flux and wavelength
        

    Examples
    --------
    >>> Outputhdulist = combine_frames(['InputFile1.fits','InputFile2.fits'],method='sum',
                                       ScaleF = 18/8)
    >>> Outputhdulist = combine_frames(['InputFile1.fits','InputFile2.fits'],method='mean',
                                       ScaleF = 18/8)

    """
    
    N = len(InputFiles)
    fitsheader = fits.getheader(InputFiles[0])
    SpecWavelF = fits.getdata(InputFiles[0], 6)
    
    SpecFluxF = []    
    if  N > 1:
        for i in range(N):
            Spechdul = fits.open(InputFiles[i]) #
            #Doing the background substraction and scaling
            SpecFlux = Spechdul[0].data - (np.mean([Spechdul[1].data, Spechdul[2].data], axis = 0)*(ScaleF))
            SpecFluxF.append(SpecFlux)
            
        if method == 'sum':
            SpecFluxF = np.sum(SpecFluxF, axis=0) #only fluxes are added
            fitsheader['SPECCOMBINE'] = (N, 'Number of frames combined')
            fitsheader['HISTORY'] = 'Sum of {0} frames'.format(N)
        
        elif method == 'mean':
            SpecFluxF = np.mean(SpecFluxF, axis=0) #mean fluxes 
            fitsheader['SPECCOMBINE'] = (N, 'Number of frames combined')
            fitsheader['HISTORY'] = 'Mean of {0} frames'.format(N)
            
        else:
            raise NotImplementedError('Unknown method {0}'.format(method))
    else:
        Spechdul = fits.open(InputFiles[0]) #
        #Doing the background substraction and scaling
        SpecFlux = Spechdul[0].data - (np.mean([Spechdul[1].data, Spechdul[2].data], axis = 0)*(ScaleF))
        if 'VARIANCE' in Spechdul:
            SpecVar = Spechdul['VARIANCE'].data + (Spechdul['BKG VARIANCE 0'].data + Spechdul['BKG VARIANCE 1'].data) * ScaleF**2/4
        SpecFluxF = SpecFlux
        VarianceF = SpecVar
        Spechdul.close()

    if len(VarianceF) == 0:
        VarianceF = None
    
    Outputhdulist = WriteSpecToFitsFile(SpecFluxF, Wavel=SpecWavelF, fitsheader=fitsheader, VarianceF=VarianceF, Name1 = 'VARIANCE', Name2 = 'WAVELENGTH')

    return Outputhdulist
        
#Spectral Extraction and Wavelneght Calibration
def SpectralExtraction_subrout(PC):
    if PC.TODO == 'SL':
        mode = 'LR'
    else:
        mode = 'XD'
        """ Extracts spectra from 2d image and also applies wavelength calibration """
    print('You are doing extraction of {} mode spectra'.format(mode))

    directories = LoadDirectories(PC,CONF=False)
    for night in directories:
        PC.currentnight = night # Updating the night directory for using PC.GetFullPath()
        print('Working on night: '+night)
        try:
            #First we load the Spectrastoextract_Lamp_BandFilter.txt
            SpecslistFILE=open(PC.GetFullPath('Spectrastoextract_Lamp_BandFilter.txt'),'r')
        except IOError:
            print('Could not find Spectrastoextract_Lamp_BandFilter.txt in this directory. Hence skipping %s directory'%(night))
            continue
        #Image to Lamp file dictionary
        Img2Lamp = []  #They are first born as lists...
        Img2NeLamp = []
        Img2Filt = []
        AllImglist = []
        for line in SpecslistFILE:
            line = shlex.split(line.rstrip())
            Img2Lamp.append((line[0],line[1]))
            Img2NeLamp.append((line[0],line[2]))
            Img2Filt.append((line[0],line[3]))
            AllImglist.append(line[0])
        SpecslistFILE.close()
        Img2Lamp = dict(Img2Lamp)
        Img2NeLamp = dict(Img2NeLamp)
        Img2Filt = dict(Img2Filt)
        #Initialising an empty dictionary for storing spectra of each filter
        Filt2finalspecs = dict()
        for filt in set(Img2Filt.values()) : Filt2finalspecs[filt] = []

        OutputObjSpecWlCaliList = [] # Initialising an empty array for storing wl calibrated spectra
        for img in AllImglist:
            print("Working on image "+night+" "+img)
            leftover = glob.glob(os.path.splitext(img)[0]+'_*.fits')  #Clean up of some previous attempts if any..
            leftover += glob.glob(os.path.splitext(img)[0]+'.ms.fits')
            if len(leftover) > 0 :
                if not os.path.isdir('Leftover'): os.makedirs('Leftover') 
                for lft in leftover :
                    shutil.move(lft,'Leftover/')

            # Spectrum extraction for object
            SpectrumFile = PC.GetFullPath(img)
            OutputObjSpec = SpectrumFile.rstrip('.fits')+'.ms.fits'
            #ReFitApertureInXD = 'p1'
            APERTUREWINDOW = PC.APERTUREWINDOW
            BKGWINDOWS = PC.BKGWINDOWS            
            ScaleFac = (APERTUREWINDOW[1] - APERTUREWINDOW[0]) / (BKGWINDOWS[1][1] - BKGWINDOWS[1][0])

            if PC.TODO == 'SX':
                OutputObjSpecWlCaliList = xdSpectralExtraction_subrout(PC,OutputObjSpecWlCaliList,SpectrumFile,OutputObjSpec,APERTUREWINDOW,BKGWINDOWS,night,
                                                                       Img2Lamp, Img2NeLamp, Img2Filt, Filt2finalspecs, img)
            elif PC.TODO == 'SL':
                OutputObjSpecWlCaliList = LrSpectralExtraction_subrout(PC,OutputObjSpecWlCaliList,SpectrumFile,OutputObjSpec,APERTUREWINDOW,BKGWINDOWS,night,
                                                                       Img2Lamp, Img2NeLamp, Img2Filt, Filt2finalspecs, img)
        N = len(OutputObjSpecWlCaliList)

        if PC.SCOMBINE == 'YES' and N> 1:
            OutputObjSpecWlCaliFinal = OutputObjSpecWlCaliList[0].rstrip('.fits') + '.final' + '.avg.fits'
            OutputObjSpecWlCaliFinalhdul = SpecMake(OutputObjSpecWlCaliList, method = 'mean',ScaleF=ScaleFac)
            OutputObjSpecWlCaliFinalhdul.writeto(OutputObjSpecWlCaliFinal, overwrite=True)
        elif PC.SCOMBINE == 'NO' and N> 1:
            for i in range(N):
                OutputObjSpecWlCaliFinal = OutputObjSpecWlCaliList[i].rstrip('.fits') + '.final.fits'
                OutputObjSpecWlCaliFinalhdul = SpecMake([OutputObjSpecWlCaliList[i]], method = None,ScaleF=ScaleFac)
                OutputObjSpecWlCaliFinalhdul.writeto(OutputObjSpecWlCaliFinal, overwrite=True)
        elif N == 1:
            OutputObjSpecWlCaliFinal = OutputObjSpecWlCaliList[0].rstrip('.fits') + '.final.fits'
            OutputObjSpecWlCaliFinalhdul = SpecMake(OutputObjSpecWlCaliList, method = None, ScaleF=ScaleFac)
            OutputObjSpecWlCaliFinalhdul.writeto(OutputObjSpecWlCaliFinal, overwrite=True)
        else:
            raise NotImplementedError('Unknown combine {0}'.format(PC.SCOMBINE))
        wlc_file = open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'Final_wlc_file.txt'),'w')
        wlc_file.write(OutputObjSpecWlCaliFinal)
        wlc_file.close()
    print('All nights over...')



def LrSpectralExtraction_subrout(PC,OutputObjSpecWlCaliList,SpectrumFile,OutputObjSpec,APERTUREWINDOW,BKGWINDOWS,night,Img2Lamp, Img2NeLamp, Img2Filt, Filt2finalspecs, img):
    # """ Extracts spectra from 2d image and also applies wavelength calibration """
                        
    #finding the configuration file for spectrum extraction as well as ContinuumFile and ApertureLabel as references
    pkgpath = os.path.split(pkgutil.get_loader('pyTANSPEC').get_filename())[0]
    ContinuumFile = os.path.join(pkgpath,'data/bg_removed_star.fits')
    ApertureLabel = os.path.join(pkgpath,'data/ApLabel_bg_removed_star.npy')
    ApertureTraceFilename = os.path.join(pkgpath,'data/TANSPEC_lr_trace.pkl')
    ConfigFileSpecExt = os.path.join(pkgpath,'config/spectrum_extractor_TANSPEC_lr.config') #Default config file

    #Finding if any 'spectrum_extractor_TANSPEC.config' file exist in working directory
    ConfigFileSpecExt_exists = os.path.exists(os.path.join(PC.RAWDATADIR,PC.SPECCONFIGFILE))

    if ConfigFileSpecExt_exists:
        ConfigFileSpecExt = PC.SPECCONFIGFILE #customized config file if the users want
        print('\n \033[1m' + 'Uses customized config file at {} users for spectrum extraction'.format(os.path.join(PC.RAWDATADIR,PC.SPECCONFIGFILE))  + '\033[0m')
    else:
        ConfigFileSpecExt = ConfigFileSpecExt
        print('\n \033[1m' + 'Uses default config file (https://github.com/astrosupriyo/pyTANSPEC/blob/main/pyTANSPEC/config/spectrum_extractor_TANSPEC_lr.config) for spectrum extraction' + '\033[0m')

    #creating a new configuration file 
    config = configparser.ConfigParser()
    config.optionxform = str  # preserve the Case sensitivity of keys
    config.read(ConfigFileSpecExt)

    #set these two variables, which are coming from the main TANSPEC confg file, into the spectrum_extractor_TANSPEC.config
    if len(config.get("tracing_settings","ContinuumFile").strip()) < 1:
        config.set("tracing_settings","ContinuumFile",str(ContinuumFile))
    if len(config.get("tracing_settings", "ApertureLabel").strip()) < 1:
        config.set("tracing_settings","ApertureLabel",str(ApertureLabel))
    if len(config.get("tracing_settings", "ApertureTraceFilename").strip()) < 1:
        config.set("tracing_settings","ApertureTraceFilename",str(ApertureTraceFilename))
    # if len(config.get("extraction_settings","ApertureWindow").strip()) < 1:
    config.set("extraction_settings","ApertureWindow",str(APERTUREWINDOW))
    # if len(config.get("extraction_settings","BkgWindows").strip()) < 1:
    config.set("extraction_settings","BkgWindows",str(BKGWINDOWS))

    #write new config file into the output directory.
    new_config_file_star = SpectrumFile.rstrip('.fits')+'.config'
    with open(new_config_file_star, 'w') as configfile:
        config.write(configfile)            
    # print(new_config_file_star) #Varghese added this
    OutputObjSpec, Avg_XD_shift, PixDomain = specextractor.main([SpectrumFile, new_config_file_star, OutputObjSpec])

    #plot the spectra for the first look
    hdulist = fits.open(OutputObjSpec)
    #Scale factor to scale the spectra
    ScaleFac = (APERTUREWINDOW[1] - APERTUREWINDOW[0]) / (BKGWINDOWS[1][1] - BKGWINDOWS[1][0])
    SigWithBkg = hdulist[0].data
    bkg = np.mean([hdulist[1].data, hdulist[2].data], axis=0)
    SigWithoutBkg = SigWithBkg - (bkg*ScaleFac)


    # Writing a configuration file for the Lamp
    ReFitApertureInXD = [tuple(Avg_XD_shift), tuple(PixDomain)]
    config.set("tracing_settings","ReFitApertureInXD",str(ReFitApertureInXD))

    #write new config file into the output directory.
    new_config_file_lamp = SpectrumFile.rstrip('.fits')+'.Lamp.config'
    with open(new_config_file_lamp, 'w') as configfile:
        config.write(configfile)

     # Spectrum extraction for Ar lamp             
    LampFile = os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,Img2Lamp[img])
    OutputArLampSpec = os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,img[:-5]+'_arc1.fits')            
    OutputArLampSpec, Avg_XD_shift, PixDomain = specextractor.main([LampFile, new_config_file_lamp,           
                                                                                                                        OutputArLampSpec])

    # Spectrum extraction for Ne lamp             
    NeLampFile = os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,Img2NeLamp[img])
    OutputNeLampSpec = os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,img[:-5]+'_arc2.fits')
    OutputNeLampSpec, Avg_XD_shift, PixDomain = specextractor.main([NeLampFile, new_config_file_lamp, 
                                                                                                                         OutputNeLampSpec])


    hdularc1 = fits.open(OutputArLampSpec)
    hdularc2 = fits.open(OutputNeLampSpec)

    #hdularc1data = hdularc1[0].data
    #hdularc2data = hdularc2[0].data
    hdularc1data = fits.getdata(OutputArLampSpec, ext=0)
    hdularc1header = fits.getheader(OutputArLampSpec, ext=0)
    hdularc2data = fits.getdata(OutputNeLampSpec, ext=0)
    hdularc2header = fits.getheader(OutputNeLampSpec, ext=0)
    #lamphdr = hdularc2[0].header
    hdularc1_slit = hdularc1header['Slit']
    hdularc2_slit = hdularc2header['Slit']

    if PC.INSTRUMENT == 'TANSPEC':
        if hdularc1_slit == hdularc2_slit: # matching slits for Ar and Ne lamps
            slit = hdularc1_slit 
        else:
            sys.exit("Slit width of lamps does not match: check header of both lamps. ")
        # for i in range(3):
        #     hdularc1data[i] = hdularc2data[i] # 1st 3-orders (from bluer end) of Ar replaced by Ne
    else:
           sys.exit("Script is for TANSPEC data reduction")     


    hdularc1[0].data = hdularc1data
    OutputCombLampSpec = OutputArLampSpec[:-9]+'combarc.fits'
    hdularc1.writeto(OutputCombLampSpec)

    #WL calibration by template matching
    pkgpath = os.path.split(pkgutil.get_loader('pyTANSPEC').get_filename())[0]
    # RefArFile = os.path.join(pkgpath,'data','LAMPIDENTDIR_LR',slit,'argon_lamp_template.npy' )
    # RefNeFile = os.path.join(pkgpath,'data','LAMPIDENTDIR_LR',slit,'neon_lamp_template.npy' )
    RefFile = os.path.join(pkgpath,'data','LAMPIDENTDIR_LR',slit,'tanspecArNeLR_template.npy' )
    # RefDispTableFile = os.path.join(pkgpath,'data','LAMPIDENTDIR',slit,'tanspecArNe' + '{}' + '.txt')
    # OutDispTableFile =  os.path.splitext(OutputCombLampSpec)[0] + '.OutDispTableFile' + '{}'
    OutputWavlFile =  os.path.splitext(OutputCombLampSpec)[0] + '.OutputWavlFile' + '{}' + '.npy'

    # ModelForDispersion =  PC.WLFITFUNC # Call recalibrate.main
    # _ = reident.main([OutputCombLampSpec, RefDispTableFile, OutDispTableFile, "--OutputWavlFile", OutputWavlFile,
    #                              "--ModelForDispersion", ModelForDispersion, "--SavePlots", "--StackOrders"])            
    Lamp_array = (fits.getdata(OutputArLampSpec)[0] + fits.getdata(OutputNeLampSpec)[0])
    # Template Matching
    # print('Template Matching')
    # print(fits.getdata(OutputArLampSpec))
    # print(RefFile)
    template = np.load(RefFile).T
    wl_soln, shift = recalibrate.ReCalibrateDispersionSolution(Lamp_array,
                                                               template.T,
                                                               method='p3')


    AllOutWlSolFile = OutputWavlFile.format('all')
    OutputObjSpechdul = fits.open(OutputObjSpec)
    AllOutWlSol = wl_soln
    #to make wavelength data as ImageHdu
    AllOutWlSolImgHDU = fits.ImageHDU(AllOutWlSol, name = 'Wavelength')
    OutputObjSpechdul.append(AllOutWlSolImgHDU)

    #wl calibrated spectra
    OutputObjSpecWlCali =  OutputObjSpec.rstrip('fits') + 'wlc.fits'
    print('wl calib file', OutputObjSpecWlCali)
    OutputObjSpechdul.writeto(OutputObjSpecWlCali)            
    OutputObjSpecWlCaliList.append(OutputObjSpecWlCali)                        
    return OutputObjSpecWlCaliList

def xdSpectralExtraction_subrout(PC,OutputObjSpecWlCaliList,SpectrumFile,OutputObjSpec,APERTUREWINDOW,BKGWINDOWS,night,Img2Lamp, Img2NeLamp, Img2Filt, Filt2finalspecs, img):
    """ Extracts spectra from 2d image and also applies wavelength calibration """

    #finding the configuration file for spectrum extraction as well as ContinuumFile and ApertureLabel as references
    pkgpath = os.path.split(pkgutil.get_loader('pyTANSPEC').get_filename())[0]
    ContinuumFile = os.path.join(pkgpath,'data/bd491296axd_s0.5Ahg-00031.Z.fits')
    ApertureLabel = os.path.join(pkgpath,'data/ApLabel_bd491296axd_s0.5Ahg-00031.Z.npy')
    ApertureTraceFilename = os.path.join(pkgpath,'data/bd491296axd_s0.5Ahg-00031.Z.fits_trace.pkl')
    ConfigFileSpecExt = os.path.join(pkgpath,'config/spectrum_extractor_TANSPEC.config') #Default config file

    #Finding if any 'spectrum_extractor_TANSPEC.config' file exist in working directory
    ConfigFileSpecExt_exists = os.path.exists(os.path.join(PC.RAWDATADIR,PC.SPECCONFIGFILE))

    if ConfigFileSpecExt_exists:
        ConfigFileSpecExt = PC.SPECCONFIGFILE #customized config file if the users want
        print('\n \033[1m' + 'Uses customized config file at {} users for spectrum extraction'.format(os.path.join(PC.RAWDATADIR,PC.SPECCONFIGFILE))  + '\033[0m')
    else:
        ConfigFileSpecExt = ConfigFileSpecExt
        print('\n \033[1m' + 'Uses default config file (https://github.com/astrosupriyo/pyTANSPEC/blob/main/pyTANSPEC/config/spectrum_extractor_TANSPEC.config) for spectrum extraction' + '\033[0m')

    #creating a new configuration file 
    config = configparser.ConfigParser()
    config.optionxform = str  # preserve the Case sensitivity of keys
    config.read(ConfigFileSpecExt)

    #set these two variables, which are coming from the main TANSPEC confg file, into the spectrum_extractor_TANSPEC.config
    if len(config.get("tracing_settings","ContinuumFile").strip()) < 1:
        config.set("tracing_settings","ContinuumFile",str(ContinuumFile))
    if len(config.get("tracing_settings", "ApertureLabel").strip()) < 1:
        config.set("tracing_settings","ApertureLabel",str(ApertureLabel))
    if len(config.get("tracing_settings", "ApertureTraceFilename").strip()) < 1:
        config.set("tracing_settings","ApertureTraceFilename",str(ApertureTraceFilename))
    # if len(config.get("extraction_settings","ApertureWindow").strip()) < 1:
    config.set("extraction_settings","ApertureWindow",str(APERTUREWINDOW))
    # if len(config.get("extraction_settings","BkgWindows").strip()) < 1:
    config.set("extraction_settings","BkgWindows",str(BKGWINDOWS))

    #write new config file into the output directory.
    new_config_file_star = SpectrumFile.rstrip('.fits')+'.config'
    with open(new_config_file_star, 'w') as configfile:
        config.write(configfile)            


    OutputObjSpec, Avg_XD_shift, PixDomain = specextractor.main([SpectrumFile, new_config_file_star, OutputObjSpec])
    
    #plot the spectra for the first look
    hdulist = fits.open(OutputObjSpec)
    #Scale factor to scale the spectra
    ScaleFac = (APERTUREWINDOW[1] - APERTUREWINDOW[0]) / (BKGWINDOWS[1][1] - BKGWINDOWS[1][0])
    SigWithBkg = hdulist[0].data
    bkg = np.mean([hdulist[1].data, hdulist[2].data], axis=0)
    SigWithoutBkg = SigWithBkg - (bkg*ScaleFac)
    plt.figure()
    plt.plot(SigWithoutBkg.flatten(), alpha=0.5)
    plt.title("Extracted Spectra from all {} orders".format((SigWithBkg.shape)[0]))
    plt.xlabel("Flattened pixels")
    plt.ylabel("Counts")
    plt.show()

    # Writing a configuration file for the Lamp
    ReFitApertureInXD = [tuple(Avg_XD_shift), tuple(PixDomain)]
    config.set("tracing_settings","ReFitApertureInXD",str(ReFitApertureInXD))
    config.set("tracing_settings","ShowPlot_Trace",str(False))

    #write new config file into the output directory.
    new_config_file_lamp = SpectrumFile.rstrip('.fits')+'.Lamp.config'
    with open(new_config_file_lamp, 'w') as configfile:
        config.write(configfile)

     # Spectrum extraction for Ar lamp             
    LampFile = os.path.join(PC.OUTDIR,night,Img2Lamp[img])
    OutputArLampSpec = os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,img[:-5]+'_arc1.fits')            
    OutputArLampSpec, Avg_XD_shift, PixDomain = specextractor.main([LampFile, new_config_file_lamp,           
                                                                                                                        OutputArLampSpec])

    # Spectrum extraction for Ne lamp             
    NeLampFile = os.path.join(PC.OUTDIR,night,Img2NeLamp[img])
    OutputNeLampSpec = os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,img[:-5]+'_arc2.fits')
    OutputNeLampSpec, Avg_XD_shift, PixDomain = specextractor.main([NeLampFile, new_config_file_lamp, 
                                                                                                                         OutputNeLampSpec])

    # Finding pixel offset
    template_filename = os.path.join(pkgpath, 'data/PIXELOFFSETTEMPLATES', 'pixeloffsettemplate_s0.5.npy')
    template_file = np.load(template_filename)
    Ar_spec = fits.getdata(OutputArLampSpec, ext=0)
    arc_lamp = Ar_spec[7:9].flatten()
    arc_filtered = ndimage.gaussian_filter(signal.medfilt(arc_lamp,3), sigma=10, radius=20)
    PixShiftGuess = recalibrate.calculate_pixshift_with_phase_cross_correlation(template_file, arc_filtered)
    print('This night have a pixel offset of', PixShiftGuess)

    #Combine Ar (7 orders from the redder end) and Ne (3 orders from the bluer end) lamps
    hdularc1 = fits.open(OutputArLampSpec)
    hdularc2 = fits.open(OutputNeLampSpec)

    #hdularc1data = hdularc1[0].data
    #hdularc2data = hdularc2[0].data
    hdularc1data = fits.getdata(OutputArLampSpec, ext=0)
    hdularc1header = fits.getheader(OutputArLampSpec, ext=0)
    hdularc2data = fits.getdata(OutputNeLampSpec, ext=0)
    hdularc2header = fits.getheader(OutputNeLampSpec, ext=0)
    #lamphdr = hdularc2[0].header
    hdularc1_slit = hdularc1header['Slit']
    hdularc2_slit = hdularc2header['Slit']

    if PC.INSTRUMENT == 'TANSPEC':
        if hdularc1_slit == hdularc2_slit: # matching slits for Ar and Ne lamps
            slit = hdularc1_slit 
        else:
            sys.exit("Slit width of lamps does not match: check header of both lamps. ")
        hdularc2data = (hdularc1data + hdularc2data)/2 #Simply adding up both the spectra. This will be used for template matching.    
        # for i in range(3):
        #     hdularc1data[i] = hdularc2data[i] # 1st 3-orders (from bluer end) of Ar 
        #   replaced by Ne
        # hdularc2data[1] = hdularc1data[1]
        # hdularc2data[3] = hdularc1data[3] # Order 1 and 3 (from bluer end) of Ne will replaced by Ar. From here, this will be hdularc1.
        # if slit == 'S-4.0':
        #     hdularc2data[5] = hdularc1data[5] # Replacing order 5 of neon with argon. Neon does not give a good fit in this case.
    else:
           sys.exit("Script is for TANSPEC data reduction")     

    #            if PC.INSTRUMENT == 'TANSPEC':
    #                for i in range(3):
    #                    hdularc1data[i] = hdularc2data[i] # 1st 3-orders (from bluer end) of Ar replaced by Ne 

    hdularc1[0].data = hdularc2data
    OutputCombLampSpec = OutputArLampSpec[:-9]+'combarc.fits'
    hdularc1.writeto(OutputCombLampSpec)

    # WL calibration by template matching
    pkgpath = os.path.split(pkgutil.get_loader('pyTANSPEC').get_filename())[0]
    RefDispTableFile = os.path.join(pkgpath,'data','LAMPIDENTDIR',slit,'tanspecArNe' + '{}' + '_template.npy')
    OutputWavlFile =  os.path.splitext(OutputCombLampSpec)[0] + '.OutputWavlFile' + '{}' + '.npy'
    ModelForDispersion = 'p3' # PC.WLFITFUNC
    OutputWavlFile = None
    for order in range(0, 10):
        print('\033[31m Doing wavelength calibration for aperture {} \033[0m'.format(order))
        lamp = hdularc2data[order]
        template = np.load(RefDispTableFile.format(order))
        soln, shift = recalibrate.ReCalibrateDispersionSolution(lamp,
                                                                template.T,
                                                                method='p3',
                                                                initial_guess=[1,
                                                                -PixShiftGuess*np.median(np.gradient(template.T[:,0]))*2/(max(template.T[:,0])-min(template.T[:,0])),
                                                                1,0,0]
                                                                )
        
        plt.figure()
        plt.plot(soln, lamp/np.max(lamp), label='Observed lamp (arg+neo)')
        plt.plot(template[0], template[1]/np.max(template[1]), label='Template')
        plt.legend()
        template_match_filename = os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'template_match_aperture{}.pdf'.format(order))
        plt.title('Aperture {}'.format(order))
        print('Sucessfully fitted, saving', template_match_filename)
        plt.savefig(template_match_filename)
        if OutputWavlFile is None:
            OutputWavlFile = soln
        else:
            OutputWavlFile = np.vstack((OutputWavlFile, soln))
    #to make wavelength data as ImageHdu
    AllOutWlSol = OutputWavlFile
    OutputObjSpechdul = fits.open(OutputObjSpec)
    AllOutWlSolImgHDU = fits.ImageHDU(AllOutWlSol, name = 'Wavelength')
    OutputObjSpechdul.append(AllOutWlSolImgHDU)

    # Plotting sky
    sky_line = os.path.join(pkgpath, 'data/Sky_spectra/Sky_spectra_xd.fits')
    sky_spec = fits.getdata(sky_line, ext=0)
    sky_wl = fits.getdata(sky_line, ext=1)
    
    sky1 = OutputObjSpechdul[1].data
    sky2 = OutputObjSpechdul[2].data

    avg_sky = (sky1 + sky2) / 2

    l = np.shape(avg_sky)[0]
    fig, axs = plt.subplots(l//2, l//(l//2))
    for i in range(l):
        sky_o = avg_sky[i] / np.median(avg_sky[i])
        sky_spec_o = sky_spec[i] / np.median(sky_spec[i])
        axs[i//2, i%2].plot(sky_wl[i], sky_spec_o, label='Std sky')
        axs[i//2, i%2].plot(OutputWavlFile[i], sky_o, label='observed sky')
        axs[i//2, i%2].set(xlabel='wavelength', ylabel='flux', ylim=(0.2, np.max(sky_spec_o)))
        
    plt.legend()
    plt.savefig(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'sky_spectra.pdf'))
                           
    #wl calibrated spectra
    OutputObjSpecWlCali =  OutputObjSpec.rstrip('fits') + 'wlc.fits'
    print(OutputObjSpecWlCali)
    OutputObjSpechdul.writeto(OutputObjSpecWlCali)            
    OutputObjSpecWlCaliList.append(OutputObjSpecWlCali)                        
    return OutputObjSpecWlCaliList   

    
def FluxCalibration_subrout(PC):
    """This will do the flux calibration using instrument response function."""
    directories = LoadDirectories(PC,CONF=False)
    pkgpath = os.path.split(pkgutil.get_loader('pyTANSPEC').get_filename())[0]
    
    # Importing the response function
    if PC.TODO == 'P':
        print('Flux calibration of photometry is yet to define in this pipeline')
    elif PC.TODO == 'SX':
        if PC.RESPFN == 'D':
            response_function_name = os.path.join(pkgpath, 'data', 'INSTRUMENT_RESPONSE','XD_Response_curve.npy')
        else:
            response_function_name = PC.RESPFN
        response_function = np.load(response_function_name)
        # if response_function.shape()[0] == :
        #     print('This is the instrument response for LR mode. You are working of XD mode spectra')
    elif PC.TODO == 'SL':
        if PC.RESPFN == 'D':
            response_function_name = os.path.join(pkgpath, 'data', 'INSTRUMENT_RESPONSE','LR_Response_curve.npy')
        else:
            response_function_name = PC.RESPFN
        response_function = np.load(response_function_name)
    for night in directories:
        with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'Final_wlc_file.txt'),'r') as wlc_data_files:
            wlc_data = list([data for data in wlc_data_files])
        for data_file in wlc_data:
#            print(data_file)
            star_data = fits.getdata(data_file, ext=0)
            star_header = fits.getheader(data_file)
            star_header['HISTORY'] = 'Flux calibration was done'

            Spechdul = fits.open(data_file)
            wl_data = Spechdul['WAVELENGTH'].data
            FLC_variance = None
            varname = None

            FLC_data = star_data / response_function
            if 'VARIANCE' in Spechdul:
                variance = Spechdul['VARIANCE'].data
                FLC_variance = variance / response_function ** 2
                varname = 'VARIANCE'

            OutputObjSpecFlCalFinal = data_file.rstrip('.fits') + '.flc.fits'

            hdulist = WriteSpecToFitsFile(FLC_data, wl_data, star_header, FLC_variance, Name1=varname, Name2='WAVELENGTH')
            hdulist.writeto(OutputObjSpecFlCalFinal)
            print('The flux calibrated spectra is saved as:',OutputObjSpecFlCalFinal)

def SpectralPairSubtraction_subrout(PC):
    """ This will display all the spectra to be extracted to choose the mutual subtraction pairs """
    directories = LoadDirectories(PC,CONF=False)
    for night in directories:
        PC.currentnight = night # Upgating the night directory for using PC.GetFullPath()
        print('Working on night: '+night)
        try:
            #First we load a dictionary of raw images to their filters
            with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects.List'),'r') as FiltrFILE :
                Filtrfiledic = dict([(filtset.split()[0],shlex.split(filtset.rstrip())[1]) for filtset in FiltrFILE])  #Dictionary of filterset for each image.
                               
            #Secondly we load a dictionary of raw images to their Calibration Lamp  file
            #with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects-BSFinalLamp.List'),'r') as LampFILE :
            with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects-FinalLamp.List'),'r') as ArLampFILE:
                ArLampfiledic = dict([(Lampset.split()[0],shlex.split(Lampset.rstrip())[-2]) for Lampset in ArLampFILE])  #Dictionary of Ar Lamp file for each image.
            with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects-FinalLamp.List'),'r') as NeLampFILE :                
                NeLampfiledic = dict([(Lampset.split()[0],shlex.split(Lampset.rstrip())[-1]) for Lampset in NeLampFILE]) #Dictionary of Ne Lamp file for each image.
            #Secondly, we load a dictionary of Dither Frame to first Raw images
            ProcessedImagesfilename = 'AllObjects-ProcessedCombinedImg.List' if PC.IMGCOMBINE == 'Y' else 'AllObjects-ProcessedImg.List'
            with open(PC.GetFullPath(ProcessedImagesfilename),'r') as DitherFILE :
                DitherFILElines = list(DitherFILE)

            Ditherfiledic = dict([(ditherset.rstrip().split()[1],ditherset.split()[0]) for ditherset in DitherFILElines if len(ditherset.split()) == 2])  #Dictionary of First image of each Dither set.
            #We shall also keep a list of images in proper order.
            Allimglist = [ditherset.rstrip().split()[1]  for ditherset in DitherFILElines if len(ditherset.split()) == 2]

        except IOError as e:
            print('Cannot open the image file list.')
            print(e)
            print('So skipping this directory.')
            print('-'*60)
            continue
        #Output file to write the table of final image, lamp and filter band
        outlog = open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'Spectrastoextract_Lamp_BandFilter.txt'),'w')
        #First Lets create the list of filters to iterate through
        FilterList = list(set(Filtrfiledic.values()))
        print("You have %d orders of spectra for this object on this night %s"%(len(FilterList),night))
        OutFilePrefix = input('Enter the prefix of you want for reduced 1d spectra: ')
        for filt in FilterList:
            print("Working on filter : "+filt)
            filtstr = filt.replace(" ","_") # Replace spaces in filter name for easier filenames
            #List of images with this filter.
            Imglist = [img for img in Allimglist if Filtrfiledic[Ditherfiledic[img]] == filt ]            
       
            for i,img in enumerate(Imglist): 
                ImagePlot(fits.getdata(PC.GetFullPath(img)), title=img) #to display the image
#                pass

            if len(Imglist) != 0:
                ABCDtoimg = dict()
                Anumb = ord('A')
                with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'ABCDtoImageTable_'+filtstr+'.txt'),'w') as AlphatoFILE :
                    for i,img in enumerate(Imglist):
                        alpha = chr(Anumb+i)
                        ABCDtoimg[alpha] = img
                        print("%s : %s"%(alpha,img))
                        AlphatoFILE.write("%s  %s"%(alpha,img)+' \n')

                print('-'*30)
                print("Enter the pairs to subtract in space separated form ")
                print("For Example an input: AB BA A")
                print("Corresponding images produced by subtraction or not are : A-B, B-A, and A")
                print("Note: the final A is not a subtracted image ")
                subpairs = input('Pairs to process: ')
                subpairs = subpairs.split()
                for instr in subpairs:
                    instr = instr.upper()
                    if len(instr) == 2 :
                        Outimg = OutFilePrefix+'_'+filtstr+'_'+instr[0]+'-'+instr[1]+'.fits'
    #                    try:                            #iraf.imarith(operand1=PC.GetFullPath(ABCDtoimg[instr[0]]),op="-",operand2=PC.GetFullPath(ABCDtoimg[instr[1]]),result=PC.GetFullPath(Outimg))
                        Outputhdulist = imarith.subtract_frames(PC.GetFullPath(ABCDtoimg[instr[0]]), PC.GetFullPath(ABCDtoimg[instr[1]]), FluxExt=[0])
                        if PC.INSTRUMENT == 'TANSPEC':
                            Outputhdulist[0].data = np.clip(Outputhdulist[0].data, -1, None) 
                        imarith.WriteFitsOutput(Outputhdulist,PC.GetFullPath(Outimg),overwrite=True)
     
#                        except iraf.IrafError as e :
#                            print(e)
#                            print('Skipping to next instruction')
#                            continue
                    elif len(instr) == 1 : 
                        Outimg = OutFilePrefix+'_'+filtstr+'_'+instr[0]+'.fits'
                        shutil.copy(PC.GetFullPath(ABCDtoimg[instr[0]]),PC.GetFullPath(Outimg))
                    else : 
                        print("Could not understand "+instr)
                        continue
                    # print(Outimg)
                    outlog.write('{0} {1} {2} "{3}" \n'.format(Outimg,ArLampfiledic[Ditherfiledic[ABCDtoimg[instr[0]]]],NeLampfiledic[Ditherfiledic[ABCDtoimg[instr[0]]]],filt))
                    
                #Now Copy the already identified Repository Lamps to this directory.
#                print("Copying already identified lines of this filter %s from Repository.."%(filt))
#                if not os.path.isdir(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'database')): os.makedirs(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'database'))
#                try:
#                    shutil.copy(PC.LAMPREPODIR+'/RepoLamp_'+filtstr+'.fits',os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'RepoLamp_'+filtstr+'.fits'))
#                    shutil.copy(PC.LAMPREPODIR+'/database/idRepoLamp_'+filtstr,os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'database','idRepoLamp_'+filtstr))
#                except IOError as e:
#                    print(e)
#                    print("ERROR: Cannot find already identified lines of this filter %s from Repository.."%(filt))
#                    print("Before you proceed to next step, do copy the identified line spectra of this filter.")
#                    print(" Or remove this image from "+os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'Spectrastoextract_Lamp_BandFilter.txt'))
#                    print('-'*10)
#
#            else:
#                print("More than 26 images (or no images!!) in single filter? Sorry, you should consider combining them somehow. Or split into two nights directory")
        #End of all filters of this night.
        outlog.close()
    print('All nights over...')

def MakeMasterFlat(PC, NormContdata, slit):
    """ This function will create a master flat using the continuum flat of each night and already generated master flat. This is
    to take care of noise for higher orders (orders 10, 11 and 12). Basically, we will use the master flat to remove noise 
    in higher orders and for the lower orders, the pipeline will use the continuum lamp observed in each night for the flat correction.
        In the end it will return the data of the new continuum flat."""
    EachNightNormContdata = NormContdata
    
    #getting data from master continuum flat
    #master continuum was created by average combine of normalised continuum flats observed on several nights.
    #finding the configuration file for spectrum extraction as well as ContinuumFile and ApertureLabel as references

    MasterContiName = 'master-cont1xd_' + slit + '.fits' 
    pkgpath = os.path.split(pkgutil.get_loader('pyTANSPEC').get_filename())[0]
    MasterContiData = fits.getdata(os.path.join(pkgpath,'data','CONTLAMPDIR',MasterContiName))
    
    #create the Mask..above this line EachNightNormContdata will be altered by MasterContiData
    ContinuumLOCData = np.load(os.path.join(pkgpath,'data','CONTLAMPDIR','ContinuumCutLine.npy'))
    Xvalue, Yvalue = ContinuumLOCData[:,0], ContinuumLOCData[:,1] 
    #fit the line...polynominal = 3
    z = np.polyfit(Xvalue, Yvalue, 3)
    p = np.poly1d(z) 

    Xnewvalue = np.arange(1, EachNightNormContdata.shape[0]+1, 1) 
    LOCarray = p(Xnewvalue)

    Ynewvalue= np.tile(np.arange(2048), (2048, 1))
    
    nrows, ncols = Ynewvalue.shape
    row, col = np.ogrid[:nrows, :ncols]
    boolmask = row < LOCarray
    
    NewFlatdata = np.where(boolmask, EachNightNormContdata, MasterContiData)
    
    return NewFlatdata
    
def DivideSmoothGradient(PC,inputimg,outputimg):
    """ This will divide the smooth gradient in image based on median filter smooting parameters in PC.
        In the end it will return the output filename."""
    hdulist = fits.open(inputimg)
    inputimgdata = hdulist[0].data
    inputimgdata = np.clip(inputimgdata,1, np.max(inputimgdata + 1))
    print('\n Generating median filtered gradient background using size : {0}'.format(PC.DVDMEDSMOOTHSIZE))
    print('It takes sometime (> 100 sec) to finish. Wait ...')
    try:
        smoothGrad = filters.median_filter(inputimgdata,size=PC.DVDMEDSMOOTHSIZE)
#        print(smoothGrad)
    except MemoryError:
        print('*** MEMORY ERROR : Skipping median filter Division ***')
        print('Try giving a smaller smooth size for median filter insted of {0}'.format(PC.DVDMEDSMOOTHSIZE))
        outputimg = inputimg  # Returning back the input file name since no subtraction was done.
    else:      
        prihdr = hdulist[0].header
        slit = prihdr['Slit'] #it is required when MakeMasterFlat will be called
        #Normalise this continuum flat using its median smoothed version
        NormContdata = inputimgdata / smoothGrad
        #generating the combined conti flat
        if PC.INSTRUMENT == 'TANSPEC' and PC.GRATINGHDR == 'GRATING1':
            hdulist[0].data = MakeMasterFlat(PC, NormContdata, slit) 
        else:
            hdulist[0].data = NormContdata      
        prihdr.add_history('Divided median filter Size:{0}'.format(PC.DVDMEDSMOOTHSIZE))
        hdulist.writeto(outputimg)
    finally:
        hdulist.close()
    # Return the name of the output filename
    return outputimg
                                
def RemoveCosmicRays(inputimg,outputimg):
    """ This will use La Cosmic algorithm to remove the cosmic rays in the image
        In the end it will return the output filename."""
    hdulist = fits.open(inputimg)
    inputimgdata = hdulist[0].data
    print('Removing cosmic rays..')
    crmask, cleanarr = astroscrappy.detect_cosmics(inputimgdata)
    prihdr= hdulist[0].header
    hdulist[0].data = cleanarr
    prihdr.add_history('Removed CosmicRays using LaCosmic')
    hdulist.writeto(outputimg)
    hdulist.close()
    # Return the name of the output filename
    return outputimg

def SubtractSmoothGradient(PC,inputimg,outputimg):
    """ This will subract smooth gradients in image based on median filter smooting parameters in PC.
        In the end it will return the output filename."""
    hdulist = fits.open(inputimg)
    inputimgdata = hdulist[0].data
    print('Calculating median filtered gradient background using size : {0}'.format(PC.MEDSMOOTHSIZE))
    try:
        smoothGrad = filters.median_filter(inputimgdata,size=PC.MEDSMOOTHSIZE)
    except MemoryError:
        print('*** MEMORY ERROR : Skipping median filter Subtraction ***')
        print('Try giving a smaller smooth size for median filter')
        outputimg = inputimg  # Returning mack the input file name since no subtraction was done.
    else:
        prihdr= hdulist[0].header
        hdulist[0].data = inputimgdata - smoothGrad
        prihdr.add_history('Subtracted median filter Size:{0}'.format(PC.MEDSMOOTHSIZE))
        hdulist.writeto(outputimg)
    finally:
        hdulist.close()
    # Return the name of the output filename
    return outputimg
                
#def FixBadPixels(PC,images,nightdir):
#    """ This will run iraf task proto.fixpix to interpolate badpixels """
#    if PC.TODO=='P' : PixelMask=nightdir+'/'+PC.PhotBadPixelMaskName
#    elif PC.TODO=='S' : PixelMask=nightdir+'/'+PC.SpecBadPixelMaskName
#    else : 
#        print("What are you doing? (S or P)")
#        return
#    if not os.path.isfile(PixelMask):
#        print("No Bad Pixel Mask file found by the name "+ PixelMask)
#        print("Hence skipping Bad pixel interpolation")
#        return

def CombDith_FlatCorr_subrout(PC,method="median"):
    """ This will combine (default=median) with avsigclip the images in single dither and also create corresponding normalized flats [,sky] and divide[,subtract] for flat [,sky] correction """
    directories = LoadDirectories(PC,CONF=False)
    for night in directories:
        PC.currentnight = night # Update the night directory for using GetFullPath()    
        print('Working on night: '+night)
        #Load all the Flat indexing file data
        with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects-FinalFlat.List'),'r') as FlatFILE :
#        with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects-Flat.List'),'r') as FlatFILE :     
            Flatfiledic=dict([(flatset.split()[0],flatset.rstrip().split()[1:]) for flatset in FlatFILE])
            
#        if SEPARATESKY=='Y' :
            #Load all the Sky files indexing file data
#            with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects-FinalSky.List'),'r') as SkyFILE :
#            with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects-Sky.List'),'r') as SkyFILE :
#                Skyfiledic=dict([(skyset.split()[0],skyset.rstrip().split()[1:]) for skyset in SkyFILE])  #Dictionary of Sky list for each image.

        #Load all the FilterSet indexing file data
        with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects.List'),'r') as FiltrFILE :
            Filtrfiledic=dict([(filtset.split()[0],shlex.split(filtset.rstrip())[1]) for filtset in FiltrFILE])  #Dictionary of filterset for each image.

        NewFiltSet='(Blah,Blah,Blah)'

        with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects2Combine.List'),'r') as Obj2CombFILE :
            #Secondly generate the list of lists of images to combine.
            ListofLists=[[]]
            for imgline in Obj2CombFILE:
                if len(imgline.split()) == 0 and len(ListofLists[-1]) != 0 :  ListofLists.append([])  #Start a new list at end
                elif len(imgline.split()) > 0 : ListofLists[-1].append(imgline.split()[0]) #Append to the last list

        if len(ListofLists[0]) == 0 : #No images this night..
            print('No images to work on this night. skipping...')
            try :
                os.remove(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects-ProcessedImg.List'))
            except OSError :
                print('Not able to remove (if any) previous '+os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects-ProcessedImg.List'))
            continue

        #Now iterate through every list of images to combine
        outlogFILE=open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects-ProcessedImg.List'),'w')
        file_used = []
        for imglist in ListofLists:
            if PC.TimeSeries == 'Y':
                print('Working on the set', imglist)
            if len(imglist) == 1 : #Single image. no need to combine
                OutCombimg=imglist[0]
                shutil.copy2(night+'/'+imglist[0],os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,OutCombimg))  #Copying the file dito..
            elif len(imglist) > 1 :
                OutCombimg=imglist[0][:-5]+'_'+method+'_'+imglist[-1][:-5]+'.fits'  #output file name

                with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,OutCombimg+'.imcombine.List'),'w') as imcombineinputFile:
                    imcombineinputFile.write('\n'.join([PC.RAWDATADIR+'/'+night+'/'+img for img in imglist])+'\n')
                    #imcombineinputFile.write('\n'.join(img for img in imglist)+'\n')
                    
                with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,OutCombimg+'.imcombine.List'),'r') as CombInpFile: 
                    CombfileList=[(combset.split()[0]) for combset in CombInpFile]
                    
                #Outputhdulist = combine_frames(os.path.join(PC.RAWDATADIR,night,CombfileList), method='median')
                Outputhdulist = imarith.combine_frames(CombfileList, method='median')
                imarith.WriteFitsOutput(Outputhdulist,os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,OutCombimg),overwrite=True)
                    
            #Now make list of flats to be combined for this image set
            Flats2Comb=[]
            for img in imglist:
                Flats2Comb+=Flatfiledic[img]  #Adding all the flat lists
            Flats2Comb=set(Flats2Comb)  #Making a set to remove duplicates
            #Write all these flat names to a file.
            imgflatlistfname=os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,OutCombimg[:-5]+'.flatlist')
            flatlist = open(imgflatlistfname, 'w')
            flatlist.close()
            flats2comblist = []
            # for fla in Flats2Comb:
            #     flats2comblist.append(fla)
            with open(imgflatlistfname,'a') as imgflatlistFILE :
                # imgflatlistFILE.write('\n'.join([PC.RAWDATADIR+'/'+night+'/'+fla for fla in Flats2Comb])+'\n')
                for fla in Flats2Comb:
                    imgflatlistFILE.write(os.path.join(PC.RAWDATADIR, PC.OUTDIR, night, fla)+'\n')

#            if SEPARATESKY=='Y' : #Now make list of skys to be combined for this image set
#                Skys2Comb=[]
#                for img in imglist:
#                    Skys2Comb+=Skyfiledic[img]  #Adding all the sky lists
#                Skys2Comb=set(Skys2Comb)  #Making a set to remove duplicates
                #Write all these Sky names to a file.
#                imgskylistfname=os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,OutCombimg[:-5]+'.skylist')
#                with open(imgskylistfname,'w') as imgskylistFILE :
#                    imgskylistFILE.write('\n'.join([night+'/'+sky for sky in Skys2Comb])+'\n')


###            print('Image section used for normalising Flat is '+FlatStatSection)
            outflatname=os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,OutCombimg[:-5]+'_flat.fits')
#            iraf.imcombine(input='@'+imgflatlistfname, output=outflatname, combine="median", scale="median",reject="sigclip", statsec=FlatStatSection)
#            Outputhdulist = combine_frames(imgflatlistfname,method='median')
#            imarith.WriteFitsOutput(Outputhdulist,outflatname,overwrite=False)

            with open(imgflatlistfname,'r') as imgflatlistFILE :
                    FlatCombfileList1=[(flatcombset.split()[0]) for flatcombset in imgflatlistFILE]
            FlatCombfileList = []
            for fname in FlatCombfileList1:
                path_file = os.path.split(fname)
                filename = path_file[-1]
                if filename[:9] == 'Biweight_':
                    FlatCombfileList.append(fname)
                #Outputhdulist = combine_frames(os.path.join(PC.RAWDATADIR,night,CombfileList), method='median')
            file_used.append(FlatCombfileList)
            Outputhdulist = imarith.combine_frames(FlatCombfileList, method='median')
            imarith.WriteFitsOutput(Outputhdulist,os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,outflatname),overwrite=True)
            
#            statout=iraf.imstatistics(outflatname+FlatStatSection,fields='mode',Stdout=1)
#            mode=float(statout[1])
            #We will normalise this flat with the mode of pixels in FlatStatSection
            Noutflatname=os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,OutCombimg[:-5]+'_Nflat.fits')
#            iraf.imarith(operand1=outflatname,op="/",operand2=mode,result=Noutflatname)
            #If we are subtracting sky, doing it before flat fielding.
#            if SEPARATESKY=='Y':
#                outskyname=os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,OutCombimg[:-5]+'_sky.fits')
#                iraf.imcombine (input='@'+imgskylistfname, output=outskyname, combine="median",reject="pclip")
                #Now subtract the sky form the science frame
#                OutSSimg=OutCombimg[:-5]+'_SS.fits'
#                iraf.imarith(operand1=os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,OutCombimg),op="-",operand2=outskyname,result=os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,OutSSimg))
 #               OutCombimg=OutSSimg  #Updating the new _SS appended input filename to continue as if nothing happened here.
            #Now divide by flat...
            if (PC.TODO == 'SX' or PC.TODO == 'SL') and (PC.CONTINUUMGRADREMOVE == 'Y'):
                # Flat corretion
                OutNormalisedFlat = DivideSmoothGradient(PC,outflatname,Noutflatname)                
           
            else:
                print('Edit the config file')
                print("Make sure that TODO= S and REMOVE_CONTINUUM_GRAD= Y")
                exit(1)
                #We will normalise this flat with the mode of pixels in FlatStatSection
                #Xmin, Xmax, Ymin, Ymax = parse_slicingstring(PC.FLATSTATSECTION)
                #mode = np.nanmedian(fits.getdata(PC.GetFullPath(outflatname))[Ymin-1:Ymax,Xmin-1:Xmax]) 
                #imarithiraf(operand1=outflatname,op="/",operand2=mode,result=Noutflatname)
                
                
            ###### Apply Flat correction to object frames
            OutFCimg=OutCombimg[:-5]+'_FC.fits'  
            Outputhdulist = imarith.divide_frames(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,OutCombimg), Noutflatname, FluxExt=[0])
            imarith.WriteFitsOutput(Outputhdulist,os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,OutFCimg),overwrite=True)
         #imarithiraf(operand1=os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,OutCombimg),op="/",operand2=Noutflatname,result=os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,OutFCimg))       
                                    
            
            #Now interpolate the bad pixels in the final image.
#            FixBadPixels(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,OutFCimg),night)

            ###### Now interpolate the bad pixels in the final image.
#            FixBadPixels(PC,PC.GetFullPath(OutFCimg),night)
            
            #If asked to do a smooth gradient removal in image, do it after everything is over now..
            if PC.GRADREMOVE == 'Y':
                print('Removing smooth Gradient from '+OutCombimg)
                OutGSobjimg = os.path.splitext(OutFCimg)[0]+'_GS.fits'
                #If Filter subtraciton was sucessfull it will return the output filename, else inputname
                OutGSobjimg = SubtractSmoothGradient(PC,PC.GetFullPath(OutFCobjimg),PC.GetFullPath(OutGSobjimg))
                #Updating the final image name with new _GS appended 
                OutFCimg = os.path.basename(OutGSobjimg)

            #If asked to do cosmic ray removal in the image after everything is over now..
            if PC.REMOVECOSMIC == 'Y':
                print('Removing Cosmic Rays from '+OutCombimg)
                OutCRobjimg = os.path.splitext(OutFCimg)[0]+'_CR.fits'
                OutCRobjimg = RemoveCosmicRays(PC.GetFullPath(OutFCimg),PC.GetFullPath(OutCRobjimg))
                #Updating the final image name with new _CR appended 00
                OutFCimg = os.path.basename(OutCRobjimg)


            if PC.TODO=='P':
                Oldfiltset=NewFiltSet
                NewFiltSet=Filtrfiledic[imglist[0]]
                if Oldfiltset != NewFiltSet : outlogFILE.write('\n')  #Entering a blank line to show change of filters in Photometry
#            elif TODO=='S':
#                outlogFILE.write('\n')  #Entering a blank line no matter what. We will ask user to change is they want to move and add.
            outlogFILE.write(imglist[0]+' '+OutFCimg+'\n')
        outlogFILE.close()
        if PC.TODO=='P': print('Edit the spaces (if required) between image sets in file '+os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'FirstoneANDcombinedImages.List')+' to align and combine them in next step.')
    print('All nights over...')        

def Manual_InspectCalframes_subrout(PC):
    """ This will display Flats, Sky and Lamps one after other, and based on user input select/reject """
    directories = LoadDirectories(PC,CONF=True)
    #filelist = ['AllObjects-Bias.List','AllObjects-Flat.List']
    filelist = ['AllObjects-Flat.List']
#    outfilelist = ['AllObjects-FinalBias.List','AllObjects-FinalFlat.List']
    outfilelist = ['AllObjects-FinalFlat.List']
    #if PC.USEALLFLATS == 'Y':
    #    filelist += ['AllFilter-Flat.List', 'AllFlat-Bias.List']
    #    outfilelist +=['AllFilter-FinalFlat.List', 'AllFlat-FinalBias.List']
    if PC.TODO == 'SX' or PC.TODO == 'SL' :
        filelist.append('AllObjects-Lamp.List')
        outfilelist.append('AllObjects-FinalLamp.List')
    #if PC.SEPARATESKY == 'Y' :
    #    filelist.append('AllObjects-Sky.List')
    #    outfilelist.append('AllObjects-FinalSky.List')

    AcceptAllThisNight = False  #Flags to skip this step
    AcceptAllEveryNight = False
    for night in directories:
        print("Working on night :\033[91m {0} \033[0m ".format(night))
        PC.currentnight = night # Upgating the night directory for using GetFullPath()
        AlwaysRemove = []
        AlwaysAccept = []

        for inpfile,outfile in zip(filelist,outfilelist):
            set_no = 0
            print('-*-'*8)
            print('Files in: \033[91m {0} \033[0m'.format(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,inpfile)))
            print('-*-'*8)
            inFILE = open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,inpfile),'r')
            ouFILE = open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,outfile),'w')
            print('-'*5)
            print('To skip all the remaining verifications, you can enter following two options')
            print('acceptall       # Accept all remaining images in this night')
            print('acceptallnights # Accept all images in all remaining nights')
            print('Use the above two options only when you are sure all the images are good. Do not use them if you are not sure.')
            print('-'*5)
            flag_list = []
            star_cal_dict = {}
            for inpline in inFILE:
                inplinelist = shlex.split(inpline.rstrip())
                if len(inplinelist) > 0 : ScienceImg = inplinelist[0]
                else : continue   #Skipping to next line
                CalImgs = [imgs for imgs in inplinelist[1:] if imgs not in AlwaysRemove]
                FinalCalImgs = CalImgs[:]
                print('For the science image: '+ScienceImg)
                if not AcceptAllThisNight :
                    for img in CalImgs:
                        if img not in AlwaysAccept:
                            #iraf.display(PC.GetFullPath(img),1)
                            #ImagePlot(fits.getdata(night+'/'+img), title=objline)
                            if inpfile == filelist[0] and PC.InspectFl == 'Y':
                                ImagePlot(fits.getdata(PC.GetFullPath(img)), title=img)
                            elif inpfile == filelist[1] and PC.InspectCal == 'Y':
                                ImagePlot(fits.getdata(PC.GetFullPath(img)), title=img)
                            print(PC.GetFullPath(img))
                            verdict=''
                            #verdict=input('Enter "r" to reject, "ra" to reject always in future, "aa" to always accept in future:')
                            verdict=input('Enter "r" to reject and "aa" to accept:')
                            if verdict == 'r' :
                                FinalCalImgs.remove(img)
                                print("Removing this image : "+img)
                            #elif verdict == 'ra' :
                            #    FinalCalImgs.remove(img)
                            #    AlwaysRemove.append(img)
                            #    print("Removing this image forever: "+img)
                            elif verdict == 'aa' :
                                AlwaysAccept.append(img)
                                print("Always accept this image forever this night : "+img)
                            elif verdict == 'acceptall' :
                                AcceptAllThisNight = True
                                print("Accepting every single remainging images of this night (Dangerous). ")
                                break
                            elif verdict == 'acceptallnights' :
                                AcceptAllEveryNight = True
                                AcceptAllThisNight = True
                                print("Accepting all remainging images from all remaining nights (Super Dangerous). ")
                                break

                if not FinalCalImgs : print('\033[91m ALERT: \033[0m No Calibration Images for {0} {1}'.format(night,ScienceImg))
                #Writing the final surviving calibration files to output file
                # Averaging have to be done here, then add file names as last 2 elements.
                input_path = night
                output_path = os.path.join(PC.RAWDATADIR,PC.OUTDIR,night)
                flag_list = FinalCalImgs
                f = False
                header_keys = ['SLIT', 'GRATING', 'OBJECT', 'CALMIR', 'ARGONL', 'NEONL', 'CONT1L','CONT2L']
#                print(star_cal_dict, FinalCalImgs, flag_list)
                for l in list(star_cal_dict.values()):
                    if FinalCalImgs == l[:-1] or FinalCalImgs == l[:-2]:
                        print('Given set already done.')
                        f = True
                        break
                if not f:
                    set_no += 1  
                    avglampfiles = imarith.combine_files(FinalCalImgs, input_path, output_path, keys=header_keys)
                star_cal_dict[ScienceImg] = flag_list
                FinalCalImgs.append(avglampfiles[0])
                if len(avglampfiles)>1:
                    FinalCalImgs.append(avglampfiles[1])
                ouFILE.write('{0} {1}\n'.format(ScienceImg,' '.join(FinalCalImgs)))
            ouFILE.close()
            inFILE.close()
        if not AcceptAllEveryNight : AcceptAllThisNight = False
    print('All nights over...') 

def DitherDetection(ObjectFile, ContWindowSelection, startLoc=None,avgHWindow=21,TraceHWidth=5):
    """identify the center of a spectrum window """
    if isinstance(ObjectFile  ,str):
        ObjectFile = fits.getdata(ObjectFile)

    if startLoc is None:
        startLoc = ObjectFile.shape[1]//2
    # Starting labelling Reference XD cut data; 
    WindowStart = ContWindowSelection[0] 
    WindowEnd = ContWindowSelection[1] 
    RefXD = np.nanmedian(ObjectFile[WindowStart:WindowEnd,startLoc-avgHWindow:startLoc+avgHWindow],axis=1) 
    Refpixels = np.arange(len(RefXD))+WindowStart
    Bkg = signal.order_filter(RefXD,domain=[True]*TraceHWidth*5,rank=int(TraceHWidth*5/10))
    Flux = np.abs(RefXD -Bkg)
    ThreshMask = RefXD > (Bkg + np.abs(mad_std(Flux))*6)
    centerpix = np.sum(Flux[ThreshMask]*Refpixels[ThreshMask])/np.sum(Flux[ThreshMask])
    
    return centerpix
 
def ImagePlot(imgdata, title=None):
    vmin = np.mean(imgdata) - 1.0*(np.std(imgdata))
    vmax = np.mean(imgdata) + 1.0*(np.std(imgdata))
    fig = plt.figure(figsize=(8,8))
    ax = plt.subplot()
    plt.imshow(imgdata, origin='lower', cmap='gray', vmin=vmin, vmax=vmax)
    plt.xlabel('Pixel column')
    plt.ylabel('Pixel row')
    
    if title:
        plt.title(title)
    fig.tight_layout()
    #plt.subplot_adjust(left=0.15)
    plt.colorbar(label='Counts')
    plt.show()     
                   

def Manual_InspectObj_subrout(PC):
    """ This will display one image after other, and based on user input classify images of each dither position """
    directories=LoadDirectories(PC,CONF=True)
#    if PC.TODO == 'P': print("Press _a_ and then _q_ over one good central star for selecting image")
#    if PC.TODO == 'S': print("Press _k_ and then _q_ over one good position on dispersed spectrum for selecting image \n IMP: Press k on some good part of star spectrum, not on the sky region around.")
    for night in directories:
        print("Working on night : "+night)
        ObjFILE = open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects.List'),'r')
        Obj2CombFILE = open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects2Combine.List'),'w')
        newX = 0
        newY = 0
        newsX = 0
        FWHM = 4.0  #Not important what number you give here...
        newfilter = 'Blah'
        newexptime = '-999'
        acceptall = False
        if PC.TimeSeries == 'Y':
            print('Since you need a timeseries data, I assume you do not want to inspect each frames')
            acceptall = True
        for objline in ObjFILE:
            try:
                img = objline.split()[0]
                if PC.INSPECTSC == 'Y' and acceptall == False:
                    ImagePlot(fits.getdata(night+'/'+img), title=objline)
            except IndexError:
                print('Blank line in '+os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects.List'))
                continue
            if acceptall:
                UserInput = 'aa'
            else:
                print('If you want to accept all the images without inspecting, enter "acceptall"')
                UserInput = input('Enter "r" to reject and "aa" to accept: ')
            if UserInput == 'r':
                print("Removing this image : "+img)
                continue
            elif UserInput == 'aa':
                print("Always accept this image forever this night : "+img)
                pass
            elif UserInput == 'acceptall':
                acceptall = True
                print("Accepting every single remainging images of this night (Dangerous). ")
            else:
                print('Enter either "r" or "aa", No other key works.')
                sys.exit()
                
           
            if PC.TODO == 'P': #If we were doing photometry (For Future)
                oldX = newX
                oldY = newY
                oldfilter = newfilter
                newfilter = shlex.split(objline)[1]
                oldexptime = newexptime
                newexptime = shlex.split(objline)[2]
#                FWHM = float(imx[3].split()[-1])
#                newX = float(imx[2].split()[0])
#                newY = float(imx[2].split()[1])
                #Print blank enter in the output file if the star has shifted
                StarShifted = np.sqrt((newX-oldX)**2 +(newY-oldY)**2) > 1*FWHM
            elif PC.TODO == 'SX' or PC.TODO == 'SL': #If doing spectroscopy
                oldsX = newsX
                oldfilter = newfilter
                newfilter = shlex.split(objline)[1]
                oldexptime = newexptime
                newexptime = shlex.split(objline)[3]
                FWHM = 4.0 # TODO try to find by python
                newsX = DitherDetection(os.path.join(PC.RAWDATADIR,night,img), PC.WindowSelection)
                #Print blank enter in the output file if the star has shifted
                StarShifted = np.abs(newsX-oldsX) > 1.5*FWHM #FWHM is defined as 4.
                
            # or filter wheels or exptime have changed.
            FiltersChanged = newfilter != oldfilter 
            ExptimeChanged = float(newexptime) != float(oldexptime)
            if StarShifted or FiltersChanged or ExptimeChanged or PC.TimeSeries == 'Y': Obj2CombFILE.write('\n')

            #Now, add this img name to dither image list
            Obj2CombFILE.write(img+' '+str(newX)+' '+str(newY)+'\n')
              
        Obj2CombFILE.close()
        ObjFILE.close()
        print('We have made the selected list of images in '+os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects2Combine.List')+' \n Add blank lines between file names to prevent them from median combining. \n Remove the blank line between file names, which you want to combine.')
        input("Press Enter to continue...")
        subprocess.call(PC.TEXTEDITOR.split()+[os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects2Combine.List')])

    print('All nights over...')
###########################################################################################################################################################################
def is_number(s):   # A function to check whether string s is a number or not.
    try:
        float(s)
        return True
    except ValueError:
        return False
        
def ListOfBiasesForGivenImg(PC,ImgIndices,BiasSizeDic):
    """ Returns the list of biase files which can be used for the input image subarray """
    # Create list of bias for this frame
    ImgXbegin, ImgXend,ImgYbegin, ImgYend, ImgXbin, ImgYbin = ImgIndices
    BiaslistForImg = []
    for key in BiasSizeDic:
        if (ImgXbegin, ImgXend,ImgYbegin, ImgYend, ImgXbin, ImgYbin) == key:
            BiaslistForImg += BiasSizeDic[key]
        elif (ImgXbegin >= key[0]) and (ImgXend <= key[1]) and (ImgYbegin >= key[2]) and (ImgYend <= key[3]) and (ImgXbin == key[4]) and (ImgYbin == key[5]):
            # We can trim out the bias frames to the image size
            slicesection = '[{0}:{1},{2}:{3}]'.format(ImgXbegin-key[0]+1,ImgXend-key[0]+1,ImgYbegin-key[2]+1,ImgYend-key[2]+1)
            for biasimg in BiasSizeDic[key]:
                outputtrimedbias = biasimg[:-5]+'_{0}-{1}_{2}-{3}.fits'.format(ImgXbegin,ImgXend,ImgYbegin,ImgYend)
                if not os.path.isfile(PC.GetFullPath(outputtrimedbias)):
                    iraf.imcopy(input=PC.GetFullPath(biasimg)+slicesection,output=PC.GetFullPath(outputtrimedbias))
                BiaslistForImg.append(outputtrimedbias)
    return BiaslistForImg
            
        
def SelectionofFrames_subrout(PC):
    """ Selects the images to reduce and create tables of corresponding Flat, Sky and Lamps """
    directories = LoadDirectories(PC,CONF=True)
    LogFilename = PC.NIGHTLOGFILE
    FiltREdic = dict()
    LampREdic = dict()
    LampNeREdic = dict()
    SkyREdic = dict()
    # Load Instrument
    Instrument = InstrumentObject(PC)
    
    # Filter/Grism Column index in log file
    if PC.TODO == 'P' :
        FiltColumn = 3
    elif PC.TODO[0] == 'S' :
        FiltColumn = 4

        
    print('*'*10) 
    print('For Regular Expression rules See: http://docs.python.org/2/howto/regex.html#regex-howto')
    print('Some examples of typical input are shown below')
    print(' .*M31.*   is the regular expression to select all the objects lines which has "M31" in it.')
    #print(' .*M31.*sky.*   is to select all the object lines which has both "M31" and then "sky" in it.')
    print('While selecting Lamp, Continuum etc, you can give a range of filenumbers to uniquely choose specific files')
    print(' .*conti1.* 89 121   is to select all filenames of _same_ filter which has "Continu" in it and has a filenumber in the range of 7 to 18')
    print('*'*10)
    ObjRE = ' '
    #Generating list of objects frames
    for night in directories:
        print("Working on night : "+night)
        PC.currentnight = night # Updating the night directory for using GetFullPath()
        print("Obs log file: file://{0}".format(os.path.join(PC.RAWDATADIR,night,LogFilename)))
        InpObjRE = input("Enter Regular Expression to select Science object frames (default: {0}): ".format(ObjRE)).strip(' ')
        if InpObjRE:
            ObjRE = InpObjRE
        regexpObj = re.compile(r''+ObjRE)

        with open(os.path.join(PC.RAWDATADIR,night,LogFilename),'r') as imglogFILE :
            # Skip blank lines and Commented out lines with #
            imglogFILElines = [imageLINE.rstrip() for imageLINE in imglogFILE if ((imageLINE.strip() != '') and (imageLINE[0] !='#'))]
            
        ObjList = [imgline.rstrip() for imgline in imglogFILElines if regexpObj.search(imgline.split()[0]) is not None ] 

        FiltList = set()  #Set to store the list of filters needed to find flat/Lamp for
        
        #Selecting Obj list
        with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects.List'),'w') as ObjFILE, open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjectsSky.List'),'w') as ObjFILEsky :
            for Objline in ObjList:
                if (PC.TODO=='P' and shlex.split(Objline)[6].upper() =="-NA-" and shlex.split(Objline)[7].upper() =="MIRROR") or (PC.TODO=='SX' and shlex.split(Objline)[4].upper() =="GRATING1" and shlex.split(Objline)[6].upper() !="-NA-" and shlex.split(Objline)[7].upper() !="MIRROR") or (PC.TODO=='SL' and shlex.split(Objline)[4].upper() =="GRATING2" and shlex.split(Objline)[6].upper() !="-NA-" and shlex.split(Objline)[7].upper() !="MIRROR") :
                    #continue
                    Name = shlex.split(Objline)[0]
                    FiltOrGrism = shlex.split(Objline)[FiltColumn]
                    Exptime = shlex.split(Objline)[2]
                    Slit = shlex.split(Objline)[7]
                    Date = shlex.split(Objline)[12]
                    #Time = shlex.split(Objline)[13]
                    FiltList.add(FiltOrGrism)
                    if (PC.TODO=='SX' or PC.TODO=='SL') and re.search('sky', Objline) :
                        ObjFILEsky.write('{0}  {1}  {2}  {3} {4}\n'.format(Name,FiltOrGrism,Slit, Exptime,Date))        
                    else:
                        ObjFILE.write('{0}  {1}  {2}  {3} {4} \n'.format(Name,FiltOrGrism,Slit, Exptime,Date))

        # Create a complete list of filters... we need to use flats from all nights
        CompleteFiltList = set()
        # print(FiltList) # Debugging
#        if PC.USEALLFLATS == 'Y':
#            for imgline in imglogFILElines:
#                ImgFrame = Instrument.IdentifyFrame(imgline)
#                if ImgFrame in ['BIAS','DARK', 'LAMP_SPEC']:  # Skip Bobvious Bias, Wavelength Calibration Lamps etc
#                    continue    #Skip these and go to the next object.
#                CompleteFiltList.add(shlex.split(imgline)[FiltColumn])


        if (not FiltList) and (not CompleteFiltList)  : #No files in this directory
            print('\033[91m ALERT: \033[0m No Images to reduce found in directory : {0}'.format(night))
            print('Please remove {0} from directory list next time.'.format(night))
            continue


        #Create dictionary for available bias frames or dark frames.
        BiasSizeDic = {}
        for imgline in imglogFILElines:
            if shlex.split(imgline)[1].upper() == 'BIAS':
                # Append to the list of biases stored in dictinary with the key (Xbegin,Xend,Ybegin,Yend,Xbin,Ybin)
                BiasSizeDic.append(imgline.split()[0])
                
        print(FiltOrGrism)        
                
        #Create dictionary of image sizes and available bias frames.
#        DarkSizeDic = {}
#        for imgline in imglogFILElines:
#            if shlex.split(imgline)[1].upper() == 'DARK':
#                # Append to the list of biases stored in dictinary with the key (Xbegin,Xend,Ybegin,Yend,Xbin,Ybin) 
#                DarkSizeDic.append(imgline.split()[0]) 

                     
#######################################################################MODIFY FOR PHOTOMETRY##############################################################################       
        #Now ask for flats in each filters
##############   Now ask for flats in each filters
##############   
        Flatlistdic = dict()
        Flatsizedic = dict()
        print("Below in addition to regexp, if needed you can enter the starting and ending filenumbers separated by space also.")
        for filt in FiltList|CompleteFiltList:
            filenumbregexp=re.compile(r'.*')
            if filt not in FiltREdic.keys() : 
                if PC.TODO=='P': FiltREdic[filt] = '.*[Ff]lat.*'  #Setting default to *[Ff]lat*
                elif PC.TODO[0]=='S': FiltREdic[filt] = '.*conti1.*'
            #Ask user again to confirm or change if he/she needs to
            InpfiltRE = input("Enter Regular Expression for the flat of filters %s (default: %s) : "%(str(filt),FiltREdic[filt])).strip(' ')
            if InpfiltRE :
                FiltREdic[filt] = InpfiltRE
                if len(InpfiltRE.split())==3 and is_number(InpfiltRE.split()[1]) and is_number(InpfiltRE.split()[2]):
                    filenumbregexp = re.compile('|'.join([str(i) for i in range(int(InpfiltRE.split()[1]), int(InpfiltRE.split()[2])+1)]))
                    FiltREdic[filt] = InpfiltRE.split()[0]
            regexpFilt = re.compile(r''+FiltREdic[filt])

            FlatList = []
            for imgline in imglogFILElines:
                ImgFrame = Instrument.IdentifyFrame(imgline)
                if ImgFrame in ['BIAS','LAMP_SPEC']:  # Skip Bobvious Bias, Wavelength Calibration Lamps etc
                    continue    #Skip these and go to the next object.
                if regexpFilt.search(' '.join(shlex.split(imgline)[0:2])) is not None :
                    if filt == shlex.split(imgline)[FiltColumn]: # Matching filter
                        if filenumbregexp.search(imgline.split()[-2]): # Second last column is filenumber
                            FlatList.append(imgline.split()[0])
 #                           Flatsizedic[imgline.split()[0]] = tuple([int(i) for i in shlex.split(imgline)[-7:-1]])  # (Xbegin,Xend,Ybegin,Yend,Xbin,Ybin) of each flat
            Flatlistdic[filt] = FlatList  #Saving flat list for this filter set

        #Now if Separate sky is being used to subtract, ask for each filter
#        if PC.SEPARATESKY == 'Y':
#            Skylistdic = dict()
#            print("Below in addition to regexp, if needed you can enter the starting and ending filenumbers separated by space also.")
#            for filt in FiltList:
#                filenumbregexp=re.compile(r'.*')
#                if filt not in SkyREdic.keys() : 
#                    SkyREdic[filt]=ObjRE+'_[Ss]ky.*'  #Setting default to object_sky
#                #Ask user again to confirm or change if he/she needs to
#                InpfiltRE=input("Enter Regular Expression for the Sky of filters %s (default: %s) : "%(str(filt),SkyREdic[filt])).strip(' ')
#                if InpfiltRE :
#                    SkyREdic[filt]=InpfiltRE
#                    if len(InpfiltRE.split())==3 and is_number(InpfiltRE.split()[1]) and is_number(InpfiltRE.split()[2]):
#                        filenumbregexp = re.compile('|'.join([str(i) for i in range(int(InpfiltRE.split()[1]), int(InpfiltRE.split()[2])+1)]))
#                        SkyREdic[filt] = InpfiltRE.split()[0]
#                        
#                regexpSky = re.compile(r''+SkyREdic[filt])
#                SkyList = []
#                for imgline in imglogFILElines:
#                    ImgFrame = Instrument.IdentifyFrame(imgline)
#                    if ImgFrame in ['BIAS','LAMP_SPEC']:  # Skip obvious Bias, Wavelength Calibration Lamps etc
#                        continue    #Skip these and go to the next object.
#                    if regexpSky.search(' '.join(shlex.split(imgline)[0:2])) is not None :
#                        if filt == shlex.split(imgline)[FiltColumn]: # Matching filter
#                            if filenumbregexp.search(imgline.split()[-1]): # Last column is filenumber
#                                SkyList.append(imgline.split()[0])
#                Skylistdic[filt] = SkyList  #Saving Sky list for this filter set

        #Now if We are doing Spectroscopy, Find the corresponding Lamp lamps also
        if PC.TODO[0] == 'S':
            Lamplistdic = dict()
            NeLamplistdic = dict()
            print("Below in addition, if needed you can enter the starting and ending filenumbers separated by space.")
            print(FiltList)
            for filt in FiltList:
                filenumbregexp = re.compile(r'.*')
                if filt not in LampREdic.keys() : 
                    LampREdic[filt] = '.*arlamp.*'  #Setting default to *Ar lamp*
                if filt not in LampNeREdic.keys() :                    
                    LampNeREdic[filt] = '.*nelamp.*'
                #Ask user again to confirm or change if he/she needs to
                InpfiltRE = input("Enter Regular Expression for the Ar-Lamp of filters %s (default: %s) : "%(str(filt),LampREdic[filt])).strip(' ')
                InpfiltNeRE = input("Enter Regular Expression for the Ne-Lamp of filters %s (default: %s) : "%(str(filt),LampNeREdic[filt])).strip(' ')
                if InpfiltRE and InpfiltNeRE:
                    LampREdic[filt] = InpfiltRE
                    LampNeREdic[filt] = InpfiltNeRE
                    if len(InpfiltRE.split())==3 and is_number(InpfiltRE.split()[1]) and is_number(InpfiltRE.split()[2]):
                        filenumbregexp = re.compile('|'.join([str(i) for i in range(int(InpfiltRE.split()[1]), int(InpfiltRE.split()[2])+1)]))
                    if len(InpfiltNeRE.split())==3 and is_number(InpfiltNeRE.split()[1]) and is_number(InpfiltNeRE.split()[2]):
                        Nefilenumbregexp = re.compile('|'.join([str(i) for i in range(int(InpfiltNeRE.split()[1]), int(InpfiltNeRE.split()[2])+1)]))
                        LampREdic[filt] = InpfiltRE.split()[0]
                        LampNeREdic[filt] = InpfiltNeRE.split()[0]
                                         
                regexpLamp = re.compile(r''+LampREdic[filt])
                regexpNeLamp = re.compile(r''+LampNeREdic[filt])                
                LampList = []
                LampListNe = []
                for imgline in imglogFILElines:
#                    ImgFrame = Instrument.IdentifyFrame(imgline)
#                    if ImgFrame != 'LAMP_SPEC':  # Skip if it is not a Wavelength Calibration Lamp.
#                        continue    #Skip these and go to the next object.
                    if regexpLamp.search(' '.join(shlex.split(imgline)[0:2])) is not None :
                        if filt == shlex.split(imgline)[FiltColumn]: # Matching filter
                            if filenumbregexp.search(imgline.split()[-2]): # 2nd last column is filenumber
                                LampList.append(imgline.split()[0])
                    if regexpNeLamp.search(' '.join(shlex.split(imgline)[0:2])) is not None :
                        if filt == shlex.split(imgline)[FiltColumn]: # Matching filter
                            if Nefilenumbregexp.search(imgline.split()[-2]): # 2nd last column is filenumber
                                LampListNe.append(imgline.split()[0])
                Lamplistdic[filt]=LampList  #Saving Ar-Lamp list for this filter set
                NeLamplistdic[filt]=LampListNe  #Saving Ne-Lamp list for this filter set                
            
        #Now, load the Object list and write to a file the Obj and corresponding flats/Lamps
        ObjFlatFILE = open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects-Flat.List'),'w')
#        ObjBiasFILE = open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects-Bias.List'),'w')
        if PC.TODO[0]=='S': ObjLampFILE = open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects-Lamp.List'),'w')
#        if PC.SEPARATESKY=='Y': ObjSkyFILE = open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects-Sky.List'),'w')

        with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects.List'),'r') as allobjectList :
            # Skip blank lines and Commented out lines with #
            ObjList=[imgline.rstrip() for imgline in allobjectList] #if regexpObj.search(imgline.split()[0]) is not None ]    #updated ObjList to be confirmed about object selection
        for Objline in ObjList:
            Name=shlex.split(Objline)[0]
            FiltOrGrism = shlex.split(Objline)[1]
            ObjFlatFILE.write(Name+' '+' '.join(Flatlistdic[FiltOrGrism])+'\n')
#
#            ObjimgIndices = tuple([int(i) for i in shlex.split(Objline)[-7:-1]])
#            ObjimgIndices = 0 #need to modify later
#            BiaslistForImg = ListOfBiasesForGivenImg(PC,ObjimgIndices,BiasSizeDic)

#            if BiaslistForImg:
#                ObjBiasFILE.write(Name+'  '+' '.join(BiaslistForImg)+'\n')
#            else:
#                print('\033[91m ERROR \033[0m: No Bias found for object img {0} with size {1}'.format(Name,ObjimgIndices))
#                if PC.SEPARATESKY=='Y':  print('Since you are subtracting seperate sky, I am ignoring this..')
#                else: print('Find bias for this night, copy here and rerun pipeline for this directory.')
#                ObjBiasFILE.write(Name+' \n') LamplistdicNe[filt]
                
            if PC.TODO[0]=='S' :ObjLampFILE.write(Name+'  '+' '.join(Lamplistdic[FiltOrGrism])+'  '+' '.join(NeLamplistdic[FiltOrGrism])+'\n')
#            if PC.SEPARATESKY=='Y': ObjSkyFILE.write(Name+'  '+' '.join(Skylistdic[FiltOrGrism])+'\n')
#        ObjFlatFILE.close()
#        ObjBiasFILE.close()
        # Create the files for all flat files
#        if PC.USEALLFLATS == 'Y':
#            FiltrFlatFILE = open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllFilter-Flat.List'),'w')
#            FlatBiasFILE = open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllFlat-Bias.List'),'w')
#            for filt in FiltList|CompleteFiltList:
#                FiltrFlatFILE.write('"{0}" {1}\n'.format(filt,' '.join(Flatlistdic[filt])))
#                for flatimg in Flatlistdic[filt]:
#                    FlatimgIndices = Flatsizedic[flatimg]
#                    BiaslistForFlat = ListOfBiasesForGivenImg(PC,FlatimgIndices,BiasSizeDic)

#                    if BiaslistForFlat:
#                        FlatBiasFILE.write('{0} {1}\n'.format(flatimg,' '.join(BiaslistForFlat)))
#                    else:
#                        print('\033[91m ERROR \033[0m: No Bias found for flat img {0} with size {1}'.format(flatimg,FlatimgIndices))
#                        print('Flat without bias subtraction is useless, and also harmful to other good flats of the night.')
#                        print('Find bias for this night, copy here and rerun pipeline. (or remove this flat!)')
#                        raise
#            FiltrFlatFILE.close()
#            FlatBiasFILE.close()
            
        print('Edit, save (if required) and close the Flat,Bias/Lamp/Sky list associations for this night :'+night)
        input("Press Enter to continue...")
        ObjFlatFILE.close()
        subprocess.call(PC.TEXTEDITOR.split()+[os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects-Flat.List')])
#        subprocess.call(PC.TEXTEDITOR.split()+[os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects-Bias.List')])
        if PC.TODO[0] == 'S': 
            ObjLampFILE.close()
            subprocess.call(PC.TEXTEDITOR.split()+[os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects-Lamp.List')])
        #if PC.SEPARATESKY == 'Y': 
        #    ObjSkyFILE.close()
        #    subprocess.call(PC.TEXTEDITOR.split()+[os.path.join(PC.RAWDATADIR,PC.OUTDIR,night,'AllObjects-Sky.List')])

    print('All nights over...') 

def CreateLogFilesFromFits_subrout(PC,hdu=0):
    """ Creates a log file of all the fits file in each directory 
        hdu = 0 is the index of hdu in hdulist to read in fits file"""

    directories = LoadDirectories(PC,CONF=True)
    LogFilename = PC.NIGHTLOGFILE
    LogColumns = [PC.OBJECTHDR, PC.EXPTIMEHDR, PC.FILTERHDR, PC.GRATINGHDR, PC.LAMPHDR, PC.GAINHDR, PC.SLITHDR, PC.ARGONLHDR, PC.NEONLHDR, PC.CONT1LHDR, PC.CONT2LHDR, PC.DATEHDR, PC.UTHDR, PC.RAHDR, PC.DECHDR, PC.COMMENTHDR, PC.FNUMHDR]
    
    #Final column heads will be: Filename LogColumns FileNumber

    RowString = ' "{'+'}" "{'.join(LogColumns)+'}"'  # String with  header keywords in curly brackets " "
    # Load Instrument
    Instrument = InstrumentObject(PC)

    for night in directories:
        print("Working on night : "+night)
        if os.path.isfile(os.path.join(PC.RAWDATADIR,night,LogFilename)):
            print('WARNING: Log file {0} already exist.'.format(os.path.join(PC.RAWDATADIR,night,LogFilename)))
            print('Skipping this directory now. If you want to create new log file, delete existing one or change filename NIGHTLOGFILE= in config file.')
            continue
        os.chdir(os.path.join(PC.RAWDATADIR,night))
        listOFimgs = [f for f in os.listdir('.') if re.match(r'^.*\.Z.fits$', f)]
        listOFimgs = sorted(listOFimgs) #, key=lambda k: int(k[:-7].split('-')[-1]))
        
        with open(LogFilename,'w') as outfile:
            for i,img in enumerate(listOFimgs):
                #hdulist = fits.open(img)
                #prihdr = hdulist[hdu].header
                prihdr = fits.getheader(img)            
                prihdr = Instrument.StandardiseHeader(prihdr,filename=img)
                for hkeys in LogColumns :   #Capture if these keywords are missing in header, and replace with -NA-
                    if hkeys not in prihdr : prihdr[hkeys]='-NA-'
                outfile.write(img+' '+RowString.format(**dict(prihdr))+' {0}\n'.format(i))
        os.chdir(PC.RAWDATADIR)
    print("{0} saved in each night's directory. Edit in manually for errors like ACTIVE filter.".format(LogFilename))


    try :
        directoriesF=open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,'directories'),'r')
    except IOError :
        #Creating a text file containing the directories which has SlopeimagesLog.txt logs to visit if it doesn't already exist
        directories=[dirs for dirs in os.walk(PC.RAWDATADIR).next()[1] if os.path.isfile(os.path.join(PC.RAWDATADIR,dirs,'SlopeimagesLog.txt'))]
        directories.sort()
        with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,'directories'),'w') as directoriesF : #Creating directories file
            directoriesF.write('\n'.join(directories)+'\n')
        #Now reopening the file to read and proceed
        directoriesF=open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,'directories'),'r')

    #Load directories list from the file
    directories = [dirs.rstrip().strip(' ').rstrip('/') for dirs in directoriesF if dirs.strip()] #Removing spaces or trailing / and ignore Blank empty lines
    directoriesF.close()

def LoadDirectories(PC,CONF=False):
    """ Loads the directories and return the list of directories to do analysis """
    try :
        directoriesF=open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,'directories'),'r')
    except IOError :
        #Creating a text file containing the directories which has fits files in it
        directories = [dirs for dirs in next(os.walk(PC.RAWDATADIR))[1] if glob.glob(os.path.join(PC.RAWDATADIR,dirs,'*.fits'))]
        directories.sort()
        with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,'directories'),'w') as directoriesF : #Creating directories file
            directoriesF.write('\n'.join(directories)+'\n')
        #Now reopening the file to read and proceed
        directoriesF=open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,'directories'),'r')

    #Load directories list from the file
    directories = [dirs.rstrip().strip(' ').rstrip('/') for dirs in directoriesF if dirs.strip()] #Removing spaces or trailing / and ignore Blank empty lines
    directoriesF.close()

    if CONF == True :
        #Ask user again to confirm or change if he/she needs to
        InpList=input('Enter the directories to analyse (default: %s) :'%', '.join(directories)).strip(' ')
        if InpList : 
            directories=[dirs.rstrip().strip(' ').rstrip('/') for dirs in InpList.split(',')] #Removing spaces or trailing /
            for dirs in list(directories): # iterating over a copy of the list
                if not os.path.isdir(os.path.join(PC.RAWDATADIR,dirs)):
                    print('Cannot find the data directory: {0} in the current directory {1}'.format(dirs,PC.RAWDATADIR))
                    print('WARNING: Removing the non-existing directory : {0} from the list'.format(dirs))
                    directories.remove(dirs)
                          
            with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,'directories'),'w') as directoriesF : #Updating directories file
                directoriesF.write('\n'.join(directories)+'\n')


    for dirs in directories:
        #Create a corresponding night directory in OUTPUT directory also if not already present.
        try:
            os.makedirs(os.path.join(PC.RAWDATADIR,PC.OUTDIR,dirs))
        except OSError:
            if os.path.isdir(os.path.join(PC.RAWDATADIR,PC.OUTDIR,dirs)) :
                print("Output directory "+os.path.join(PC.RAWDATADIR,PC.OUTDIR,dirs)+" already exists.\n Everything inside it will be overwritten.")
            else:
                raise
        else:
            print('Created directory :'+os.path.join(PC.RAWDATADIR,PC.OUTDIR,dirs))
        
    if len(directories) == 0 : 
        print('ERROR: No valid directories to reduce data found.')
        print('Atleast one directory containing data should be given as input.')
        sys.exit(1)
    else :
        return directories

def KeyboardInterrupt_handler():
    print('\nYou pressed Ctrl+C!')
    print('Stoping the reduction abruptly...')
    sys.exit(2)

def InitialTest(PC):
    """ Few initial tests to see all settings are correct . Input PC is the PipelineConfiguration object"""
            
    #Check for SExtractor installation
    if PC.TODO == 'P':
        with open('/dev/null','w') as devnull:
            try :
                subprocess.call(['sex','--version'],stdout=devnull)
            except OSError:
                try :
                    subprocess.call(['sextractor','--version'],stdout=devnull)
                except OSError:            
                    print('Cannot find the command: sex OR sextractor')
                    print('SExtractor needs to be installed from http://www.astromatic.net/software/sextractor')
                    raise

    #Check LampRepo directory exists.
#    if PC.TODO == 'S':
#    	 print("")
#        if not os.path.exists(PC.LAMPREPODIR):
#            print('Lamp Repository directory not found: {0}'.format(PC.LAMPREPODIR))
#            print('You can obtain {0} LampRepo directory by extracting data.tar.gz'.format(PC.INSTRUMENT))
#            print('Please add the correct path to Lamp Repository before proceeding.')
#            raise IOError(PC.LAMPREPODIR)


class PipelineConfig(object):
    """ This class is just a collection of variables required to run the pipeline """
    def __init__(self,ConfigFilename = None):
        if ConfigFilename is not None:
            self.LoadFromFile(ConfigFilename)

        self.RAWDATADIR = os.getcwd()
        # variables for helper functions
        self._FrameDatabaseInitialised = False
        self.currentnight = None

    def LoadFromFile(self,ConfigFilename):
        """ Loads the configuration from the input file """
        with open(ConfigFilename,'r') as configfile:
            for con in configfile:
                con = con.rstrip()
                if len(con.split()) >= 2 :
                    
                    if con.split()[0] == "INSTRUMENT=" :
                        self.INSTRUMENT = con.split()[1]
                    elif con.split()[0] == "TODO=" :
                        self.TODO = con.split()[1]
                    elif con.split()[0] == "OUTPUTDIR=" :
                        self.OUTDIR = con.split()[1]
                    #elif con.split()[0] == "VERBOSE=" :
                        #self.VER = con.split()[1]
                    elif con.split()[0] == "TEXTEDITOR=" :
                        self.TEXTEDITOR = shlex.split(con)[1]
                    #elif con.split()[0] == "OVERSCAN=" :
                    #    self.OVERSCAN = con.split()[1]
                    elif con.split()[0] == "Selected_Window_Dither_Find=" :
                        self.WindowSelection = (int(con.split()[1]), int(con.split()[2]))

                    elif con.split()[0] == "TIMESERIES=":
                        self.TimeSeries = con.split()[1][0].upper()
                    elif con.split()[0] == "INSPECTSC=":
                        self.INSPECTSC = con.split()[1][0].upper()
                    elif con.split()[0] == "INSPECTFL=":
                        self.InspectFl = con.split()[1][0].upper()
                    elif con.split()[0] == "INSPECTCAL=":
                        self.InspectCal = con.split()[1][0].upper()

                    elif con.split()[0] == "COMBINEIMGS=" :
                        self.IMGCOMBINE = con.split()[1][0].upper()
                    #elif con.split()[0] == "IMGCOMBMETHOD=" :
                   #     self.IMGCOMBMETHOD = con.split()[1]
                    #elif con.split()[0] == "FLATSTATSECTION=" :
                    #    self.FLATSTATSECTION = con.split()[1]
                    #elif con.split()[0] == "USEALLFLATS=" :
                    #    self.USEALLFLATS = con.split()[1][0].upper()

                    #elif con.split()[0] == "SEPARATE_SKY=" :
                    #    self.SEPARATESKY = con.split()[1][0].upper()
                    elif con.split()[0] == "REMOVE_COSMIC=" :
                        self.REMOVECOSMIC = con.split()[1][0].upper()

                    elif con.split()[0] == "GRADIENT_REMOVE=" :
                        self.GRADREMOVE = con.split()[1][0].upper()
                    elif con.split()[0] == "GRAD_FILT_SIZE=" :
                        self.MEDSMOOTHSIZE = (int(con.split()[1]), int(con.split()[2]))

                    #elif con.split()[0] == "BPMPHOTO=" :
                        #self.PhotBadPixelMaskName = con.split()[1]
                    #elif con.split()[0] == "BPMSPEC=" :
                        #self.SpecBadPixelMaskName = con.split()[1]
        
                    elif con.split()[0] == "THRESHOLD=" :
                        self.THRESHOLD = float(con.split()[1])
                    elif con.split()[0] == "EPADU=" :
                        self.EPADU = float(con.split()[1])
                    elif con.split()[0] == "READNOISE=" :
                        self.READNOISE = float(con.split()[1])
                    elif con.split()[0] == "DATAMAX=" :
                        self.DATAMAX = con.split()[1]
                    elif con.split()[0] == "XYMATCHMIN=" :
                        self.XYMATCHMIN = int(con.split()[1])
            
                    elif con.split()[0] == "APERTURE=" :
                        self.APERTURE = con.split()[1]
                    elif con.split()[0] == "ANNULUS=" :
                        self.ANNULUS = con.split()[1]
                    elif con.split()[0] == "DANNULUS=" :
                        self.DANNULUS = con.split()[1]

                    elif con.split()[0] == "EXPTIME=" :
                        self.EXPTIMEHDR = con.split()[1]
                    elif con.split()[0] == "FILTER=" :
                        self.FILTERHDR = con.split()[1]
                    elif con.split()[0] == "GRATING=" :
                        self.GRATINGHDR = con.split()[1]
                    elif con.split()[0] == "SLIT=" : 
                        self.SLITHDR = con.split()[1]
                    elif con.split()[0] == "CALMIR=" : 
                        self.LAMPHDR = con.split()[1]
                    elif con.split()[0] == "GAIN=" : 
                        self.GAINHDR = con.split()[1]          
                    elif con.split()[0] == "ARGONL=" : 
                        self.ARGONLHDR = con.split()[1]  
                    elif con.split()[0] == "NEONL=" : 
                        self.NEONLHDR = con.split()[1]
                    elif con.split()[0] == "CONT1L=" : 
                        self.CONT1LHDR = con.split()[1]
                    elif con.split()[0] == "CONT2L=" : 
                        self.CONT2LHDR = con.split()[1]
                            

                    elif con.split()[0] == "UT=" :
                        self.UTHDR = con.split()[1]
                    elif con.split()[0] == "ODATE=" :
                        self.DATEHDR = con.split()[1]

                    elif con.split()[0] == "OBJECT=" : 
                        self.OBJECTHDR = con.split()[1]
                    elif con.split()[0] == "COMMENT=" :
                        self.COMMENTHDR = con.split()[1]
                    elif con.split()[0] == "XSTART=" : 
                        self.XSTARTHDR = con.split()[1]
                    elif con.split()[0] == "XEND=" : 
                        self.XENDHDR = con.split()[1]
                    elif con.split()[0] == "YSTART=" : 
                        self.YSTARTHDR = con.split()[1]
                    elif con.split()[0] == "YEND=" : 
                        self.YENDHDR = con.split()[1]
                    elif con.split()[0] == "XBIN=" : 
                        self.XBINHDR = con.split()[1]
                    elif con.split()[0] == "YBIN=" : 
                        self.YBINHDR = con.split()[1]

                    elif con.split()[0] == "RA_HDR=" : 
                        self.RAHDR = con.split()[1]
                    elif con.split()[0] == "DEC_HDR=" :
                        self.DECHDR = con.split()[1]


                    elif con.split()[0] == "OUTPUT=" :
                        self.OUTPUTFILE = con.split()[1]
                    #elif con.split()[0] == "BACKUP=" :
                        #self.BACKUPDIR = con.split()[1]

                    #elif con.split()[0] == "CONVOLVEIMG=" :
                     #   self.CONVOLVEIMG = con.split()[1]
                    #elif con.split()[0] == "DOPSF=" :
                    #    self.DOPSF = con.split()[1]

                    #elif con.split()[0] == "LAMPDIRECTORY=" :
                    #    self.LAMPREPODIR = con.split()[1]
                    elif con.split()[0] == "SPECAPERTURE=" :
                        self.SPECAPERTURE = con.split()[1]
                    elif con.split()[0] == "SPECAPERTURE_LLIMIT=" :
                        self.SPECAPERTURE_LLIMIT = con.split()[1]
                    elif con.split()[0] == "SPECAPERTURE_ULIMIT=" :
                        self.SPECAPERTURE_ULIMIT = con.split()[1]


                    elif con.split()[0] == "BACKGROUND=" :
                        self.BACKGROUND = con.split()[1]
                    elif con.split()[0] == "TRACEFUNC=" :
                        self.TRACEFUNC = con.split()[1]
                    elif con.split()[0] == "TRACEORDER=" :
                        self.TRACEORDER = con.split()[1]
                    elif con.split()[0] == "NORMFUNC=" :  # Not used yet
                        self.NORMFUNC = con.split()[1]
                    elif con.split()[0] == "NORMORDER=" :  # Not used yet
                        self.NORMORDER = con.split()[1]
                    elif con.split()[0] == "SCOMBINE=" :
                        self.SCOMBINE = con.split()[1]
                    elif con.split()[0] == "DISPAXIS=" :
                        self.DISPAXIS = con.split()[1]
                    elif con.split()[0] == "REMOVE_CONTINUUM_GRAD=" :
                        self.CONTINUUMGRADREMOVE = con.split()[1][0].upper()
                    elif con.split()[0] == "CONT_GRAD_FILT_SIZE=" :
                        self.DVDMEDSMOOTHSIZE = (int(con.split()[1]), int(con.split()[2]))
                    elif con.split()[0] == "ApertureWindow=" :    
                        self.APERTUREWINDOW = (int(con.split()[1]), int(con.split()[2]))
                    elif con.split()[0] == "BkgWindows=" :
                        self.BKGWINDOWS = [(int(con.split()[1]),int(con.split()[2])), (int(con.split()[3]), int(con.split()[4]))]
                    elif con.split()[0] == "SPECCONFIGFILE=" :
                        self.SPECCONFIGFILE = con.split()[1]
                    elif con.split()[0] == "WLFitFunc=" :
                        self.WLFITFUNC = con.split()[1]                        
                    elif con.split()[0] == "RESPONSE_FUNCTION=":
                        self.RESPFN = con.split()[1]

                    elif con.split()[0] == "NIGHTLOGFILE=" :
                        self.NIGHTLOGFILE = con.split()[1]
                    elif con.split()[0] == "FNUM=" : # Varghese added this
                        self.FNUMHDR = con.split()[1]

    # Some helper functions for obtaining the full path to files.
    # They are not necessory for many steps.
    def GetFullPath(self,filename):
        """ Returns full path to filename.
        If the filename is a raw image, then it will give adress to the raw data directory.
        If the filename is new, then it will give address to the output directory of PC.currentnight""" 
        if not self._FrameDatabaseInitialised:
            #Initialise the database. # Don't call this if the Logfiles are not already created (in 0th step)
            self.FrameDatabase = {}
            for night in LoadDirectories(self,CONF=False):
                with open(os.path.join(self.RAWDATADIR,night,self.NIGHTLOGFILE),'r') as imglogFILE:
                    # Skip blank lines and Commented out lines with #
                    images = [imageLINE.split()[0] for imageLINE in imglogFILE if \
                              ((imageLINE.strip() != '') and (imageLINE[0] !='#'))]
                self.FrameDatabase[night] = set(images)
            self._FrameDatabaseInitialised = True

        if self.currentnight is None:
            # Raise error
            print('ERROR: Current directory not defined PC.currentnight')
            raise ValueError('ERROR: PC.currentnight not defined')

        # All fine now,  give full path intelegently.
        if filename in self.FrameDatabase[self.currentnight]:
            return os.path.join(self.RAWDATADIR,self.currentnight,filename)
        else:
            return os.path.join(self.RAWDATADIR,self.OUTDIR,self.currentnight,filename)
            

################################# Instrument Definitions  #################################################################
# Details as well as the name of each function calls of each instrument
# each step in pipeline is a function call with one input argument PC

# Ascii art courtesy : http://patorjk.com/software/taag/

class InstrumentObject(object):
    """ All the definitions and procedures of each instrument is defined inside this class."""
    def __init__(self,PC):
        self.PC = PC
        self.Name = self.PC.INSTRUMENT
        self.Banner = InstrumentDictionaries[self.Name]['Banner']
        self.About = InstrumentDictionaries[self.Name]['About']
        self.Steps = InstrumentDictionaries[self.Name]['Steps']
        self.StepsToRun = self.GetListofStepsToRun()

    def GetListofStepsToRun(self):
        """ Returns the list of steps to run based onPC, the PipelineConfiguration object """
        StepsToRun = []
        if self.Name in ['TANSPEC']:
#            StepsToRun += [0, 1, 2, 3, 4, 5, 11, 12]
            StepsToRun += [0, 1, 2, 3, 4, 5, 6, 7]
#            if self.PC.IMGCOMBINE == 'Y':
#                StepsToRun += [6]
#            if self.PC.TODO == 'P':
#                StepsToRun += [7, 8, 9, 10]
#            elif self.PC.TODO == 'S':
#                StepsToRun += [11, 12]

        else:
            print('Unknown Instrument')

        return StepsToRun

    def StandardiseHeader(self,prihdr,filename=None):
        """ Return the header object after standardising the values in it, specific to instrument."""
        
        if self.Name == 'TANSPEC':
            FNUM = prihdr['FNUM']
            prihdr['FNUM'] = int(FNUM)
            prihdr[self.PC.COMMENTHDR] = prihdr['COMMENT'][-1].split('/')[0][1:].strip()
            del prihdr['COMMENT']
            try:
                del prihdr['']
            except KeyError:
                pass
            if  prihdr[self.PC.GRATINGHDR] == '-1':
                prihdr[self.PC.GRATINGHDR]= 'grating1'
        
        elif self.Name == 'HFOSC':
            ut_sec = float(prihdr[self.PC.UTHDR])
            prihdr[self.PC.UTHDR] = str(datetime.timedelta(seconds=ut_sec))  # Converting to HH:MM:SS
            # XSTART, XEND, YSTART, YEND, XBIN, YBIN etc not in HFOSC header. so we shall add dummy values
            prihdr[self.PC.XSTARTHDR] = 1
            prihdr[self.PC.XENDHDR] = prihdr['NAXIS1']
            prihdr[self.PC.YSTARTHDR] = 1
            prihdr[self.PC.YENDHDR] = prihdr['NAXIS2']
            # If the Xsize is 250 then it must be most likely binned by 2 in Xaxis for HFOSC!!
            prihdr[self.PC.XBINHDR] = 2 if prihdr['NAXIS1'] == 250 else 1
            prihdr[self.PC.YBINHDR] = 1

        elif self.Name == 'IFOSC':
            prihdr[self.PC.EXPTIMEHDR] = float(prihdr[self.PC.EXPTIMEHDR])/1000.0  # Convert ms to sec for IFOSC

        elif self.Name == 'ARIES1.3m512':
            date_ut = prihdr[self.PC.UTHDR]  # Is in the form YYYY-MM-DDTHH:MM:SS
            date = date_ut.split('T')[0]
            ut = date_ut.split('T')[1]
            prihdr[self.PC.UTHDR] = ut  # in HH:MM:SS
            prihdr[self.PC.DATEHDR] = date
            
            # Filter and object name has to be read form filename
            fsearch = re.search('(.*)([ubvriUBVRI]).*',os.path.splitext(filename)[0])
            if fsearch:
                prihdr[self.PC.OBJECTHDR] = fsearch.group(1)
                prihdr[self.PC.FILTERHDR] = fsearch.group(2).upper()
            # Add the X start and Ystart header objects
            subimageformat = prihdr['SUBRECT']
            prihdr[self.PC.XSTARTHDR] = int(subimageformat.split(',')[0])
            prihdr[self.PC.XENDHDR] = int(subimageformat.split(',')[1])
            prihdr[self.PC.YSTARTHDR] = int(subimageformat.split(',')[2])
            prihdr[self.PC.YENDHDR] = int(subimageformat.split(',')[3])
        return prihdr

    def IdentifyFrame(self,ObjectLine):
        """ Returns what object the input frame is based on the ObjectLine 
        'Filename, PC.OBJECTHDR, PC.EXPTIMEHDR, PC.FILTERHDR, PC.GRATINGHDR, PC.LAMPHDR, PC.GAINHDR, PC.SLITHDR, PC.ARGONLHDR, PC.NEONLHDR, PC.CONT1LHDR, PC.CONT2LHDR, PC.DATEHDR, PC.UTHDR, PC.RAHDR, PC.DECHDR, PC.COMMENTHDR, PC.XSTARTHDR, PC.XENDHDR, PC.YSTARTHDR, PC.YENDHDR, PC.XBINHDR, PC.YBINHDR, FileNumber' """
        Frame = 'UNKNOWN'
        if self.Name == 'HFOSC':
            if float(shlex.split(ObjectLine)[2]) == 0:
                Frame = 'BIAS'
            elif 'GRISM' in shlex.split(ObjectLine)[4].upper():
                if 'FE-' in shlex.split(ObjectLine)[5].upper():
                    Frame = 'LAMP_SPEC'
                elif 'HALOGEN' in shlex.split(ObjectLine)[5].upper():
                    Frame = 'FLAT_SPEC'
                else:
                    Frame = 'OBJECT_SPEC'
            else:
                Frame = 'OBJECT_IMG'
            

        return Frame
        
        
                             

InstrumentDictionaries = {'TANSPEC':{
    'Banner' :"""

                     )  (    (                
  *   )    (      ( /(  )\ ) )\ )        (    
` )  /(    )\     )\())(()/((()/( (      )\   
 ( )(_))((((_)(  ((_)\  /(_))/(_)))\   (((_)  
(_(_())  )\ _ )\  _((_)(_)) (_)) ((_)  )\___  
|_   _|  (_)_\(_)| \| |/ __|| _ \| __|((/ __| 
  | |     / _ \  | .` |\__ \|  _/| _|  | (__  
  |_|    /_/ \_\ |_|\_||___/|_|  |___|  \___| 
                                              Data Reduction Pipeline...
""",
    'About' : "TANSPEC on 3.6-m DOT, ARIES ...",
    'Steps' : {
        0 : {'Menu': 'Generate Log files of fits files in each directory.',
             'RunMessage': "RUNNING TASK:0  Generating log files of fits files in each directory..",
             'function': CreateLogFilesFromFits_subrout },
        1 : {'Menu': 'Selection of object frames, Flats/Sky/Lamps etc. to reduce',
             'RunMessage': "RUNNING TASK:1 Selecting object frames, Flats/Sky/Lamps etc..",
             'function': SelectionofFrames_subrout },
        2 : {'Menu': 'Visually inspect and/or reject object images one by one.',
             'RunMessage': "RUNNING TASK:2  Visual inspection and/or rejection of object frames..",
             'function': Manual_InspectObj_subrout},
        3 : {#'Menu': 'Visually inspect and/or reject Bias/Flats/Sky/Lamps one by one.',
             'Menu': 'Visually inspect and/or reject Flats/Lamps one by one.',
             'RunMessage': "RUNNING TASK:3  Visual inspection and/or rejection of Flats/Lamps frames..",
             'function': Manual_InspectCalframes_subrout},
#        4 : {'Menu': 'Subtract Overscan,Bias/sky',
#             'RunMessage': "RUNNING TASK:4  Subtracting biases or sky..",
#             'function': Bias_Subtraction_subrout},
       4  : {#'Menu': 'Apply Flat Correction and/or Bad pixel interpolation, and/or CR removal',
             'Menu': 'Apply Flat Correction and/or CR removal',
             'RunMessage': "RUNNING TASK:4  Flat correction etc..",
             'function': CombDith_FlatCorr_subrout},
#        6 : {'Menu': 'Align and combine images.',
#             'RunMessage': "RUNNING TASK:6  Aligning and combining images..",
#             'function': AlignNcombine_subrout },
#        7 : {'Menu': 'Make the list of images, Images4Photo.in to do Photometry.',
#             'RunMessage': "RUNNING TASK:7  Makeing list of images, Images4Photo.in to do Photometry..",
#             'function': Createlist_subrout },
#        8 : {'Menu': 'Select Stars and Sky region of the field on first image',
#             'RunMessage': "RUNNING TASK:8  Selecting Stars and Sky region of the field from first image..",
#             'function': Star_sky_subrout },
#        9 : {'Menu': 'Create Sextracter config file & coordinate output of first image.',
#             'RunMessage': "RUNNING TASK:9  Create Sextracter config file & coordinate output of first image..",
#             'function': Sextractor_subrout },
#        10 : {'Menu': 'Do Photometry',
#             'RunMessage': "RUNNING TASK:10 Doing Photometry..",
#             'function': Photometry },
        5: {'Menu': 'Input Spectrum pair subtraction and filenames.',
             'RunMessage': "RUNNING TASK:5  Inputing Spectrum pair subtraction and filenames...",
             'function': SpectralPairSubtraction_subrout },
        6: {'Menu': 'Extract wavelength calibrated 1D spectra from image.',
             'RunMessage': "RUNNING TASK:6  Extracting wavelength calibrated 1D spectra..",
             'function': SpectralExtraction_subrout },
        7: {'Menu': 'Flux calibration of extracted 1D spectra.',
            'RunMessage': "RUNNING TASK:7 Flux calibration using Instrument Response..",
            'function': FluxCalibration_subrout }

    }
}
} # END of Instrument Definitions..#######
 
def parse_args(raw_args=None):
    """ Parses the command line input arguments """
    parser = argparse.ArgumentParser(description="Spectrum Extraction Pipeline")
    parser.add_argument('ConfigFile', type=str,
                         help="Configuration file which contains settings for TANSPEC XD-mode extraction. Hint: You can download sample config file from https://github.com/astrosupriyo/pyTANSPEC/blob/main/pyTANSPEC/config/TANSPECscript.config")
                         
    args = parser.parse_args(raw_args)  
    return args     
          
#-----Main Program Calls Begins here........................
# If you are using the functions seperately from module, see the Class PipelineConfig 
# to create the PipelineConfiguration instance which has to be passed on to many function calls.

def main(raw_args=None):
    """ The Pipeline main funtion to run while run from terminal as a standalone pipeline . """
    """ Extracts 2D spectrum image into 1D spectrum """
    args = parse_args(raw_args)
    
#    if len(sys.argv) < 2 :
#        print('-'*10)
#        print('Usage : {0} ScriptSettings.conf'.format(sys.argv[0]))
#        print('where,')
#        print('     ScriptSettings.conf is the configuration file for this run of reduction pipeline')
#        print(' ')
#        print("Note: This script should be run from the directory containing all the night's data directories.")
#        print('-'*10)
#        sys.exit(1)

    #pkgpath = os.path.split(pkgutil.get_loader('pyTANSPEC').get_filename())[0]
    #config_file_Main = os.path.join(pkgpath,'config/TANSPECscript.conf')
    
    try :
#        PC = PipelineConfig(ConfigFilename = sys.argv[1])
        PC = PipelineConfig(ConfigFilename = args.ConfigFile)
    except IOError :
        print ("Error: Cannot open the file "+args.ConfigFile+". Setup the config file based on ScriptSettings.conf file correctly, before running the script.")
        sys.exit(1)

    # Load Instrument
    Instrument = InstrumentObject(PC)
    # Print Welcome Banner.. 
    print("Time: {0}".format(time.strftime("%c")))
    print(Instrument.Banner)
    print(Instrument.About)

    print('Running Initial Tests...')
    try :
        InitialTest(PC)  # Some important Initial tests for Sanity before running the code.
    except Exception as e :
        print(e)
        sys.exit(1)

    print("\n *** Very Very Important: Backup your RAW data first. Don't proceed without backup. ***\n")
    input("  Hit Enter to continue...")
    print('-'*10)
    print('IMP: All outputs of this run will be written to the directory '+os.path.join(PC.RAWDATADIR,PC.OUTDIR)+'\n')
    try:
        os.makedirs(os.path.join(PC.RAWDATADIR,PC.OUTDIR))
    except OSError:
        if os.path.isdir(os.path.join(PC.RAWDATADIR,PC.OUTDIR)) :
            try:
                with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,'StepsFinished'),'r') as stepsoverFLS:
                    StepsOver = stepsoverFLS.read()
            except IOError:
                    StepsOver = 'Nothing...'
            print('Steps you have already finished : '+StepsOver)        
        
            print("\n WARNING : Output directory "+os.path.join(PC.RAWDATADIR,PC.OUTDIR)+" already exists.\n\n         If you run the steps that have already been finished, everything inside it will be overwritten. Be warned...")
            while True:
                UserInput = input('\n Wanna proceed (hit Y/N): ').strip(' ')
                if UserInput.upper() == 'Y':
                    break 
                elif UserInput.upper() == 'N':
                    sys.exit()
                else:
                    continue
        else:
            raise
    else:
        print('Created directory :'+os.path.join(PC.RAWDATADIR,PC.OUTDIR))
        
        
    if PC.TODO == 'P' : todoinwords = 'Photometry' # For future
    elif PC.TODO == 'SX' : todoinwords = 'XDSpectroscopy'
    elif PC.TODO == 'SL' : todoinwords = 'LRSpectroscopy'


    print(" ---------------- Welcome to \033[91m {0} {1} \033[0m Pipeline --------------- \n".format(PC.INSTRUMENT,todoinwords))
    print("Enter the Serial numbers (space separated if more than one task in succession)")
    for i in Instrument.StepsToRun:
        print("\n {0}  {1}".format(i,Instrument.Steps[i]['Menu']))
    print("--------------------------------------------------------------- \n")

    try:
        todo = input('Enter the list : ')
        todo = todo.split()

        for task in todo :
            print("Time now: {0}".format(time.strftime("%c")))
            CalledTheTask=True
            try:
                task = int(task)
                if task not in Instrument.StepsToRun : 
                    raise IndexError('Task asked to run is not consistent with .conf file') # Raise Index error if a task not supposed to be run was given.
                print(Instrument.Steps[task]['RunMessage'])
            except (ValueError, IndexError) :
                print('Cannot understand the input task: {0}'.format(task))
                print('Skipping task :{0}'.format(task))
                CalledTheTask=False
            else: # Run the task..
                Instrument.Steps[task]['function'](PC)

            if CalledTheTask :
                with open(os.path.join(PC.RAWDATADIR,PC.OUTDIR,'StepsFinished'),'a') as stepsoverFLS:
                    stepsoverFLS.write(str(task)+' ')
    except KeyboardInterrupt :
        KeyboardInterrupt_handler()

    print("All tasks over... Enjoy!!!")
    print('Reduced data are at '+os.path.join(PC.RAWDATADIR,PC.OUTDIR))       
if __name__ == "__main__":
    main()           

