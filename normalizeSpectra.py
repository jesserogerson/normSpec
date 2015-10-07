#!/usr/bin/env python
'''
--------------------------------------------------------------------------------
Author: Jesse A. Rogerson, jesserogerson.com, rogerson@yorku.ca
Additional Credits: Patrick B. Hall,
                    Paola Rodriguez Hidalgo

These code takes a list of raw spectra from the command line and allows the user
to interactively normalize them. The program was specifically desinged to be
used for the main author's PhD thesis, and thus isn't necessarily general to
all raw spectra.

The program can broken down into roughly two separate functions:

1. normalize()
This function allows the user to interactively set the parameters of
normalization and then execute the normalization.

2. plotNorm()
This method is called automatically after the users has asked these code to
calculate the normalized spectra. plotNorm() allows the user to interactively
augment a plot of the normalized spectra in order to both check if the
normalization parameters were good/not, as well as have a publication-worthy
figure if desired.

After the users elects to quit plotNorm(), they are sent back to the original
normalize() function where they may augment the parameters and redo the
normalization, or quit normalizeSpectra.py altogether.

To Run:
----------

$> ./normalizeSpectra.py JHHMMSS.card

The *.card file is both the list of raw spectra and where the major Information
of the object is held. In order for these code to run the *.card file must
be structured in the following way:

The *.card file:
----------
JHHMMSS.card <--- MUST be called this
contents:
line0: SDSS Jhhmmss.ss+/-ddmmss.s
line1: RA Dec
line2: gmag
line3: redshift
line4: label1 MJD1 /path/to/raw/spec1/
line5: label1 MJD2 /path/to/raw/spec2/
line6: label1 MJD3 /path/to/raw/spec3/
line7: label1 MJD4 /path/to/raw/spec4/
                .
lineF: labelN MJDN /path/to/raw/specN/

The *.card file may have an arbitrary number of raw spectra for normalization
but the first four lines MUST be name/RA DEC/gmag/redshift.
The raw spectra lines MUST be label/MJD/path (space separated).

HISTORY
--------------------------------------------------------------------------------
2014-09-16 - JAR - created
2015-02-05 - JAR - changed 'data' to accept dictionary
                 - made plotting part of the program
2015-02-10 - JAR - removed plotting, doesn't need to be in the 'normalize'
                 - function
2015-04-17 - JAR - converted to a general program, for use by ALL
2015-04-20 - JAR - added 'plotNorm()' method, which makes a normalized
                 - spectral plot
2015-04-22 - JAR - added objInfo{} dictionary... to be passed around
                 - merged mjd{} dicionary into objInfo{}
2015-08-26 - JAR - changed the output normalized ascii spectra to have
                   suffix's equal to their name (SDSS, GEM1, etc.)
                 - validation for yscaling added, in case there is no
                   spectrum labeled 'SDSS' given to the program
                 - created a clearer separation between normalization
                   routine and the plotting routine, which both have
                   commmand pages.
                 - filled out the docstring above.
2015-09-04 - JAR - original commit to github.com/jesserogerson
2015-09-15 - JAR - multiple bug fixes, especially around parameter files
                 - added date/time to each writing of the param files
2015-09-16 - JAR - multiple bug fixes, more user-entry validation
                 - in plotNorm, made legend lines thicker
                 - added smoothing to plotNorm
2015-09-21 - JAR - addd SNRreg to normJHHMMSS.parm file
                 - changed defaults of lw, smooth, xlimits
                 - made parmfiles automatically written when quit or normalize
2015-09-24 - JAR - fixed the *.card read-in to account for the newly added gmag
                 - bug fix: lw was reading in as str() not float()
2015-09-25 - JAR - added an SNR output file, which takes all the SNRs calculated
                   and outputs them, for use later.
2015-09-28 - JAR - changed SNRreg default from [1600,1650] to [1600,1700]
2015-09-29 - JAR - changed defaults for various passed parameters
                 - changed the SNR output file's format
2015-09-30 - JAR - added a print statement so user knows SNRreg from param file
2015-10-01 - JAR - added dashed line at 1.0 (continuum) for normalized plot
                   (also added minor tick marks to this plot on y-axis)
                 - bug fix: adding/removing spectra from plot in plotNorm()
                 - added SNRreg warning, is the SNRreg inside all lambda cover?
--------------------------------------------------------------------------------
'''
#Libraries used
import scipy.optimize as spot
import numpy as np
from sys import argv
import sys
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import copy as cp
import os.path
import argparse
import jarTools
import datetime

####declare methods() and functions()
def plotNorm(spectra,
             normList,
             RLF,
             colourDict,
             objInfo,
             smooth=True,
             xlimits=[1200,1600],
             ylimits=[0,2.5],
             lw=1.0,
             annotations=False):
    '''
    Plotting Program, build a normalized spectra plot from the plotlist
    '''
    yes=set(['yes','y','YES','Y',True,1])
    no=set(['no','n','NO','N',False,0])

    normOrig=cp.deepcopy(spectra)
    #Search for normJHHMMSS.parm file?
    parmFile='plot'+objInfo['shortObjName']+'.parm'
    if os.path.exists(parmFile):
        print '############################################################'
        print '*** Detected a plotting parameter file:',parmFile
        user_input=raw_input('*** Would you like to use it? [y/n]:')
        if user_input in yes:
            parmDict={}
            print '*** Reading in previously used parameters'
            print '***'
            with open(parmFile,'r') as f:
                s=f.readlines()
                for line in s[-6:-1]:
                    listedline=line.strip().split('=')
                    parmDict[listedline[0]]=listedline[1]
            f.close()
            #pulling out the parm values
            if parmDict['annotations']=='True':
                annotations=True
            else:
                annotations=False
            lw=float(parmDict['lw'])
            xlimits=map(float,parmDict['xlimits'].split(','))
            ylimits=map(float,parmDict['ylimits'].split(','))
            temp=map(float,parmDict['RLF'].split(','))
            RLF=[[temp[0],temp[1]]]
            for t in range(2,len(temp)-1,2):
                RLF.insert(0,[temp[t],temp[t+1]])
            print '*** xlimits'+'='+str(xlimits)
            print '*** ylimits'+'='+str(ylimits)
            print '*** RLF'+'='+str(RLF)
            print '*** annotations'+'='+str(annotations)
            print '*** lw'+'='+str(lw)
        else:
            print '*** Using default parameters.'
        print '############################################################'

    filename='norm'+objInfo['shortObjName']+'.eps'
    print '#######----------------------------------------------#######'
    print '#######---------Normalized Spectra Plotter-----------#######'
    print '#######----------------------------------------------#######'
    print '### I made a plot for you. See-->',filename
    print '### be sure to refresh to see your changes take effect.'
    #plotList=cp.deepcopy(spectra.keys()) #the list that will be plotted
    if smooth==True:
        print '*** smoothing spectrum'
        for spec in spectra:
            spectra[spec][:,1]=np.array(jarTools.boxcarSmooth(spectra[spec][:,1]))
    escape=False
    first=False
    windows=False
    annotations=True
    user_input='commands'
    while escape==False:
        #build a plot to play with
        fig = plt.figure()
        ax1=fig.add_subplot(111)
        ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
        plt.rc('text',usetex=True)
        plt.rc('font',family='sans-serif')
        plt.plot([100,10000],[1.0,1.0],'--',color='k')

        plt.xlim(xlimits[0],xlimits[1])
        plt.ylim(ylimits[0],ylimits[1])

        #sort plotList by smallest to largest MJD
        plotList=sorted(normList, key=objInfo.get)

        #plot all normalized spectra in plotlist
        #calculate rest-frame time between observations on the fly
        for i,spec in enumerate(plotList):
            if i==0:
                deltaT=0
            else:
                deltaT=round((objInfo[plotList[i]]-objInfo[plotList[i-1]])/(1+objInfo['zem']),2)
            plt.plot(spectra[spec][:,0],(spectra[spec][:,1]),colourDict[spec],linewidth=lw,label=str(round(objInfo[spec],2))+' '+spec+' '+str(deltaT))

        #turns on/off the RLF gray'd out regions
        if windows==True:
            for w in RLF:
                plt.axvspan(w[0],w[1],facecolor='0.8',linewidth=0)

        #Setting labels, ticks, limits on y-axis and bottom x-axis
        ax1.set_xlabel('Rest-frame Wavelength (\AA)')
        ax1.set_ylabel('Normalized Flux Density (10$^{-17}$ erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$)')
        ax1.set_xlim(xlimits[0],xlimits[1])
        ax1.set_xticks([1250,1350,1450,1550,1650])
        ax1.set_ylim(ylimits[0],ylimits[1])
        ax1.xaxis.set_minor_locator(MultipleLocator(25))

        #The 2nd axis (which is really just the top x-axis
        ax2=ax1.twiny() #copies everything from the y
        ax2.set_xlim(xlimits[0]*(1+objInfo['zem']),xlimits[1]*(1+objInfo['zem'])) #set the observed frame
        ax2.set_xticks([4000,4500,5000,5500]) #the ticks I want
        ax2.set_xlabel('Observed-frame Wavelength (\AA)')
        ax2.xaxis.set_minor_locator(MultipleLocator(100))

        if annotations==True:
            #Adding Annotations
            ax1.annotate(objInfo['objName'],xy=(1275,(ylimits[1]*0.95)))
            ax1.annotate('z='+str(objInfo['zem']),xy=(1450,(ylimits[1]*0.95)))
            #ax1.annotate('CIV',xy=(1542,1.73))
            #ax1.annotate('SIV',xy=(1392,1.73))
            #ax1.plot([1550,1550],[1.6,1.7],'k',linewidth=1)
            #ax1.plot([1400,1400],[1.6,1.7],'k',linewidth=1)
            #Adding the legend
            leg=ax1.legend(loc='lower left',prop={'size':12})
            for legobj in leg.legendHandles:
                legobj.set_linewidth(2.5)
        plt.savefig(filename)

        #let it play the first 'options' command first
        if first==True:
            user_input=raw_input('Enter a command:')
        first=True
        if user_input=='commands':
            print '############################################################'
            print 'You can modify the plot of your normalized spectra'
            print 'using the commands below. Be sure to refersh the'
            print 'plot '+filename+' everytime you make a change:'
            print 'commands       : displays list of all command options'
            print 'q,Q            : to quit the EW measurement'
            print 'xlimits        : create a new xrange by entering [x1,x2]'
            print 'ylimits        : create a new yrange by entering [y1,y2]'
            print 'RLF            : turn plotting for these regions on/off'
            print 'filename       : change name of image file'
            print 'plotlist       : add or remove spectra from final plot.'
            print 'annotations    : turn on/off annotations/legend/etc'
            print 'lw             : change linewidth for plotted spectra'
            print 'smooth         : smooth the spectra.'
            print '############################################################'
        elif user_input=='smooth':
            print '############################################################'
            user_input=raw_input('Turn on smoothing? [y,n]:')
            if user_input in yes:
                smooth=True
                for spec in spectra:
                    spectra[spec][:,1]=np.array(jarTools.boxcarSmooth(spectra[spec][:,1]))
            elif user_input in no:
                smooth=False
                spectra=cp.deepcopy(normOrig)
            else:
                print user_input+': Not a valid entry. Back to command page.'
            print 'Smoothing:'+str(smooth)
            print '############################################################'
        elif user_input=='lw':
            print '############################################################'
            try:
                print 'Current linewidth:',lw
                user_input=raw_input('Enter new linewidth (float):')
                lw=float(user_input)
                print 'New linewidth:',lw
            except ValueError:
                print 'That didnt make any sense!'
            print '############################################################'
        elif user_input=='xlimits':
            print '############################################################'
            try:
                user_input=raw_input('Enter new x-axis limits (comma separated):')
                xlimits=map(float,user_input.split(','))
                print 'Reset figure xlimits to:'+str(xlimits)
            except ValueError:
                print 'Nope, no good, maybe you typed something wrong?'
            print '############################################################'
        elif user_input=='ylimits':
            print '############################################################'
            try:
                user_input=raw_input('Enter new y-axis limits (comma separated):')
                ylimits=map(float,user_input.split(','))
                print 'Reset figure ylimits to:'+str(ylimits)
            except ValueError:
                print 'Not a real set of numbers... somethings wrong'
            print '############################################################'
        elif user_input=='annotations':
            print '############################################################'
            print 'Annotations: object name, redshift, RA, Dec, legend, etc.'
            user_input=raw_input('Plot the annotations? [y,n]:')
            if user_input in yes:
                annotations=True
            elif user_input in no:
                annotations=False
            else:
                print user_input+': Not a valid entry. Back to command page.'
            print 'Plotting all annotations',annotations
            print '############################################################'
        elif user_input=='q' or user_input=='Q':
            print '############################################################'
            escape=True
            print 'Quitting'
            print ''
            print 'Writing current plotting parameters to:',parmFile
            print 'xlimits'+'='+str(xlimits)
            print 'ylimits'+'='+str(ylimits)
            print 'RLF'+'='+str(RLF)
            print 'annotations'+'='+str(annotations)
            print 'lw'+'='+str(lw)
            now = datetime.datetime.now()
            outfile=open(parmFile,'a')
            outfile.write('-------------------------'+now.strftime("%Y-%m-%d %H:%M")+'------------------------\n')
            outfile.write('annotations'+'='+str(annotations)+'\n')
            outfile.write('lw'+'='+str(lw)+'\n')
            outfile.write('xlimits'+'='+str(xlimits[0])+','+str(xlimits[1])+'\n')
            outfile.write('ylimits'+'='+str(ylimits[0])+','+str(ylimits[1])+'\n')
            outfile.write('RLF'+'=')
            for r in RLF[:-1]:
                outfile.write(str(r[0])+','+str(r[1])+',')
            else:
                outfile.write(str(RLF[-1][0])+','+str(RLF[-1][1])+'\n')
            outfile.write('----------------------------------------------------------------------------\n')
            outfile.close()
            print '############################################################'
        elif user_input=='filename':
            print '############################################################'
            print 'Current output filename:',filename
            user_input=raw_input('Enter a filename ( .eps will be added to end):')
            filename=user_input+'.eps'
            print 'Reset output filename to:',filename
            print '############################################################'
        elif user_input=='plotlist':
            print '############################################################'
            print 'Current list of spectra to plot:',plotList
            print 'To add OR remove, enter the name.'
            user_input=raw_input('Enter name of spectra:')
            if user_input in plotList:
                normList.remove(user_input)
            elif user_input not in plotList:
                normList.append(user_input)
            print 'New list of spectra to plot:',plotList
            print '############################################################'
        elif user_input=='RLF':
            print '############################################################'
            user_input=raw_input('Plot the RLF windows? [y,n]:')
            if user_input in yes:
                windows=True
            elif user_input in no:
                windows=False
            else:
                print user_input+': Not a valid entry. Back to command page.'
            print 'Plotting Relatively Line Free Windows:'+str(windows)
            print '############################################################'
        else:
            print '############################################################'
            print 'What chu talkin bout Willis'
            print '############################################################'
    print 'Plotted normalized spectra in:','norm'+objInfo['shortObjName']+'.eps'
    print '#######------------------EXITING---------------------#######'
    print '#######---------Normalized Spectra Plotter-----------#######'
    print '#######----------------------------------------------#######'

def normalize(spectra,
              objInfo,
              smooth=True,
              funcType='plaw',
              RLF=[[1300,1320],[1590,1620],[1700,1750]],
              xlimits=[1100,1800],
              ylimits=[0,40],
              SNRreg=[1600,1700]):
    '''
    Normalization Routine
    '''
    #For validation of responses coming up
    yes=set(['yes','y','YES','Y',True,1])
    no=set(['no','n','NO','N',False,0])

    #Search for spectraJHHMMSS.parm file?
    parmFile='norm'+objInfo['shortObjName']+'.parm'
    if os.path.exists(parmFile):
        print '------------------------------------------------------------'
        print '*** Detected a normalization parameter file:',parmFile
        user_input=raw_input('*** Would you like to use it? [y/n]:')
        if user_input in yes:
            print '*** Reading in previously used parameters...'
            print '*'
            parmDict={}
            with open(parmFile,'r') as f:
                s=f.readlines()
                for line in s[-7:-1]:
                    listedline=line.strip().split('=')
                    parmDict[listedline[0]]=listedline[1]
            f.close()
            #pulling out the parm values
            if parmDict['smooth']=='True':
                smooth=True
            else:
                smooth=False
            funcType=str(parmDict['funcType'])
            SNRreg=map(float,parmDict['SNRreg'].split(','))
            xlimits=map(float,parmDict['xlimits'].split(','))
            ylimits=map(float,parmDict['ylimits'].split(','))
            temp=map(float,parmDict['RLF'].split(','))
            RLF=[[temp[0],temp[1]]]
            for t in range(2,len(temp)-1,2):
                RLF.insert(0,[temp[t],temp[t+1]])
            print '*** SNRreg'+'='+str(SNRreg)
            print '*** smooth'+'='+str(smooth)
            print '*** funcType'+'='+str(funcType)
            print '*** xlimits'+'='+str(xlimits)
            print '*** ylimits'+'='+str(ylimits)
            print '*** RLF'+'='+str(RLF)
        else:
            print '*** Using default parameters.'
        print '------------------------------------------------------------'
    #Constants
    lightspeed=299792.458 #km/s
    civ_0=1548.202 #Ang

    filename='spectra'+objInfo['shortObjName']+'.eps'
    #Pull out keys of the incoming dictionaries
    spectraOriginal=cp.deepcopy(spectra) #keeping a real copy of the original
    spectraNormalized={} #will be populated by the normalized data arrays
    keyList=spectra.keys() #just to have a keylist, cause why not
    normList=cp.deepcopy(keyList) #the list that will be normalized/plotted
    colourDict={'SDSS':'k','SDSS1':'k','SDSS2':'0.70',
    'BOSS':'r','BOSS1':'r','BOSS2':'b',
    'GEM':'c','GEM1':'c','GEM2':'g','GEM3':'orange'}
    #J022143
    #colourDict={'SDSS1':'b','SDSS2':'g','SDSS3':'r','SDSS4':'c',
    #'SDSS5':'m','SDSS6':'y','BOSS':'r','BOSS1':'r','BOSS2':'b',
    #'GEM':'c','GEM1':'c','GEM2':'g','GEM3':'orange'}
    #J073232, J083546, J083017
    #colourDict={'SDSS':'k','SDSS1':'k','SDSS2':'0.70','BOSS':'r',
    #'GEM1':'c','GEM2':'g','GEM3':'orange','GEM4':'m','GEM5':'b'}
    #J015017
    #colourDict={'SDSS1':'k','SDSS2':'0.70','SDSS3':'g',
    #'BOSS1':'r','BOSS2':'b',
    #'GEM':'c',}



    #YSCALE - this is a dictionary that is used to automatically scale
    #the y-axis fluxes to be near eachother. (scaled to SDSS value)
    #there is some validation done here that allows the user to Choose
    #which spectrum they want to scale to if an SDSS spectrum is not
    #sent in (just makes it a bit more general)
    print '*** scaling raw spectra to match SDSS'
    print '*** scaling using the mean flux value between 1270 < lambda < 1350'
    scaleToName='SDSS'
    validate=False
    if scaleToName not in keyList:
        print '-----ASIDE:'
        print '-----No SDSS spectrum!'
        while validate==False:
            print '-----which spectrum would you like to scale to?'
            print '-----'+str(keyList)
            scaleToName=raw_input('-----enter spectrum:')
            if scaleToName not in keyList:
                print '-----That spectrum is not available, try again'
            else:
                print '*** scaling raw spectra to match '+scaleToName
                validate=True
    yscale={}
    for spec in spectra:
        #find the region specified by 1270->1350
        w=[index for index,value in enumerate(spectra[spec][:,0]) if value > 1590 and value < 1650]
        #find the mean flux value in that region
        yscale[spec]=np.mean(spectra[spec][w,1])
    for spec in yscale:
        #calculate the ratio between the flux of any given spectrum
        #and the scaleToName spectrum, this ratio will be the
        #value you we scale the unnormalized spectra by for easy plotting
        if spec==scaleToName:
            continue
        yscale[spec]=yscale[scaleToName]/yscale[spec]
    yscale[scaleToName]=1.0
    #while loop only escapes when asked
    #'first' is designed to make the useability easier the while loop first
    #        plots a spectrum, THEN ask the user for input. But 'first' allows
    #        it to execute the initially programmed command of
    #        user_input='commands' so the user knows
    if smooth==True:
        print '*** smoothing spectrum'
        for spec in spectra:
            spectra[spec][:,1]=np.array(jarTools.boxcarSmooth(spectra[spec][:,1]))
    escape=False
    first=False
    user_input='commands'
    while escape==False:
        if first==False:
            print '------------------------------------------------------------'
            print '-----------------------Normalizer---------------------------'
            print '------------------------------------------------------------'
            print 'The following spectra have been detected...'
            print 'labels:'+str(keyList)
            print '*** Plot built, see '+filename
        plt.clf()
        plt.rc('text',usetex=True)
        plt.rc('font',family='sans-serif')
        plt.xlim(xlimits[0],xlimits[1])
        plt.ylim(ylimits[0],ylimits[1])
        for spec in normList:
            plt.plot(spectra[spec][:,0],(spectra[spec][:,1]*yscale[spec]),color=colourDict[spec])
        for spec in normList:
            if spec not in spectraNormalized:
                continue
            plt.plot(spectra[spec][:,0],yscale[spec]*(spectra[spec][:,1]/spectraNormalized[spec][:,1]),color=colourDict[spec],linestyle='--')
        for w in RLF:
            plt.axvspan(w[0],w[1],facecolor='0.9',linewidth=0)
        plt.xlabel('Rest-frame Wavelength (\AA)')
        plt.ylabel('Flux Density (10$^{-17}$ erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$)')
        plt.savefig(filename,transparent=False)

        #SNRregion validation - do all spectra have coverage for this SNRreg?
        for spec in normList:
            if max(spectra[spec][:,0])<SNRreg[1]:
                print '*** WARNING The SNR region youve selected doesnt'
                print '    is not fully covered by'+spec

        #see above - skips the first question
        if first==True:
            user_input=raw_input('Enter a command:')
        first=True

        #What to do with each user_command given
        if user_input=='smooth':
            print '------------------------------------------------------------'
            user_input=raw_input('Turn on smoothing? [y,n]:')
            if user_input in yes:
                smooth=True
                for spec in spectra:
                    spectra[spec][:,1]=np.array(jarTools.boxcarSmooth(spectra[spec][:,1]))
            elif user_input in no:
                smooth=False
                spectra=cp.deepcopy(spectraOriginal)
            else:
                print user_input+': Not a valid entry. Back to command page.'
            print 'Smoothing:'+str(smooth)
            print '------------------------------------------------------------'
        elif user_input=='xlimits':
            print '------------------------------------------------------------'
            try:
                user_input=raw_input('Enter new x-axis limits (comma separated):')
                xlimits=map(float,user_input.split(','))
                print 'Reset figure xlimits to:'+str(xlimits)
            except ValueError:
                print 'That aint valid, yo..'
            print '------------------------------------------------------------'
        elif user_input=='ylimits':
            print '------------------------------------------------------------'
            try:
                user_input=raw_input('Enter new y-axis limits (comma separated):')
                ylimits=map(float,user_input.split(','))
                print 'Reset figure ylimits to:'+str(ylimits)
            except ValueError:
                print 'try not sucking at typing... that might work...'
            print '------------------------------------------------------------'
        elif user_input=='SNRreg':
            print '------------------------------------------------------------'
            print 'Current region to calculate SNR over:',SNRreg
            try:
                user_input=raw_input('Enter new region to calculate SNR (comma separated):')
                SNRreg=map(float,user_input.split(','))
                print 'Reset figure SNRreg to:'+str(SNRreg)
            except ValueError:
                print 'Try again mo fo... that didnt work'
            print '------------------------------------------------------------'
        elif user_input=='q' or user_input=='Q':
            escape=True
            print '------------------------------------------------------------'
            print 'Writing current normalization parameters to:',parmFile
            print 'SNRreg'+'='+str(SNRreg)
            print 'smooth'+'='+str(smooth)
            print 'funcType'+'='+str(funcType)
            print 'xlimits'+'='+str(xlimits)
            print 'ylimits'+'='+str(ylimits)
            print 'RLF'+'='+str(RLF)
            now = datetime.datetime.now()
            outfile=open(parmFile,'a')
            outfile.write('-------------------------'+now.strftime("%Y-%m-%d %H:%M")+'-------------------------\n')
            outfile.write('SNRreg'+'='+str(SNRreg[0])+','+str(SNRreg[1])+'\n')
            outfile.write('smooth'+'='+str(smooth)+'\n')
            outfile.write('funcType'+'='+str(funcType)+'\n')
            outfile.write('xlimits'+'='+str(xlimits[0])+','+str(xlimits[1])+'\n')
            outfile.write('ylimits'+'='+str(ylimits[0])+','+str(ylimits[1])+'\n')
            outfile.write('RLF=')
            for r in RLF[:-1]:
                outfile.write(str(r[0])+','+str(r[1])+',')
            else:
                outfile.write(str(RLF[-1][0])+','+str(RLF[-1][1])+'\n')
            outfile.write('----------------------------------------------------------------------------\n')
            outfile.close()
            print '------------------------------------------------------------'
        elif user_input=='commands':
            print '------------------------------------------------------------'
            print 'You can change the parameters which your spectra'
            print 'are normalized by here.'
            print 'Choose one of the commands below, and be sure to refresh'
            print 'the plot '+filename+' to see the changes.'
            print ''
            print 'commands       : displays list of all command options'
            print 'q,Q            : to quit the EW measurement'
            print 'smooth         : smooth the continuum [y,n]'
            print 'funcType       : choose from plaw or poly'
            print 'xlimits        : create a new xrange by entering [x1,x2]'
            print 'ylimits        : create a new yrange by entering [y1,y2]'
            print 'RLF            : add/remove RLF windows[x1,x2]'
            print 'SNRreg         : change the region SNR is calculated over'
            print 'filename       : change name of image file'
            print 'normalize      : execute normalization.'
            print 'normlist       : add or remove spectra from final plot.'
            print '------------------------------------------------------------'
        elif user_input=='filename':
            print '------------------------------------------------------------'
            try:
                user_input=raw_input('Enter a filename ( .eps will be added to end):')
                filename=user_input+'.eps'
                print 'Reset output filename to:'+str(filename)
            except ValueError:
                print 'wha happen?'
            print '------------------------------------------------------------'
        elif user_input=='normlist':
            print '------------------------------------------------------------'
            print 'Current list of spectra to plot:',normList
            print 'To add OR remove, enter the name.'
            user_input=raw_input('Enter name of spectra:')
            if user_input in normList:
                normList.remove(user_input)
            elif user_input not in normList:
                normList.append(user_input)
            print 'New list of spectra to plot:',normList
            print '------------------------------------------------------------'
        elif user_input=='RLF':
            print '------------------------------------------------------------'
            print 'Current RLF windows:',RLF
            print 'To add OR remove, enter the windows beginning/ending.'
            try:
                user_input=raw_input('Enter a RLF window (comma separated):')
                temp=map(float,user_input.split(','))
                w=[round(temp[0],0),round(temp[1],0)]
                found=False
                for i,bounds in enumerate(RLF):
                    if w==bounds:
                        found=True
                        if len(RLF)==1:
                            print 'Sorry, you cannot have an empty RLF array'
                        elif len(RLF)==2:
                            RLF.pop(i)
                            print 'Removed:',w
                            print '[WARNING]: Your RLF array now has only one entry.'
                            print '[WARNING]: You are about to do bad science.'
                        elif len(RLF)>=3:
                            RLF.pop(i)
                            print 'Removed:',w
                if found==False:
                    print 'Added:',w
                    RLF.append(w)
                w,temp=0,0
                print 'New RLF windows:',RLF
            except ValueError:
                print 'That didnt make any sense, back to command page.'
            print '------------------------------------------------------------'
        elif user_input=='normalize':
            print '------------------------------------------------------------'
            print 'Writing current normalization parameters to:',parmFile
            print 'SNRreg'+'='+str(SNRreg)
            print 'smooth'+'='+str(smooth)
            print 'funcType'+'='+str(funcType)
            print 'xlimits'+'='+str(xlimits)
            print 'ylimits'+'='+str(ylimits)
            print 'RLF'+'='+str(RLF)
            now = datetime.datetime.now()
            outfile=open(parmFile,'a')
            outfile.write('-------------------------'+now.strftime("%Y-%m-%d %H:%M")+'-------------------------\n')
            outfile.write('SNRreg'+'='+str(SNRreg[0])+','+str(SNRreg[1])+'\n')
            outfile.write('smooth'+'='+str(smooth)+'\n')
            outfile.write('funcType'+'='+str(funcType)+'\n')
            outfile.write('xlimits'+'='+str(xlimits[0])+','+str(xlimits[1])+'\n')
            outfile.write('ylimits'+'='+str(ylimits[0])+','+str(ylimits[1])+'\n')
            outfile.write('RLF=')
            for r in RLF[:-1]:
                outfile.write(str(r[0])+','+str(r[1])+',')
            else:
                outfile.write(str(RLF[-1][0])+','+str(RLF[-1][1])+'\n')
            outfile.write('----------------------------------------------------------------------------\n')
            outfile.close()
            print '------------------------------------------------------------'
            print '***Normalizing the following spectra:'
            print normList
            print '(If all the spectra are not in the list above, it is because'
            print 'you took some out of the normlist)'
            SNRoutput=objInfo['objName'][6:]+' '+str(SNRreg[0])+' '+str(SNRreg[1])
            for spec in normList:
                data,normalized,lam,flux,flux_err=[],[],0,0,0
                #validate: make sure the datacube is shape 3
                #must move to next one if not
                if np.shape(spectra[spec])[1]!=3:
                    print '-----ASIDE:'
                    print '-----Data Array associated with label -'+spec+'- is INCORRECT shape.'
                    print '-----Requires 3 columns: lambda,flux,flux_err.'
                    print '-----...Exiting entire program'
                    sys.exit()
                #starting analysis
                data=spectra[spec]
                normalized=np.zeros(np.shape(data)) #numpy return array
                lam=data[:,0]
                flux=data[:,1]
                flux_err=data[:,2]
                print '----------------------------------------------------'
                print '***Normalizing spectrum: '+spec
                SNR=0
                SNR=np.array([])
                for i in range(len(spectraOriginal[spec][:,0])):
                    if spectraOriginal[spec][i,0] >= SNRreg[0] and spectraOriginal[spec][i,0] <= SNRreg[1]:
                        SNR=np.concatenate((SNR,([spectraOriginal[spec][i,1]/spectraOriginal[spec][i,2]])))
                print '*** SNR in range '+str(SNRreg)+'is '+str(np.median(SNR))
                #prepping for SNR writeout
                SNRoutput=SNRoutput+' '+spec+' '+str(np.median(SNR))
                #identify the indicies that reflect the given RLF windows
                w=0
                w=np.array([],dtype=int)
                for bounds in RLF:
                    temp_w=[index for index,value in enumerate(lam) if value > bounds[0] and value < bounds[1]]
                    w=np.concatenate((w,temp_w))
                    temp_w=[]
                print '*** Windows used for function fitting:'
                print RLF
                #Choose Continuum fitting function
                if funcType=='poly':
                    #NOTE: There is a BUILT-IN polyfit function in numpy
                    print '*** Normalizing using a Polynomial Fit'
                    print '*** Fitting Function to data: y = mx + b'
                    fit=np.polyfit(lam[w],flux[w],1)
                    yfit=fit[1]+(fit[0]*lam)
                    print '*** Solution Found: y = ('+str(fit[0])+')x + ('+str(fit[1])+')'
                elif funcType=='plaw':
                    #NOTE: was required to BUILD MY OWN power-law function
                    #it is defined in 'jarTools.powerlaw()'
                    print '*** Normalizing using a Power-law Fit'
                    print '*** Fitting function to data: y = a*x^b'
                    fit=jarTools.powerfit(lam[w],flux[w],flux_err[w])
                    yfit=fit[1]*lam**fit[0]
                    print '*** Solution Found: y = ('+str(fit[1])+')x^('+str(fit[0])+')'
                else:
                    print '*** Do not recognize specified fitting function'
                normalized[:,0]=lam
                normalized[:,1]=cp.deepcopy(spectraOriginal[spec][:,1])/yfit
                normalized[:,2]=cp.deepcopy(spectraOriginal[spec][:,2])/yfit
                spectraNormalized[spec]=normalized
                outfile=open(normFileList[spec],'w')
                for i in range(len(normalized)):
                    outfile.write(str(normalized[i,0]*(1+objInfo['zem']))+' '+
                                  str(normalized[i,1])+' '+
                                  str(normalized[i,2])+'\n')
                outfile.close()
                print '*** Spectrum Normalized: '+spec
                print '*** Written to file:',normFileList[spec]
                print '*** NOTE: normalized the UNsmoothed spectrum'
                print '*** Finished with: '+spec
            print '------------------------------------------------------------'
            print '*** All spectra are normalized'
            print '*** writing Signal-to-Noise ratios to file...'
            outfile=open('SNR_outfile.dat','a')
            outfile.write(SNRoutput+'\n')
            outfile.close
            print '*** calling plotting program'
            plotNorm(spectraNormalized,normList,RLF,colourDict,objInfo)
            user_input='commands'
            first=False
    print '-------------------------EXITING----------------------------'
    print '------------------------Normalizer--------------------------'
    print '------------------------------------------------------------'
    print '-----------------------..EXITING..--------------------------'
    print '-------------------------Program----------------------------'
    return spectraNormalized
#
#--------------------------------------------------------------------------#
#
#Main program begins here, calls the above functions

#arguments from the command line
script,filename=argv
print '----------------------------------------------------'
print '***Working on:',filename
if filename[-4:] !='card':
    print 'File must be a *.card file containing:'
    print 'SDSS Jhhmmss.ss+/-ddmmss'
    print 'RA Dec'
    print 'redshift'
    print 'name & MJD & location of SDSS spectrum'
    print 'name & MJD & location of BOSS spectrum'
    print 'name & MJD & location of GEM spectrum'
    print '***EXITING'
    sys.exit()

#read in contents of filename
f=open(filename,'r')
lines=[line.rstrip('\n') for line in f]
f.close()

objInfo={}
objInfo['objName']=lines[0]
objInfo['shortObjName']=filename[-12:-5]
coords=lines[1].split()
objInfo['RA']=float(coords[0])
objInfo['Dec']=float(coords[1])
objInfo['gmag']=float(lines[2])
objInfo['zem']=float(lines[3])

spectra={}
normFileList={}
#run a loop from 4th line to end of lines
for l in lines[4:]:
    temp=l.split()
    key=temp[0] ### spectrum name must be FIRST!
    spectra[key]=np.genfromtxt(temp[2],usecols=(0,1,2))
    normFileList[key]='norm'+objInfo['shortObjName']+'.'+key.lower()
    objInfo[key]=float(temp[1])

print 'Information in card file:'
print 'objName:',objInfo['objName']
print 'redshift:',objInfo['zem']
print 'g_mag:',objInfo['gmag']
print 'RA:',objInfo['RA'],'Dec:',objInfo['Dec']
print 'Spectra:',spectra.keys()
print 'HK: scaling to rest-frame.'

for spec in spectra:
    spectra[spec][:,0]=spectra[spec][:,0]/(1.+objInfo['zem'])
print '*** heading into normalization routine, follow commands to normalize.'
print '----------'

normspec=normalize(spectra,objInfo)
