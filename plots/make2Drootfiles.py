# make2Drootfiles.py



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import mplhep as hep
hep.style.use("CMS")
from coffea import util
import itertools
import os, sys
import glob
import copy
import uproot
import time

sys.path.append('../python/')
import functions

# functions.makeSaveDirectories()
# directory='../outputs/loosenottight/'
directory='../outputs/tight_jetid1_ak8pftrigger/'

label_map = functions.getLabelMap(directory)
label_to_int = {label: i for i, label in label_map.items()}
signal_cats = [ i for label, i in label_to_int.items() if '2t' in label]
pretag_cats = [ i for label, i in label_to_int.items() if 'pre' in label]
antitag_cats = [ i for label, i in label_to_int.items() if 'at' in label]

plt.rcParams["font.size"] = 20



toc = time.time()

if (len(sys.argv) > 1) and (sys.argv[1] in ['2016all', '2017', '2018', 'all']):
    
    year = sys.argv[1]

else:

    year = '2016all'

systematics = ['nominal', 'jes', 'jer', 'pileup', 'pdf', 'q2', 'btag', 'prefiring']
syst_labels = ['nominal']
if '2018' in year: 
    systematics = ['nominal', 'jes', 'jer', 'pileup', 'pdf', 'q2', 'btag']

for s in systematics:
    if not 'nominal' in s and not 'hem' in s:
        syst_labels.append(s+'Down')
        syst_labels.append(s+'Up')
        
print(syst_labels)


yearLabel = year.replace('20', '').replace('all','')

imagesavedir = 'images/tight_jetid1_ak8pftrigger/png/'

# tag = '_blind'
tag = ''

dataOnly = False
inclusive = False

if inclusive:
    
    cats, cat_labels = [''], ['']
    
else:

    cats = ['', 'cen', 'fwd', '0bcen', '0bfwd', '1bcen', '1bfwd', '2bcen', '2bfwd']
    cat_labels = ['', 'cen', 'fwd', 'cen0b', 'fwd0b', 'cen1b', 'fwd1b', 'cen2b', 'fwd2b']
    
#     cats = ['cen', 'fwd']
#     cat_labels = ['cen', 'fwd']

signals = []
signals = ['RSGluon2000', 'ZPrime2000_1', 'ZPrime2000_10', 'ZPrime2000_30', 'ZPrime2000_DM']
signals = ['RSGluon'+str(int(b*100)) for b in [10,15,20,25,30,35,40,45,50,55,60]]
signals += ['ZPrime'+str(int(b*100))+'_10' for b in [10,12,14,16,18,20,25,30,35,40,45,50,60,70]]
signals += ['ZPrime'+str(int(b*100))+'_30' for b in [10,12,14,16,18,20,25,30,35,40,45,50,60,70]]
signals += ['ZPrime'+str(int(b*100))+'_DM' for b in [10,15,20,25,30,35,40,45,50]]
signals += ['ZPrime'+str(int(b*100))+'_1' for b in [10,12,14,16,18,20,25,30,35,40,45]]


savefileheader = directory+'twodalphabet/TTbarAllHad{}_'.format(year.replace('20', '').replace('all',''))
                                                                
fdata  = uproot.recreate(savefileheader+'Data'+tag+'.root')

if not dataOnly:
    
    fttbar = uproot.recreate(savefileheader+'TTbar'+tag+'.root')
    sigfiles = [uproot.recreate(savefileheader+'signal'+sig+tag+'.root') for sig in signals ]

for cat, catname in zip(cats, cat_labels):
    
    if cat == '':
        
        signal_cats = [ i for label, i in label_to_int.items() if '2t' in label]
        antitag_cats = [ i for label, i in label_to_int.items() if 'at' in label]
        sum_axes = ['anacat']

    elif 'b' in cat :
        
        print('b cats')
        
        signal_cats = label_to_int['2t'+cat]
        antitag_cats = label_to_int['at'+cat]
        sum_axes = []
        
    else:
        
        
        signal_cats = []
        antitag_cats = []
        
        for label, i in label_to_int.items():
            
            if '2t' in label and cat in label:
                signal_cats.append(i)
            if 'at' in label and cat in label:
                antitag_cats.append(i)
                
        
        sum_axes = ['anacat']
        
        
    
    for syst in syst_labels:
        print(syst, cat)

        integrate_pass = {'anacat':signal_cats, 'systematic': syst}
        integrate_fail = {'anacat':antitag_cats, 'systematic': syst}

        systname = syst.upper()[:-2] + 'up' if 'Up' in syst else syst.upper()[:-4] + 'down'

        if 'nominal' in syst:
            
            print('getting files from ', directory+'/scale/')

            systname = ''
            hdata_pass = functions.getHist2('mtt_vs_mt', 'JetHT', year, sum_axes=sum_axes, integrate_axes=integrate_pass, tag=tag, coffea_dir=directory+'/scale/')
            hdata_fail = functions.getHist2('mtt_vs_mt', 'JetHT', year, sum_axes=sum_axes, integrate_axes=integrate_fail, tag=tag, coffea_dir=directory+'/scale/') 

            fdata["MttvsMt"+catname+yearLabel+"Pass"+systname] = hdata_pass
            fdata["MttvsMt"+catname+yearLabel+"Fail"+systname] = hdata_fail

            fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(nrows=2, ncols=2)
            hep.hist2dplot(hdata_pass, ax=ax3, label='hdata_pass_'+systname)
            hep.histplot(hdata_pass[{'jetmass':sum}], histtype='step', color='k', ax=ax1, label='hdata_pass_projy'+systname)
            hep.histplot(hdata_pass[{'ttbarmass':sum}], histtype='step', color='k', ax=ax2, label='hdata_pass_projx'+systname)
            plt.savefig(imagesavedir+'hdata_pass_'+catname+yearLabel+systname+'.png')
            plt.close()
            del fig, ax1
            
            fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(nrows=2, ncols=2)
            hep.hist2dplot(hdata_fail, color='k', ax=ax3, label='hdata_fail_'+systname)
            hep.histplot(hdata_fail[{'jetmass':sum}], histtype='step', color='k', ax=ax1, label='hdata_fail_projy'+systname)
            hep.histplot(hdata_fail[{'ttbarmass':sum}], histtype='step', color='k', ax=ax2, label='hdata_fail_projx'+systname)
            plt.savefig(imagesavedir+'hdata_fail_'+catname+yearLabel+systname+'.png')
            plt.close()
            del fig, ax1
            
        if not dataOnly:
            
            sig_pass = [functions.getHist2('mtt_vs_mt', sig, year, sum_axes=sum_axes, integrate_axes=integrate_pass, tag=tag, coffea_dir=directory+'/scale/') for sig in signals]
            sig_fail = [functions.getHist2('mtt_vs_mt', sig, year, sum_axes=sum_axes, integrate_axes=integrate_fail, tag=tag, coffea_dir=directory+'/scale/') for sig in signals]

            httbar_pass = functions.getHist2('mtt_vs_mt', 'TTbar', year, sum_axes=sum_axes, integrate_axes=integrate_pass, tag=tag, coffea_dir=directory+'/scale/') 
            httbar_fail = functions.getHist2('mtt_vs_mt', 'TTbar', year, sum_axes=sum_axes, integrate_axes=integrate_fail, tag=tag, coffea_dir=directory+'/scale/') 


            # save hists

            fttbar["MttvsMt"+catname+yearLabel+"Pass"+systname] = httbar_pass
            fttbar["MttvsMt"+catname+yearLabel+"Fail"+systname] = httbar_fail

            fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(nrows=2, ncols=2)
            hep.hist2dplot(httbar_pass, ax=ax3, label='httbar_pass_'+systname)
            hep.histplot(httbar_pass[{'jetmass':sum}], histtype='step', color='k', ax=ax1, label='httbar_pass_projy'+systname)
            hep.histplot(httbar_pass[{'ttbarmass':sum}], histtype='step', color='k', ax=ax2, label='httbar_pass_projx'+systname)
            plt.savefig(imagesavedir+'httbar_pass_'+catname+yearLabel+systname+'.png')
            plt.close()
            del fig, ax1
            
            fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(nrows=2, ncols=2)
            hep.hist2dplot(httbar_fail, color='k', ax=ax3, label='httbar_fail_'+systname)
            hep.histplot(httbar_fail[{'jetmass':sum}], histtype='step', color='k', ax=ax1, label='httbar_fail_projy'+systname)
            hep.histplot(httbar_fail[{'ttbarmass':sum}], histtype='step', color='k', ax=ax2, label='httbar_fail_projx'+systname)
            plt.savefig(imagesavedir+'httbar_fail_'+catname+yearLabel+systname+'.png')
            plt.close()
            del fig, ax1
            
            for i, file in enumerate(sigfiles):

                
                file["MttvsMt"+catname+yearLabel+"Pass"+systname] = sig_pass[i]
                file["MttvsMt"+catname+yearLabel+"Fail"+systname] = sig_fail[i]

                fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(nrows=2, ncols=2)
                hep.hist2dplot(sig_pass[i], ax=ax3, label=signals[i]+'_pass_'+systname)
                hep.histplot(sig_pass[i][{'jetmass':sum}], histtype='step', color='k', ax=ax1, label=signals[i]+'_pass_projy'+systname)
                hep.histplot(sig_pass[i][{'ttbarmass':sum}], histtype='step', color='k', ax=ax2, label=signals[i]+'_pass_projx'+systname)
                plt.savefig(imagesavedir+signals[i]+'_pass_'+catname+yearLabel+systname+'.png')
                plt.close()
                del fig, ax1
                
                fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(nrows=2, ncols=2)
                hep.hist2dplot(sig_fail[i], color='k', ax=ax3, label=signals[i]+'fail_'+systname)
                hep.histplot(sig_fail[i][{'jetmass':sum}], histtype='step', color='k', ax=ax1, label=signals[i]+'_fail_projy'+systname)
                hep.histplot(sig_fail[i][{'ttbarmass':sum}], histtype='step', color='k', ax=ax2, label=signals[i]+'_fail_projx'+systname)
                plt.savefig(imagesavedir+signals[i]+'_fail_'+catname+yearLabel+systname+'.png')
                plt.close()
                del fig, ax1
    
                    


fdata.close()
                                                                
print('saving '+savefileheader+'Data'+tag+'.root')
 
if not dataOnly:
    
    fttbar.close()
    for file in sigfiles:
        
        file.close()
        
    
    print('saving '+savefileheader+'TTbar'+tag+'.root')
    
    for sig in signals:
        print('saving '+savefileheader+sig+tag+'.root')

        
        
tic = time.time()
print()
functions.printTime(tic-toc)
