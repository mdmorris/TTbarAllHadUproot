# scaleCoffeaFiles.py


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



label_map = functions.getLabelMap(coffea_dir='../outputs/plot_3d/')
label_to_int = {label: i for i, label in label_map.items()}
signal_cats = [ i for label, i in label_to_int.items() if '2t' in label]
pretag_cats = [ i for label, i in label_to_int.items() if 'pre' in label]
antitag_cats = [ i for label, i in label_to_int.items() if 'at' in label]




toc = time.time()




if (len(sys.argv) > 1) and (sys.argv[1] in ['2016', '2016APV', '2017', '2018']):
    
    year = sys.argv[1]

else:

    year = '2016'


IOVs = [year]
    
# IOVs = [
#     '2016APV',
#     '2016',
#     '2017', 
#     '2018'
# ]

for IOV in IOVs:
    
    coffeafiles = functions.getCoffeaFilenames()

    datasets = [
        'QCD',
        'TTbar', 
        'JetHT', 
        'RSGluon',  
        'ZPrime10', 
        'ZPrime30', 
        'ZPrimeDM',
        'ZPrime1'
    ]
    

    hasBkgEst = False
    blind = False #False if '2016' in IOV else True


    bkgest_str = '_bkgest' if hasBkgEst else ''
    blind_str = '_blind' if blind else ''
    filetype = 'weighted' if hasBkgEst else 'unweighted'

    lumifactor = 0.1 if blind else 1.0
    
    
#     scaledir = 'outputs/tight_fixed/scale/'
    
#     outputdir = 'outputs/loosenottight/'
    # outputdir = 'outputs/tight_jetid1_ak8pftrigger_plots/'
    outputdir = 'outputs/plot_3d/'

    scaledir = outputdir+'/scale/'
    
    for ds in datasets:


        try:

            coffeafiles[ds][filetype][IOV].keys()
            sections = coffeafiles[ds][filetype][IOV].keys()

            files = []


            for s in sections:

                filename = coffeafiles[ds][filetype][IOV][s]
                
                
                
                filename = filename.replace('outputs/', outputdir)




                if 'JetHT' in ds and blind:
                    filename = filename.replace('.coffea', '_blind.coffea')
                                
                original_file = util.load(filename)

                file = copy.deepcopy(original_file)

                if 'JetHT' in ds:

                    files.append(file)

                else:

                    factor = 1.0 #functions.toptag_sf**2 if 'TTbar' in ds else 1.0
                    sf = functions.lumi[IOV] * lumifactor * functions.xs[ds][s] * factor / file['cutflow']['sumw']
                    sf_events = functions.lumi[IOV] * lumifactor * functions.xs[ds][s] / file['cutflow']['all events']
                    
                    for key in file.keys():

                        if 'hist' in str(type(file[key])) and not 'nocut' in key:
                            file[key] = file[key] * sf

                        elif 'accumulator' in str(type(file[key])):
                            for cut in file[key].keys():
                                file[key][cut] = file[key][cut] * sf_events



                    files.append(file)

                    if 'RSGluon' in ds: 

                        util.save(file, f'../{scaledir}{ds}{s}_{IOV}{blind_str}.coffea')
                        print(f'saving ../{scaledir}{ds}{s}_{IOV}{blind_str}.coffea')

                    elif 'ZPrime' in ds:

                        util.save(file, f'../{scaledir}ZPrime{s}_{ds.replace("ZPrime","")}_{IOV}{blind_str}.coffea')
                        print(f'saving ../{scaledir}ZPrime{s}_{ds.replace("ZPrime","")}_{IOV}{blind_str}.coffea')


            if 'RSGluon' not in ds and 'ZPrime' not in ds:

                file = files[0]

                for f in files[1:]:
                    for key in file.keys():

                        if 'hist' in str(type(file[key])) and not 'nocut' in key:
                            file[key] = file[key] + f[key]
                            
                        elif 'accumulator' in str(type(file[key])):
                            for cut in f[key].keys():
                                f[key][cut] = f[key][cut] + file[key][cut]

                savefilename = f'../{scaledir}{ds}_{IOV}{bkgest_str}{blind_str}.coffea'
                util.save(file, savefilename)
                print(f'saving {savefilename}')

        except:

            filename = coffeafiles[ds][filetype][IOV]
            filename = filename.replace('outputs/', outputdir)

            original_file = util.load(filename)
            file = copy.deepcopy(original_file)


            sf = functions.lumi[IOV] * lumifactor * functions.xs[ds] / file['cutflow']['sumw']
            sf_events = functions.lumi[IOV] * lumifactor * functions.xs[ds] / file['cutflow']['all events']


            for key in file.keys():
                

                if 'hist' in str(type(file[key])):
                    file[key] = file[key] * sf

                elif 'accumulator' in str(type(file[key])):

                    for cut in file[key].keys():

                        file[key][cut] = file[key][cut] * sf_events

            savefilename = f'../{scaledir}{ds}_{IOV}{bkgest_str}{blind_str}.coffea'
            util.save(file, savefilename)
            print(f'saving {savefilename}')


if year=='2016APV':
    
    print('combing 2016 and 2016noAPV')
    
    datasets = []

    IOVs = ['2016', '2016APV']
    datasets = [
        'TTbar', 
        'QCD', 
        'JetHT'
               ]


    hasBkgEst = False
    bkgest_str = '_bkgest' if hasBkgEst else ''


    # add all signal samples
    datasets += ['RSGluon'+str(int(b*100)) for b in [10,15,20,25,30,35,40,45,50,55,60]]
    datasets += ['ZPrime'+str(int(b*100))+'_1' for b in [10,12,14,16,18,20,25,30,35,40,45]]
    datasets += ['ZPrime'+str(int(b*100))+'_10' for b in [10,12,14,16,18,20,25,30,35,40,45,50,60,70]]
    datasets += ['ZPrime'+str(int(b*100))+'_30' for b in [10,12,14,16,18,20,25,30,35,40,45,50,60,70]]
    datasets += ['ZPrime'+str(int(b*100))+'_DM' for b in [10,15,20,25,30,35,40,45,50]]

    for ds in datasets:

        files = []
        for IOV in IOVs:
            

            file = util.load(f'../{scaledir}{ds}_{IOV}{bkgest_str}.coffea')
            files.append(file)
            
        


        file = files[0]
        for f in files[1:]:
            for key in file.keys():

                if 'hist' in str(type(file[key])) and not 'nocut' in key:

                    file[key] = file[key] + f[key]

                elif 'cutflow' in key:
                    for cut in f[key].keys():
                        f[key][cut] = f[key][cut] + file[key][cut]  

        savefilename = f'../{scaledir}{ds}_2016all{bkgest_str}.coffea'
        util.save(file, savefilename)
        print(f'saving {savefilename}')

