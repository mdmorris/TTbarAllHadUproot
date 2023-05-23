# ttbaranalysis.py

import time
from coffea import nanoevents, util
import hist
import coffea.processor as processor
import awkward as ak
import numpy as np
import glob as glob
import itertools
from coffea.nanoevents.methods import vector
from coffea.nanoevents import NanoAODSchema, BaseSchema
import os

import warnings
warnings.filterwarnings("ignore")

import inspect
import json
import glob
from dask.distributed import Client, performance_report
from lpcjobqueue import LPCCondorCluster
from dask.distributed.diagnostics.plugin import UploadDirectory
import argparse

ak.behavior.update(vector.behavior)


savedir = 'outputs/'

from ttbarprocessor import TTbarResProcessor

if __name__ == "__main__":
    
    tic = time.time()
        
    parser = argparse.ArgumentParser(
                    prog='ttbaranalysis.py',
                    description='Run ttbarprocessor',
                    epilog='help')
    
    
    # run options
    parser.add_argument('--dask', action='store_true')
    parser.add_argument('--test', action='store_true')
    
    # datasets to run
    parser.add_argument('-d', '--dataset', choices=['JetHT', 'QCD', 'TTbar'], default='QCD')
    parser.add_argument('--iov', choices=['2016APV', '2016', '2017', '2018'], default='2016APV')
    parser.add_argument('--era', choices=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'], action='append', default='')

    # platforms
    parser.add_argument('-C', '--casa', action='store_true', help='Use Coffea-Casa redirector: root://xcache/')
    parser.add_argument('-L', '--lpc', action='store_true', help='Use CMSLPC redirector: root://cmsxrootd.fnal.gov/')
    parser.add_argument('-W', '--winterfell', action='store_true', help='Get available files from UB Winterfell /mnt/data/cms')
    
    args = parser.parse_args()
    print('args', args)

    
    

    if args.casa: redirector = 'root://xcache/'
    elif args.winterfell: redirector = '/mnt/data/cms/'
    else: redirector = 'root://cmsxrootd.fnal.gov/' # default LPC
    
    sample = args.dataset
    IOV = args.iov  
    jsonfiles = {
        "JetHT": 'data/JetHT.json',
        "QCD": 'data/QCD.json',
        "TTbar": 'data/TTbar.json',
    }
        
    
    # get root files
    
    inputfile = jsonfiles[sample]
    files = []
    with open(inputfile) as json_file:
        
        data = json.load(json_file)
        
        rfiles = []
        if 'JetHT' in sample:
            
            # loop through all eras
            if len(args.era) < 1:
                for era in data[IOV].keys():
                    rfiles.append(data[IOV][era])
                    
            # loop through specified eras
            else:
                for era in args.era:
                    if era in data[IOV].keys():
                        rfiles.append(data[IOV][era])
                    else:
                        print(f'no era {era} in {IOV}')
        else:
            rfiles = data[IOV]
            
        if args.test: rfiles = [rfiles[0]]

        files = [redirector + f for f in rfiles]
    
        
    fileset = {sample: files}
    
    
    
    
    # run uproot jobs
    
    testString = ''
    eraString = ''
    
    if args.test: testString = '_test'
    for era in args.era: eraString = eraString + era
    
    savefilename = f'{savedir}{sample}_{IOV}{eraString}{testString}.coffea'  
    print(f'Running {IOV}{eraString} {sample} {testString[1:]}')
    
    if args.dask:
        
        exe_args = {
            "client": client,
            "skipbadfiles": True,
            "savemetrics": True,
            "schema": NanoAODSchema,
        }

        cluster = LPCCondorCluster(memory='6GB')
        cluster.adapt(minimum=1, maximum=100)
        client = Client(cluster)
        client.upload_file('ttbarprocessor.py')
        client.upload_file('corrections/corrections.py')
        client.upload_file('corrections/btagCorrections.py')
        client.upload_file('corrections/functions.py')
        client.upload_file('corrections/subjet_btagging.json.gz')

        hists, metrics = processor.run_uproot_job(
            fileset,
            treename="Events",
            processor_instance=TTbarResProcessor(iov=IOV),
            executor=processor.dask_executor,
            executor_args=exe_args,
            chunksize=100000,
        )

    else:
        
        exe_args = {
            "skipbadfiles": True,
            "savemetrics": True,
            "schema": NanoAODSchema,
            "workers":4
        }

        print("Waiting for at least one worker...")
        
        hists, metrics = processor.run_uproot_job(
            fileset,
            treename="Events",
            processor_instance=TTbarResProcessor(iov=IOV),
            executor=processor.iterative_executor,
            executor_args=exe_args,
            chunksize=100000,
        )


    util.save(hists, savefilename)
    print('saving', savefilename)

    elapsed = time.time() - tic
    print(f"\nFinished in {elapsed:.1f}s")
    print(f"Events/s: {metrics['entries'] / elapsed:.0f}")



