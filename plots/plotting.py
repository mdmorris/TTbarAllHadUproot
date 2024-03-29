# plotting.py

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
from scipy.optimize import curve_fit
import uproot

sys.path.append('../python/')
import functions

# suppress warnings
import warnings
warnings.filterwarnings("ignore")


# get IOV from command line

if (len(sys.argv) > 1) and (sys.argv[1] in ['2016', '2016APV', '2016all', '2017', '2018', 'all']):
    
    IOV = sys.argv[1]

else:
    
    IOV = '2016'

    
print(IOV, type(IOV))

# scale factors and luminosity

lumi = functions.lumi
print(lumi)
rsgluon_xs = functions.rsgluon_xs



# initialize


print('making images')
functions.makeSaveDirectories(coffea_dir='../outputs/tight_jetid1_ak8pftrigger/')


blind = False
lumifactor = 0.1 if blind else 1.0

label_map = functions.getLabelMap(coffea_dir='../outputs/tight_jetid1_ak8pftrigger/')
label_to_int = {label: i for i, label in label_map.items()}

signal_cats = [ i for label, i in label_to_int.items() if '2t' in label]
signal_cen_cats = [ i for label, i in label_to_int.items() if '2t' in label and 'cen' in label]
signal_fwd_cats = [ i for label, i in label_to_int.items() if '2t' in label and 'fwd' in label]

pretag_cats = [ i for label, i in label_to_int.items() if 'pre' in label]
antitag_cats = [ i for label, i in label_to_int.items() if 'at' in label]
antitag_cen_cats = [ i for label, i in label_to_int.items() if 'at' in label and 'cen' in label]
antitag_fwd_cats = [ i for label, i in label_to_int.items() if 'at' in label and 'fwd' in label]


lines_dict = {'solid': 'solid',
 'dotted': (0, (1, 1)),
 'dashed': (0, (5, 5)),
 'dashdot': 'dashdot',
 'loosely dotted': (0, (1, 10)),
 'densely dotted': (0, (1, 1)),
 'long dash with offset': (5, (10, 3)),
 'loosely dashed': (0, (5, 10)),
 'densely dashed': (0, (5, 1)),
 'loosely dashdotted': (0, (3, 10, 1, 10)),
 'dashdotted': (0, (3, 5, 1, 5)),
 'densely dashdotted': (0, (3, 1, 1, 1)),
 'dashdotdotted': (0, (3, 5, 1, 5, 1, 5)),
 'loosely dashdotdotted': (0, (3, 10, 1, 10, 1, 10)),
 'densely dashdotdotted': (0, (3, 1, 1, 1, 1, 1))}

lines = list(lines_dict.values())


# categories and systematics labels

# cats = ['0bcen', '0bfwd', '1bcen', '1bfwd', '2bcen', '2bfwd']
# cat_labels = ['cen0b', 'fwd0b', 'cen1b', 'fwd1b', 'cen2b', 'fwd2b']


cats = ['cen', 'fwd']
cat_labels=['cen', 'fwd']


def getUncertainy(hbkg, hUnc):
    
    axes_systematics = hUnc.axes['systematic']
    
    nomvals = hbkg.values()

    for syst in axes_systematics:

        if 'Up' in syst:

            upvals = (1 + np.abs(nomvals - hUnc[{'systematic':syst}].values())/nomvals)

        elif 'Down' in syst:

            downvals = (1 - np.abs(nomvals - hUnc[{'systematic':syst}].values())/nomvals)
            
    return hbkg*upvals, hbkg*downvals




        
# systematics plots #

def plotSystematics(IOV, dataset='TTbar', dirstring=''):

    print('\nPlotting systematics\n')
    
    print(IOV)


    coffea_dir = f'../outputs/{dirstring}/scale/'
    
    

    for cat, catname in zip(['']+cats, ['inclusive']+cat_labels):

        
        
        if catname == 'inclusive':
            
            signal_cat = signal_cats

        elif catname == 'cen':
            
            signal_cat = signal_cen_cats

        elif catname == 'fwd':
            
            signal_cat = signal_fwd_cats
            
        else:
            signal_cat = label_to_int['2t'+cat]
            
        
        
        systsUpDn = list(functions.getHist2('ttbarmass', dataset, IOV, sum_axes=['anacat'], coffea_dir=coffea_dir).axes['systematic'])
        
        systematics = [syst.replace('Up', '') for syst in systsUpDn if 'Down' not in syst]        
        

        for syst in systematics[1:]:

            fig, (ax1, ax2) = plt.subplots(nrows=2, height_ratios=[3, 1])

            text = r'MC TTbar'+'\n'+syst+' systematic variations'

            dytext = ''
            if 'cen' in cat:
                dytext = r'$\bf{central}$'
            elif 'fwd' in cat:
                dytext = r'$\bf{foward}$'

            btext = ''
            if '0b' in cat:
                btext = '0 b-tags'
            elif '1b' in cat:
                btext = '1 b-tag'
            elif '2b' in cat:
                btext = '2 b-tags'

                
            if catname == 'inclusive':
                text = f'MC {dataset}\n{syst} systematic variations\n' + r'b-tag, $\Delta y$ inclusive'
            elif catname == 'cen':
                text = f'MC {dataset}\n{syst} systematic variations\n{dytext}'
            elif catname == 'fwd':
                text = f'MC {dataset}\n{syst} systematic variations\n{dytext}'
            else:
                text = f'MC {dataset}\n{syst} systematic variations\n{btext}, {dytext}'


            hep.cms.label('Preliminary', data=True, lumi='{0:0.1f}'.format(lumi[IOV]*lumifactor/1000.), year=IOV.replace('all',''), loc=1, fontsize=20, ax=ax1)
            hep.cms.text(text, loc=2, fontsize=20, ax=ax1)


            
            if 'hem' in syst:
                
                httbar = functions.getHist2('ttbarmass', dataset, IOV,
                     sum_axes=['anacat'],
                     integrate_axes={'systematic':'nominal', 'anacat':signal_cats},
                     coffea_dir=coffea_dir
                    )
                httbarDn = functions.getHist2('ttbarmass', dataset, IOV,
                     sum_axes=['anacat'],
                     integrate_axes={'systematic':syst, 'anacat':signal_cats},
                     coffea_dir=coffea_dir
                    )


                hep.histplot(httbar, histtype='step', color='k', ax=ax1, label='Nominal')
                hep.histplot(httbarDn, histtype='step', color='red', ax=ax1, label=syst)


                ratioDn = httbarDn / httbar.values()

                hep.histplot(ratioDn, histtype='step', color='red', ax=ax2, label=syst)
                
                
            else:
                httbar = functions.getHist2('ttbarmass', dataset, IOV,
                     sum_axes=['anacat'],
                     integrate_axes={'systematic':'nominal', 'anacat':signal_cats},
                     coffea_dir=coffea_dir
                    )
                httbarUp = functions.getHist2('ttbarmass', dataset, IOV,
                     sum_axes=['anacat'],
                     integrate_axes={'systematic':syst+'Up', 'anacat':signal_cats},
                     coffea_dir=coffea_dir
                    )
                httbarDn = functions.getHist2('ttbarmass', dataset, IOV,
                     sum_axes=['anacat'],
                     integrate_axes={'systematic':syst+'Down', 'anacat':signal_cats},
                     coffea_dir=coffea_dir
                    )

                hep.histplot(httbar, histtype='step', color='k', ax=ax1, label='Nominal')
                hep.histplot(httbarUp, histtype='step', color='green', ax=ax1, label='Up')
                hep.histplot(httbarDn, histtype='step', color='red', ax=ax1, label='Down')


                ratioUp = httbarUp / httbar.values()
                ratioDn = httbarDn / httbar.values()

                hep.histplot(ratioUp, histtype='step', color='green', ax=ax2)
                hep.histplot(ratioDn, histtype='step', color='red', ax=ax2)
            ax2.axhline(1, color='black', ls='--')
            
            ymax = np.max(httbar.values()) * 1.5

            ax2.set_ylabel('Syst/Nom')
            ax2.set_xlabel(ax1.get_xlabel())
            ax1.set_xlabel('')
            ax1.set_ylim(1e-1,ymax)
            ax2.set_ylim(0.5,1.5)


            ax1.legend()

            imagefile = f'images/{dirstring}/png/systematics/{IOV}/{dataset}_{catname}_{syst}.png'

            plt.savefig(imagefile)
            plt.savefig(imagefile.replace('png','pdf'))
            print('saving ', imagefile)
            print('saving ', imagefile.replace('png','pdf'))

            ax1.plot()



                
def plotClosureTest():

    
    print('\nPlotting Closure Test\n')
    

    fig, (ax1, ax2) = plt.subplots(nrows=2, height_ratios=[3, 1])



    hbkg = functions.getHist2('ttbarmass', 'JetHT', IOV,
             sum_axes=['anacat'],
             integrate_axes={'systematic':'nominal', 'anacat':antitag_cats},
             tag = '_bkgest'
            )

    httbar = functions.getHist2('ttbarmass', 'TTbar', IOV,
             sum_axes=['anacat'],
             integrate_axes={'systematic':'nominal', 'anacat':signal_cats},        
            )

    hsig = functions.getHist2('ttbarmass', 'JetHT', IOV,
             sum_axes=['anacat'],
             integrate_axes={'systematic':'nominal', 'anacat':signal_cats}        
            )


    # 2D histogram with all uncertainties for getting background uncertainty
    hUnc = functions.getHist2('ttbarmass', 'JetHT', IOV,
             sum_axes=['anacat'],
             integrate_axes={'anacat':antitag_cats},
             tag = '_bkgest')


    hUp, hDn = getUncertainy(hbkg, hUnc)


    text = 'Preliminary'+'\n'+r'$\Delta y$ inclusive'+', '+r'b-tag inclusive'

    hep.cms.label('', data=True, lumi='{0:0.1f}'.format(functions.lumi[IOV]*lumifactor/1000), year=IOV.replace('all',''), loc=2, fontsize=20, ax=ax1)
    hep.cms.text(text, loc=2, fontsize=20, ax=ax1)

    hep.histplot(hsig, histtype='errorbar', color='black', label='Data', ax=ax1)
    hep.histplot(hbkg+httbar, histtype='fill', color='xkcd:pale gold', label='NTMJ Bkg Est', ax=ax1)
    hep.histplot(httbar, histtype='fill', color='xkcd:deep red', label='SM TTbar', ax=ax1)

    height = hUp.values() + hDn.values()
    bottom = hbkg.values() - hDn.values()
    edges  = hbkg.axes['ttbarmass'].edges


    ax1.bar(x = edges[:-1],
               height=height,
               bottom=bottom,
               width = np.diff(edges), align='edge', hatch='//////', edgecolor='gray',
               linewidth=0, facecolor='none', alpha=0.7,
               zorder=10, label='Unc.')


    ratio_plot =  hsig / hbkg.values()
    ratioUp = hUp.values() / hbkg.values()
    ratioDn = hUp.values() / hbkg.values()


    ax2.bar(x = edges[:-1],
               height=(np.ones_like(ratio_plot.values()) + ratioUp + ratioDn),
               bottom=(np.ones_like(ratio_plot.values()) - ratioDn),
               width = np.diff(edges), align='edge', edgecolor='gray',
               linewidth=0, facecolor='gray', alpha=0.3,
               zorder=10, label='Unc.')


    ratio_plot =  hsig / hbkg.values()
    hep.histplot(ratio_plot, ax=ax2, histtype='errorbar', color='black')
    ax2.set_ylim(-10,10)
    ax2.axhline(1, color='black', ls='--')
    ax2.set_ylabel('Data/Bkg')


    ax1.legend()
    ax1.set_ylabel(f'Events / Bin GeV'.replace('j',''))
    ax1.set_yscale('log')
    ax1.set_ylim(1e-1, 1e6)
    ax1.set_xlim(900,6000)
    ax2.set_xlim(900,6000)
    ax1.set_xlabel('')


    savefigname = f'images/{dirstring}/png/closureTest/{IOV}/closure_inclusive.png'
    plt.savefig(savefigname)
    plt.savefig(savefigname.replace('png', 'pdf'))

    print('saving '+savefigname)
    print('saving '+savefigname.replace('png', 'pdf'))






    # plot regions

    for cat in cats:


        fig, (ax1, ax2) = plt.subplots(nrows=2, height_ratios=[3, 1])


        signal_cat = label_to_int['2t'+cat]
        antitag_cat = label_to_int['at'+cat]


        hbkg = functions.getHist2('ttbarmass', 'JetHT', IOV,
             sum_axes=[],
             integrate_axes={'systematic':'nominal', 'anacat':antitag_cat},
             tag = '_bkgest'

                             )

        httbar = functions.getHist2('ttbarmass', 'TTbar', IOV,
                 sum_axes=[],
                 integrate_axes={'systematic':'nominal', 'anacat':signal_cat},        
                )

        hsig = functions.getHist2('ttbarmass', 'JetHT', IOV,
                 sum_axes=[],
                 integrate_axes={'systematic':'nominal', 'anacat':signal_cat}        
                )
        

        # 2D histogram with all uncertainties for getting background uncertainty
        hUnc = functions.getHist2('ttbarmass', 'JetHT', IOV,
                 sum_axes=[],
                 integrate_axes={'anacat':antitag_cat},
                 tag = '_bkgest')


        hUp, hDn = getUncertainy(hbkg, hUnc)


        dytext = ''
        if 'cen' in cat:
            dytext = r'$\Delta y$ < 1.0'
        elif 'fwd' in cat:
            dytext = r'$\Delta y$ > 1.0'

        btext = ''
        if '0b' in cat:
            btext = '0 b-tags'
        elif '1b' in cat:
            btext = '1 b-tag'
        elif '2b' in cat:
            btext = '2 b-tags'

    
        text = 'Preliminary'+'\n'+dytext+', '+ btext
        hep.cms.label('', data=True, lumi='{0:0.1f}'.format(functions.lumi[IOV]*lumifactor/1000), year=IOV.replace('all',''), loc=2, fontsize=20, ax=ax1)
        hep.cms.text(text, loc=2, fontsize=20, ax=ax1)

        hep.histplot(hsig, histtype='errorbar', color='black', label='Data', ax=ax1)
        hep.histplot(hbkg+httbar, histtype='fill', color='xkcd:pale gold', label='NTMJ Bkg Est', ax=ax1)
        hep.histplot(httbar, histtype='fill', color='xkcd:deep red', label='SM TTbar', ax=ax1)

        height = hUp.values() + hDn.values()
        bottom = hbkg.values() - hDn.values()
        edges  = hbkg.axes['ttbarmass'].edges


        ax1.bar(x = edges[:-1],
                   height=height,
                   bottom=bottom,
                   width = np.diff(edges), align='edge', hatch='//////', edgecolor='gray',
                   linewidth=0, facecolor='none', alpha=0.7,
                   zorder=10, label='Unc.')


        ratio_plot =  hsig / hbkg.values()
        ratioUp = hUp.values() / hbkg.values()
        ratioDn = hUp.values() / hbkg.values()


        ax2.bar(x = edges[:-1],
                   height=(np.ones_like(ratio_plot.values()) + ratioUp + ratioDn),
                   bottom=(np.ones_like(ratio_plot.values()) - ratioDn),
                   width = np.diff(edges), align='edge', edgecolor='gray',
                   linewidth=0, facecolor='gray', alpha=0.3,
                   zorder=10, label='Unc.')


        ratio_plot =  hsig / hbkg.values()
        hep.histplot(ratio_plot, ax=ax2, histtype='errorbar', color='black')
        ax2.set_ylim(-10,10)
        ax2.axhline(1, color='black', ls='--')
        ax2.set_ylabel('Data/Bkg')


        ax1.legend()
        ax1.set_ylabel(f'Events / Bin GeV'.replace('j',''))
        ax1.set_yscale('log')
        ax1.set_ylim(1e-1, 1e6)
        ax1.set_xlim(900,6000)
        ax2.set_xlim(900,6000)
        ax1.set_xlabel('')

        savefigname = f'images/{dirstring}/png/closureTest/{IOV}/closure_{cat}.png'
        plt.savefig(savefigname)
        plt.savefig(savefigname.replace('png', 'pdf'))

        print('saving '+savefigname)
        print('saving '+savefigname.replace('png', 'pdf'))


    
    

def plotClosureTestQCD():
    
    
    print('\nPlotting QCD Closure Test\n')


    fig, (ax1, ax2) = plt.subplots(nrows=2, height_ratios=[3, 1])

    hbkg = functions.getHist2('ttbarmass', 'QCD', IOV,
             sum_axes=['anacat'],
             integrate_axes={'systematic':'nominal', 'anacat':antitag_cats},
             tag = '_bkgest'

            )
    
    hsig = functions.getHist2('ttbarmass', 'QCD', IOV,
             sum_axes=['anacat'],
             integrate_axes={'systematic':'nominal', 'anacat':signal_cats},
             tag = ''
            )

    # 2D histogram with all uncertainties for getting background uncertainty
    hUnc = functions.getHist2('ttbarmass', 'QCD', IOV,
             sum_axes=['anacat'],
             integrate_axes={'anacat':antitag_cats},
             tag = '_bkgest')


    hUp, hDn = getUncertainy(hbkg, hUnc)


    text = 'Preliminary'+'\n'+r'$\Delta y$ inclusive'+', '+r'b-tag inclusive'

    hep.cms.label('', data=True, lumi='{0:0.1f}'.format(functions.lumi[IOV]*lumifactor/1000), year=IOV.replace('all','').replace('all',''), loc=2, fontsize=20, ax=ax1)
    hep.cms.text(text, loc=2, fontsize=20, ax=ax1)

    hep.histplot(hsig, histtype='errorbar', color='black', label='QCD SR', ax=ax1)
    hep.histplot(hbkg, histtype='fill', color='xkcd:pale gold', label='QCD Bkg Est', ax=ax1)
    herr = np.sqrt(hbkg.variances())

    height = hUp.values() + hDn.values()
    bottom = hbkg.values() - hDn.values()
    edges  = hbkg.axes['ttbarmass'].edges


    ax1.bar(x = edges[:-1],
               height=height,
               bottom=bottom,
               width = np.diff(edges), align='edge', hatch='//////', edgecolor='gray',
               linewidth=0, facecolor='none', alpha=0.7,
               zorder=10, label='Unc.')


    ratio_plot =  hsig / hbkg.values()
    ratioUp = hUp.values() / hbkg.values()
    ratioDn = hUp.values() / hbkg.values()


    ax2.bar(x = edges[:-1],
               height=(np.ones_like(ratio_plot.values()) + ratioUp + ratioDn),
               bottom=(np.ones_like(ratio_plot.values()) - ratioDn),
               width = np.diff(edges), align='edge', edgecolor='gray',
               linewidth=0, facecolor='gray', alpha=0.3,
               zorder=10, label='Unc.')

    hep.histplot(ratio_plot, ax=ax2, histtype='errorbar', color='black')
    ax2.set_ylim(-10,10)
    ax2.axhline(1, color='black', ls='--')
    ax2.set_ylabel('Data/Bkg')


    ax1.legend()
    ax1.set_ylabel(f'Events / Bin GeV'.replace('j',''))
    ax1.set_yscale('log')
    ax1.set_ylim(1e-1, 1e13)
    ax1.set_xlim(900,6000)
    ax2.set_xlim(900,6000)
    ax1.set_xlabel('')

    savefigname = f'images/{dirstring}/png/closureTest/{IOV}/closureQCD_inclusive.png'
    plt.savefig(savefigname)
    plt.savefig(savefigname.replace('png', 'pdf'))

    print('saving '+savefigname)
    print('saving '+savefigname.replace('png', 'pdf'))


    
    # plot regions
    
    for cat in cats:
        
        
        signal_cat = label_to_int['2t'+cat]
        antitag_cat = label_to_int['at'+cat]        

        fig, (ax1, ax2) = plt.subplots(nrows=2, height_ratios=[3, 1])

        hbkg = functions.getHist2('ttbarmass', 'QCD', IOV,
                 sum_axes=[],
                 integrate_axes={'systematic':'nominal', 'anacat':antitag_cat},
                 tag = '_bkgest'

                )

        hsig = functions.getHist2('ttbarmass', 'QCD', IOV,
                 sum_axes=[],
                 integrate_axes={'systematic':'nominal', 'anacat':signal_cat}        
                )


        
        # 2D histogram with all uncertainties for getting background uncertainty
        hUnc = functions.getHist2('ttbarmass', 'QCD', IOV,
                 sum_axes=[],
                 integrate_axes={'anacat':antitag_cat},
                 tag = '_bkgest')


        hUp, hDn = getUncertainy(hbkg, hUnc)



        dytext = ''
        if 'cen' in cat:
            dytext = r'$\Delta y$ < 1.0'
        elif 'fwd' in cat:
            dytext = r'$\Delta y$ > 1.0'

        btext = ''
        if '0b' in cat:
            btext = '0 b-tags'
        elif '1b' in cat:
            btext = '1 b-tag'
        elif '2b' in cat:
            btext = '2 b-tags'

    
        text = 'Preliminary'+'\n'+dytext+', '+ btext

        hep.cms.label('', data=True, lumi='{0:0.1f}'.format(functions.lumi[IOV]*lumifactor/1000), year=IOV.replace('all',''), loc=2, fontsize=20, ax=ax1)
        hep.cms.text(text, loc=2, fontsize=20, ax=ax1)

        hep.histplot(hsig, histtype='errorbar', color='black', label='QCD SR', ax=ax1)
        hep.histplot(hbkg, histtype='fill', color='xkcd:pale gold', label='QCD Bkg Est', ax=ax1)
        herr = np.sqrt(hbkg.variances())

        height = hUp.values() + hDn.values()
        bottom = hbkg.values() - hDn.values()
        edges  = hbkg.axes['ttbarmass'].edges


        ax1.bar(x = edges[:-1],
                   height=height,
                   bottom=bottom,
                   width = np.diff(edges), align='edge', hatch='//////', edgecolor='gray',
                   linewidth=0, facecolor='none', alpha=0.7,
                   zorder=10, label='Unc.')
        


        ratio_plot =  hsig / hbkg.values()
        ratioUp = hUp.values() / hbkg.values()
        ratioDn = hUp.values() / hbkg.values()


        ax2.bar(x = edges[:-1],
                   height=(np.ones_like(ratio_plot.values()) + ratioUp + ratioDn),
                   bottom=(np.ones_like(ratio_plot.values()) - ratioDn),
                   width = np.diff(edges), align='edge', edgecolor='gray',
                   linewidth=0, facecolor='gray', alpha=0.3,
                   zorder=10, label='Unc.')

        hep.histplot(ratio_plot, ax=ax2, histtype='errorbar', color='black')
        ax2.set_ylim(-10,10)
        ax2.axhline(1, color='black', ls='--')
        ax2.set_ylabel('Data/Bkg')


        ax1.legend()
        ax1.set_ylabel(f'Events / Bin GeV'.replace('j',''))
        ax1.set_yscale('log')
        ax1.set_ylim(1e-1, 1e13)
        ax1.set_xlim(900,6000)
        ax2.set_xlim(900,6000)
        ax1.set_xlabel('')


        savefigname = f'images/{dirstring}/png/closureTest/{IOV}/closureQCD_{cat}.png'
        plt.savefig(savefigname)
        plt.savefig(savefigname.replace('png', 'pdf'))

        print('saving '+savefigname)
        print('saving '+savefigname.replace('png', 'pdf'))




        
def plotKinematics(IOV, dataset='TTbar', dirstring=''):

    
    print('\nPlotting Kinematics \n')
    

    hists = [
        'jeteta',
        'jetphi',
        'ttbarmass',
        'jetmass',
        'jetmsd',
        'ht'
    ]

    for histname in hists:
        
        coffea_dir = f'../outputs/{dirstring}/'
        
        fttbar = util.load(coffea_dir+'/scale/TTbar_'+IOV+'.coffea')
        fqcd = util.load(coffea_dir+'/scale/QCD_'+IOV+'.coffea')
        fdata = util.load(coffea_dir+'/scale/JetHT_'+IOV+'.coffea')
        fsig = util.load(coffea_dir+'/scale/RSGluon2000_'+IOV+'.coffea')
        
        
        # h_ht = fttbar['jetpt']
        h_ttbar = fttbar[histname][{'systematic':'nominal'}][{'anacat':sum}]
        h_qcd = fqcd[histname][{'systematic':'nominal'}][{'anacat':sum}]
        
        h_data = fdata[histname][{'systematic':'nominal'}][{'anacat':sum}]
        
        h_bkg = (h_qcd + h_ttbar)
        
        h_sig = fsig[histname][{'systematic':'nominal'}][{'anacat':sum}]
        
        
        nomvals = fttbar[histname][{'systematic':'nominal'}][{'anacat':sum}].values()
        systUp2 = np.zeros_like(nomvals)
        systDn2 = np.zeros_like(nomvals)
    
    
            # rebinning and axes
    
    
    
        # h_bkg = h_bkg
        # h_ttbar = h_ttbar
        # h_data = h_data
        # h_sig = h_sig
        # h_qcd = h_qcd
    
    
    
    
        
        for syst in fttbar[histname].axes['systematic']:
            nomvals = h_ttbar.values()
        
            if 'Up' in syst:
                systvals = np.where(fttbar[histname][{'systematic':syst}][{'anacat':sum}].values()[0] > 0, fttbar[histname][{'systematic':syst}][{'anacat':sum}].values() - nomvals, 0) 
                systUp2 = systUp2 + systvals * systvals
            if 'Down' in syst:
                systvals = np.where(fttbar[histname][{'systematic':syst}][{'anacat':sum}].values()[0] > 0, fttbar[histname][{'systematic':syst}][{'anacat':sum}].values() - nomvals, 0) 
                systDn2 = systDn2 + systvals * systvals
        
        for syst in fqcd[histname].axes['systematic']:
            # nomvals = fqcd[histname][{'systematic':'nominal'}][{'anacat':sum}].values()
            nomvals = h_qcd.values()
        
            if 'Up' in syst:
                systvals = np.where(fqcd[histname][{'systematic':syst}][{'anacat':sum}].values()[0] > 0, fqcd[histname][{'systematic':syst}][{'anacat':sum}].values() - nomvals, 0) 
                systUp2 = systUp2 + systvals * systvals
            if 'Down' in syst:
                systvals = np.where(fqcd[histname][{'systematic':syst}][{'anacat':sum}].values()[0] > 0, fqcd[histname][{'systematic':syst}][{'anacat':sum}].values() - nomvals, 0) 
                systDn2 = systDn2 + systvals * systvals
        
        
        systUp = np.sqrt(systUp2)
        systDn = np.sqrt(systDn2)
        
        
        
        fig, (ax1, ax2) = plt.subplots(nrows=2, height_ratios=[3, 1], sharex=True)
        # fig, ax1 = plt.subplots(nrows=2, )
        
        
        text=''
    
    
        
        
        hep.cms.label('Work In Progress', data=True, lumi='{0:0.1f}'.format(lumi[IOV]/1000.), year=IOV.replace('all',''), loc=1, fontsize=20, ax=ax1)
        hep.cms.text(text, loc=2, fontsize=20, ax=ax1)
        hep.histplot(h_bkg, density=False, histtype='fill', color='xkcd:pale gold', label='SM QCD', ax=ax1)
        hep.histplot(h_ttbar, density=False, histtype='fill', color='xkcd:deep red', label='SM TTbar', ax=ax1)
        hep.histplot(h_data, density=False, histtype='errorbar', color='black', label='Data', ax=ax1)
        hep.histplot(h_sig, density=False, histtype='step', color='black', label='$g_{KK}$ (2 TeV)', ax=ax1)
    
    
        if 'jetmsd' in histname:
            edges = h_ttbar.axes['jetmass'].edges
        else:
            edges = h_ttbar.axes[histname].edges
    
        height = systUp + systDn
        bottom = h_bkg.values() - systDn
        
        
        
        ax1.bar(x = edges[:-1],
                   height=height,
                   bottom=bottom,
                   width = np.diff(edges), align='edge', hatch='///', edgecolor='gray',
                   linewidth=0, facecolor='none', alpha=0.8,
                   zorder=10, label='Syst. Unc.')
        
        
        ax2.bar(x = edges[:-1],
                   height=height/h_bkg.values(),
                   bottom=bottom/h_bkg.values(),
                   width = np.diff(edges), align='edge',
                   linewidth=0, facecolor='lightgrey', alpha=0.8,
                   zorder=0, label='Syst. Unc.')
        ax2.axhline(1, color='black', ls='--')
        
        hep.histplot(h_data/h_bkg.values(), density=False, histtype='errorbar', color='black', ax=ax2)
        
        
        
        
        # print((h_data/h_bkg.values()).values())
        # print(np.average((h_data/h_bkg.values()).values()))
        
        # ax2.bar(x = edges[:-1],
        #            height=1.5,
        #            bottom=0.5,
        #            width = np.diff(edges), align='edge', edgecolor='gray',
        #            linewidth=0, facecolor='gray', alpha=0.8,
        #            zorder=10, label='Syst. Unc.')
        
    
        
        # elif 'jetphi' in histname:
        #     plt.savefig('jetphi_'+IOV+'.pdf')
    
    
        if 'jetmass' in histname:
            ax1.set_xlim(25,500)
            ax2.set_xlim(25,500)
            ax1.set_xlabel('Jet mass [GeV]')
    
        elif 'jetmsd' in histname:
            ax1.set_xlim(25,500)
            ax2.set_xlim(25,500)   
    
        elif 'ttbarmass' in histname:
            ax1.set_xlim(800,6800)
            ax2.set_xlim(800,6800)
    
    
        
        elif 'jetpt' in histname:
            ax1.set_xlim(400,2000)
            ax2.set_xlim(400,2000)
        
        
        elif 'ht' in histname:
            ax1.set_xlim(500,4500)
            ax2.set_xlim(500,4500)
        
        elif 'jeteta' in histname:
            ax1.set_xlim(-2.4,2.4)
            ax2.set_xlim(-2.4,2.4)
    
        elif 'jetphi' in histname:
            ax1.set_xlim(-3,3)
            ax2.set_xlim(-3,3)
        
        
        ax1.set_ylim(1e-1, 1e7)
        ax2.set_ylim(1-0.5,1+0.5)
        
        
        ax1.legend(loc=1, fontsize=10)
        # plt.axvline(-1, color='darkred', ls='--')
        # plt.axvline(1, color='darkred', ls='--')
        # ax1.set_xlim(400,None)
        ax1.set_yscale('log')
        ax1.set_xlabel('')
        
        ax1.set_ylabel('Events/Bin')
        ax2.set_ylabel('Data/MC')
        
        
        savefigname = f'images/{dirstring}/png/kinematics/{IOV}/{histname}.png'
        print('saving', savefigname)
        plt.savefig(savefigname)
        plt.savefig(savefigname.replace('png', 'pdf'))


plotKinematics(IOV, dataset='TTbar', dirstring='tight_jetid1_ak8pftrigger')
# plotKinematics(IOV,  dataset='RSGluon2000', dirstring='tight_jetid1_ak8pftrigger')
# plotKinematics(IOV, dataset='ZPrime2000_1', dirstring='tight_jetid1_ak8pftrigger')
# plotKinematics(IOV, dataset='ZPrime2000_10', dirstring='tight_jetid1_ak8pftrigger')
# plotKinematics(IOV, dataset='ZPrime2000_30', dirstring='tight_jetid1_ak8pftrigger')
# plotKinematics(IOV, dataset='ZPrime2000_DM', dirstring='tight_jetid1_ak8pftrigger')   
# plotKinematics(IOV, dataset='RSGluon4000', dirstring='tight_jetid1_ak8pftrigger')
# plotKinematics(IOV, dataset='ZPrime4000_1', dirstring='tight_jetid1_ak8pftrigger')
# plotKinematics(IOV, dataset='ZPrime4000_10', dirstring='tight_jetid1_ak8pftrigger')
# plotKinematics(IOV, dataset='ZPrime4000_30', dirstring='tight_jetid1_ak8pftrigger')
# plotKinematics(IOV, dataset='ZPrime4000_DM', dirstring='tight_jetid1_ak8pftrigger')   
      

# plotSystematics(IOV, dataset='QCD', dirstring='tight_jetid1_ak8pftrigger')
# plotSystematics(IOV, dataset='TTbar', dirstring='tight_jetid1_ak8pftrigger')
# plotSystematics(IOV, dataset='RSGluon2000', dirstring='tight_jetid1_ak8pftrigger')
# plotSystematics(IOV, dataset='ZPrime2000_1', dirstring='tight_jetid1_ak8pftrigger')
# plotSystematics(IOV, dataset='ZPrime2000_10', dirstring='tight_jetid1_ak8pftrigger')
# plotSystematics(IOV, dataset='ZPrime2000_30', dirstring='tight_jetid1_ak8pftrigger')
# plotSystematics(IOV, dataset='ZPrime2000_DM', dirstring='tight_jetid1_ak8pftrigger')

# plotSystematics(IOV, dataset='QCD')

# plotClosureTest()
# plotClosureTestQCD()
# plotMtt()

   
