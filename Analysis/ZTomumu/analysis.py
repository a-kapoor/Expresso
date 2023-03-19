import awkward as ak, coffea, copy
from coffea import hist, analysis_tools
from modules.selection import *
import modules.ExpressoTools as ET
import modules.objects as obj
from coffea.analysis_tools import PackedSelection
import numpy.ma as ma
import numpy as np
from modules.hcoll import binning
    ############################# KEEP THIS BLOCK ##################
histograms = {
    ############################# KEEP THIS BLOCK ##################

    'Mmmg':hist.Hist("Events",hist.Cat('process', 'process'),hist.Bin("Mmmg", "Mmmg", 20,50,250)),
    'Mmm':hist.Hist("Events",hist.Cat('process', 'process'),hist.Bin("Mmm", "Mmm", 20,50,200)),
    'Mmm_zoom':hist.Hist("Events",hist.Cat('process', 'process'),hist.Bin("Mmm", "Mmm", 20,50,120)),
    'MET':hist.Hist('Events',hist.Cat('process', 'process'),hist.Bin('MET', '$MET$', 100,0,200)),
    'PhotonpT':hist.Hist("Events",hist.Cat('process', 'process'),hist.Bin("PhotonpT", "PhotonpT", 150,0,150)),
    'LeadingMuonpT':hist.Hist("Events",hist.Cat('process', 'process'),hist.Bin("LeadingMuonpT", "LeadingMuonpT", 150,0,150)),

    ############################# KEEP THIS BLOCK ##################
}
    ############################# KEEP THIS BLOCK ##################


############################# KEEP THIS BLOCK ##################
def myanalysis(pars, logger, h, ev, doweight=True):

    dataset,isData,histAxisName,year=pars['dataset'],pars['isData'],pars['histAxisName'],pars['year']
    xsec,sow,pass_options,analysis_point=pars['xsec'],pars['sow'],pars['passoptions'],pars['analysispoint']
    nEvents=pars['nEvents']

    from modules.hcoll import hcoll,binning
    hists = hcoll(h, isData, xsec, sow, doweight, process=histAxisName)
    ET.autolog(f"{len(ev)} Events at the start of your analysis", logger, 'i')
############################# KEEP THIS BLOCK ##################

    # Start your analysis``
    #-------------------------------------------------------------------------------------------------------

    if not isData:
         if pass_options=='Xsecweight':
            genw = ev["genWeight"]
            ev["weight_norm"] = (xsec / sow)* genw
            #ev["weight_norm"] =1
         else:   
            ev["weight_norm"]=1
    else:
        ev["weight_norm"]=1
            
    ev['M_uu']=(ev.recoMu[:,0]+ev.recoMu[:,1]).mass
    #ev['M_uug']=(ev.recoMu[:,0]+ev.recoMu[:,1]+ev.recoPho[:,0]).mass

    #ev['PhotonpT']=ev.recoPho[:,0].pt
    ev['LeadingMuonpT']=ev.recoMu[:,0].pt

    ev=ev[ev.M_uu>55]
    
    #hists.fill('Mmmg',ev.weight_norm,(ev.event>0), ev, Mmmg="M_uug")
    hists.fill('Mmm',ev.weight_norm,(ev.event>0), ev, Mmm="M_uu")
    #hists.fill('Mmm_zoom',ev.weight_norm,(ev.event>0), ev, Mmm="M_uu")
    #hists.fill('PhotonpT',ev.weight_norm,(ev.event>0), ev, PhotonpT="PhotonpT")
    hists.fill('LeadingMuonpT',ev.weight_norm,(ev.event>0), ev, LeadingMuonpT="LeadingMuonpT")
    

    ############################# KEEP THIS BLOCK ##################
    ET.autolog(f"{len(ev)} Events at the end of your analysis", logger, 'i')
    return hists.get()
    ############################# KEEP THIS BLOCK ##################
