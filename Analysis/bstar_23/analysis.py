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

    #'Mmmg':hist.Hist("Events",hist.Cat('process', 'process'),hist.Bin("Mmmg", "Mmmg", 20,50,250)),
    #'Mmm':hist.Hist("Events",hist.Cat('process', 'process'),hist.Bin("Mmm", "Mmm", 75,50,200)),
    #'Mmm_zoom':hist.Hist("Events",hist.Cat('process', 'process'),hist.Bin("Mmm", "Mmm", 20,50,120)),

    #'Mtaumg':hist.Hist("Events",hist.Cat('process', 'process'),hist.Bin("Mtaumg", "Mtaumg", 20,50,250)),
    #'Mtaum':hist.Hist("Events",hist.Cat('process', 'process'),hist.Bin("Mtaum", "Mtaum", 75,50,200)),
    #'Mtaum_zoom':hist.Hist("Events",hist.Cat('process', 'process'),hist.Bin("Mtaum", "Mtaum", 20,50,120)),
    
    'MET':hist.Hist('Events',hist.Cat('process', 'process'),hist.Bin('MET', '$MET$', 100,0,200)),
    #'PhotonpT':hist.Hist("Events",hist.Cat('process', 'process'),hist.Bin("PhotonpT", "PhotonpT", 150,0,150)),
    #'LeadingMuonpT':hist.Hist("Events",hist.Cat('process', 'process'),hist.Bin("LeadingMuonpT", "LeadingMuonpT", 150,0,150)),
    
    #'LeadingMuonpT_mtau':hist.Hist("Events",hist.Cat('process', 'process'),hist.Bin("LeadingMuonpT", "LeadingMuonpT", 150,0,150)),

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
         else:   
            ev["weight_norm"]=1
    else:
        ev["weight_norm"]=1
    
    hists.fill('MET',ev.weight_norm,(ev.event>0), ev.MET, MET="pt")
    #-------------------------------------------------------------------------------------------------------
    ############################# KEEP THIS BLOCK ##################
    ET.autolog(f"{len(ev)} Events at the end of your analysis", logger, 'i')
    return hists.get()
    ############################# KEEP THIS BLOCK ##################
