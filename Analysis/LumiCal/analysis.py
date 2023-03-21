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
    'LumiCal':hist.Hist("root files",hist.Cat('process', 'process'),hist.Bin("Lumi", "Lumi (pb)", 100,0,100)),

    ############################# KEEP THIS BLOCK ##################
}
    ############################# KEEP THIS BLOCK ##################


############################# KEEP THIS BLOCK ##################
def myanalysis(pars, logger, h, ev, doweight=True):

    dataset,isData,histAxisName,year=pars['dataset'],pars['isData'],pars['histAxisName'],pars['year']
    xsec,sow,pass_options,analysis_point=pars['xsec'],pars['sow'],pars['passoptions'],pars['analysispoint']
    nEvents=pars['nEvents']
    from modules.paths import IHEP_path,golden_json_path    
    from coffea.lumi_tools import LumiData, LumiMask, LumiList
    from modules.hcoll import hcoll,binning
    hists = hcoll(h, isData, xsec, sow, doweight, process=histAxisName)
    ET.autolog(f"{len(ev)} Events at the start of your analysis", logger, 'i')
############################# KEEP THIS BLOCK ##################

    # Start your analysis``
    #-------------------------------------------------------------------------------------------------------
    lumidata=LumiData(f'modules/lumidata/lumi{year}.csv')

    ev["weight_norm"]=1/len(ev)
    #-------------------------------------------------------------------------------------------------------#### Di Muon events
    lumi=0
    if 'mask' in pass_options:
        lumilist=LumiMask(golden_json_path(year))(ev.run,ev.luminosityBlock)
        lumi=lumidata.get_lumi(LumiList(ev[lumilist].run,ev[lumilist].luminosityBlock))
    
    lumi=lumidata.get_lumi(LumiList(ev.run,ev.luminosityBlock))
    f= open('Analysis/LumiCal/lumi.txt',"a+")
    f.write(str(lumi)+'\n')
    #print(f'lumi is :{lumi}')
    #hists.fill('LumiCal',ev.weight_norm,(ev.event>0), ev, Lumi="lumi")
    #-------------------------------------------------------------------------------------------------------
    ############################# KEEP THIS BLOCK ##################
    ET.autolog(f"{len(ev)} Events at the end of your analysis", logger, 'i')
    return hists.get()
    ############################# KEEP THIS BLOCK ##################
