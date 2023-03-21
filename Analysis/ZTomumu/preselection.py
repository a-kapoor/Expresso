import awkward as ak
from modules.paths import IHEP_path,golden_json_path
from coffea.lumi_tools import LumiMask

def preselection(pars, events, selections):
    #-----------Add your pre selection here----------------------#
    year,isData,analysispoint=pars['year'],pars['isData'],pars['analysispoint']
    '''
    if isData:
              lumi_mask = LumiMask(golden_json_path(year))(events.run,events.luminosityBlock)
              selections.add("is_good_lumi",lumi_mask)
    '''
    #selections.add('passed HLT_Mu50', events.HLT.Mu50>0)
    #selections.add('passed HLT_DoubleIsoMu17_eta2p1', events.HLT.DoubleIsoMu17_eta2p1>0)
    #selections.add('passed HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ', (events.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ>0) |(events.HLT.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ>0))
    selections.add('passed HLT_IsoMu27', events.HLT.IsoMu27>0)
    selections.add("2 muons", (ak.num(events.recoMu) == 2))
    #selections.add("leading Muon > 60", events.recoMu[:,0].pt>60)
    #selections.add("1 photon", (ak.num(events.recoPho) == 1))

    #-----------Add your pre selection here----------------------#
    return events, selections


