import awkward as ak
from modules.paths import IHEP_path,golden_json_path
from coffea.lumi_tools import LumiMask

def preselection(pars, events, selections):
    #-----------Add your pre selection here----------------------#
    year,isData,analysispoint=pars['year'],pars['isData'],pars['analysispoint']
    # if isData:
    #     lumi_mask = LumiMask(golden_json_path(year))(events.run,events.luminosityBlock)
    #     selections.add("is_good_lumi",lumi_mask)
    selections.add('passed HLT_DoubleIsoMu17_eta2p1', events.HLT.DoubleIsoMu17_eta2p1>0)
    selections.add("2 muons", (ak.num(events.recoMu) == 2))
    selections.add("2 photons", (ak.num(events.recoPho) == 1))

    #-----------Add your pre selection here----------------------#
    return events, selections


