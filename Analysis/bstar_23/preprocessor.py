from asyncio import events
from modules.corrections import SFevaluator, GetBTagSF, ApplyJetCorrections, GetBtagEff, AttachMuonSF, AttachElectronSF, AttachPerLeptonFR, GetPUSF, ApplyRochesterCorrections, ApplyJetSystematics, AttachPSWeights, AttachPdfWeights, AttachScaleWeights, GetTriggerSF
from modules.objects import isTightTau
from modules.cut import isClean
import awkward as ak

def preprocess(pars, events, AttachSF=False):
    import modules.ExpressoTools as ET
    import awkward as ak
    import numpy as np
    ###################################
    dataset,isData,histAxisName,year,xsec,sow=pars['dataset'],pars['isData'],pars['histAxisName'],pars['year'],pars['xsec'],pars['sow']
    #events['nEvents']=nEvents
    ################################### Start writing your preprocessor below -> for creating new branches on top of the ones already present in the NanoAOD
    ## some variables are auto available picked up from the sample json
    ## The isData variable is true or false depending on if you are operating on data or mc
    ## The xsec has the xsec value, sow has sum of gen weights
    ## histAxisName is the name of the sample from json (will be used as legends for plots)
    
    events["recoMu"]=events.Muon[(events.Muon.tightId==True) & (events.Muon.pt > 53) & (abs(events.Muon.eta)<2.4)]

    events["vetoMu"]=events.Muon[(events.Muon.looseId==True) & (events.Muon.pt > 30) & (abs(events.Muon.eta)<2.4)]

    vidNestedWPBitmap_withoutiso= (1<<2) | (1<<5) | (1<<8) | (1<<11) | (1<<14) | (1<<17) | (1<<20) | (1<<26) |(1<<29)
    events["Electron","tight_withoutiso"]= (vidNestedWPBitmap_withoutiso == events.Electron.vidNestedWPBitmap & vidNestedWPBitmap_withoutiso)
    events["recoEle"]=events.Electron[(events.Electron.tight_withoutiso) & (events.Electron.pt > 53) & (abs(events.Electron.eta)<2.4)] ## Remove isolation requirements

    events["vetoEle"]=events.Electron[(events.Electron.cutBased>1) & (events.Electron.pt > 30) & (abs(events.Electron.eta)<2.4)]
    
    events["jets"] = events.Jet[(events.Jet.pt > 30) & (abs(events.Jet.eta) < 2.4)]
    events["fatjets"] = events.FatJet[(events.FatJet.pt > 200) & (abs(events.FatJet.eta) < 2.4)]
    
    events["recoMu","dr_jets"]=ak.firsts(events.recoMu.metric_table(events.jets), axis=-1)
    events["recoEle","dr_jets"]=ak.firsts(events.recoEle.metric_table(events.jets), axis=-1)

    events["leptons"] = ak.concatenate([events["recoMu"], events["recoEle"]], axis=1)

    events["jets","dr_lep"]=ak.firsts(events.jets.metric_table(events.leptons), axis=-1)
    events["fatjets","dr_lep"]=ak.firsts(events.fatjets.metric_table(events.leptons), axis=-1)

    events["jets"] = events.jets[events["jets"].dr_lep > 0.4]
    events["fatjets"] = events.fatjets[events["fatjets"].dr_lep > 0.8]

    # Apply b-tagging
    events["bjets"] = events["jets"][events["jets"].btagDeepFlavB > 0.6321]

    events["wjets"] = events["fatjets"][(events["fatjets"].msoftdrop > 65) & (events["fatjets"].msoftdrop < 105)]

    events["bjets","dr_lep"]=ak.firsts(events.bjets.metric_table(events.leptons), axis=-1)
    events["bjets","dr_wjets"]=ak.firsts(events.bjets.metric_table(events.wjets), axis=-1)
    
    events["has_atleast_onegood_PV"]=events.PV.npvsGood > 0
    events["has_one_goodlepton"]=ak.num(events["leptons"])==1
    
    events["passesMETfilters"] =(
        (events.Flag.goodVertices) &
        (events.Flag.globalSuperTightHalo2016Filter) &
        (events.Flag.HBHENoiseFilter) &
        (events.Flag.HBHENoiseIsoFilter) &
        (events.Flag.EcalDeadCellTriggerPrimitiveFilter) &
        (events.Flag.BadPFMuonFilter) &
        (events.Flag.eeBadScFilter) &
        (events.Flag.ecalBadCalibFilter)
    )

    events["bjet_selection"] = (ak.num(events["bjets"]) >= 1 & (ak.num(events["bjets"].dr_lep < 2) >= 1))
    events["wjet_selection"] = (ak.num(events["wjets"]) == 1)
    
    events["trigger_selection"] = (
        events.HLT["Ele35_WPTight_Gsf"] | events.HLT["Photon200"] | events.HLT["Ele115_CaloIdVT_GsfTrkIdT"] | events.HLT["Mu50"]
    )

    ################################### Keep the return line as is
    return events
