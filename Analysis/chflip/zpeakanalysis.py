#!/usr/bin/env python
import awkward as ak
import coffea as coffea
from coffea import hist
from modules.selection import *
import modules.ExpressoTools as ET
import modules.objects as obj
from modules.paths import IHEP_path
from coffea.analysis_tools import PackedSelection
#####################################################################################################################

def myanalysis(logger,h,ev,dataset,isData,histAxisName,year,xsec,sow):
    ET.autolog(f'{len(ev)} Events at the start of your analysis',logger,'i')

    #------------------------Start your analysis below this
    
    #----------
    el=ev.ele_fo
    mu=ev.mu_fo
    jets=ev.goodJets
    #----------    
    
    el["isTightLep"] = obj.tightSelElec(el.isFO, el.mvaTTH)
    el = el[el.isPres & el.isLooseE & el.isFO & el.isTightLep & (el.tightCharge>=2)]

    h["Nele"].fill(Nele=ak.num(el),sam=histAxisName)
    
    isflip = (el.gen_pdgId == -el.pdgId)
    noflip = (el.gen_pdgId == el.pdgId)
    isprompt = ((el.genPartFlav==1) | (el.genPartFlav == 15))
    truthFlip_mask   = ak.fill_none((isflip & isprompt),False)
    truthNoFlip_mask = ak.fill_none((noflip & isprompt),False)
    
    dense_objs_flat = ak.flatten(el[truthFlip_mask])
    h["ptabseta_flip"].fill(pt=dense_objs_flat.pt,abseta=abs(dense_objs_flat.eta),sam=histAxisName)

    dense_objs_flat = ak.flatten(el[truthNoFlip_mask])
    h["ptabseta_noflip"].fill(pt=dense_objs_flat.pt,abseta=abs(dense_objs_flat.eta),sam=histAxisName)

    #######################  Work with dielectron events ##########################################

    el=el[ak.num(el)==2]
    
    el0=el[:,0]
    el1=el[:,1]
    

    isflip = (el0.gen_pdgId == -el0.pdgId)
    noflip = (el0.gen_pdgId == el0.pdgId)
    isprompt = ((el0.genPartFlav==1) | (el0.genPartFlav == 15))
    truthFlip_mask   = ak.fill_none((isflip & isprompt),False)
    truthNoFlip_mask = ak.fill_none((noflip & isprompt),False)
    

    histmasks={"ptabseta_flip_el0":truthFlip_mask,
               "ptabseta_noflip_el0":truthNoFlip_mask}

    for histname in histmasks.keys():
        dense_objs_flat = el0[histmasks[histname]]
        h[histname].fill(pt=dense_objs_flat.pt,abseta=abs(dense_objs_flat.eta),sam=histAxisName)
    
    isflip = (el1.gen_pdgId == -el1.pdgId)
    noflip = (el1.gen_pdgId == el1.pdgId)
    isprompt = ((el1.genPartFlav==1) | (el1.genPartFlav == 15))
    truthFlip_mask   = ak.fill_none((isflip & isprompt),False)
    truthNoFlip_mask = ak.fill_none((noflip & isprompt),False)

    histmasks={"ptabseta_flip_el1":truthFlip_mask,
               "ptabseta_noflip_el1":truthNoFlip_mask}

    for histname in histmasks.keys():
        dense_objs_flat = el1[histmasks[histname]]
        h[histname].fill(pt=dense_objs_flat.pt,abseta=abs(dense_objs_flat.eta),sam=histAxisName)
    
    #------------------------End your analysis here

    ET.autolog(f'{len(ev)} Events at the end of your analysis',logger,'i')
    return h
    #######################  Define your analysis here ##########################################


histograms={
    "Nele" : hist.Hist(
        "Events",
        hist.Cat("sam", "sam"),
        hist.Bin("Nele", "Nele", [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]),
    ),
    "ptabseta_noflip" : hist.Hist(
        "Events",
        hist.Cat("sam", "sam"),
        hist.Bin("pt", "pt", [0, 30.0, 45.0, 60.0, 100.0, 200.0]),
        hist.Bin("abseta", "abseta", [0, 0.4, 0.8, 1.1, 1.4, 1.6, 1.9, 2.2, 2.5]),
    ),
    "ptabseta_flip" : hist.Hist(
        "Events",
        hist.Cat("sam", "sam"),
        hist.Bin("pt", "pt", [0, 30.0, 45.0, 60.0, 100.0, 200.0]),
        hist.Bin("abseta", "abseta", [0, 0.4, 0.8, 1.1, 1.4, 1.6, 1.9, 2.2, 2.5]),
    ),
    "ptabseta_noflip_el0" : hist.Hist(
        "Events",
        hist.Cat("sam", "sam"),
        hist.Bin("pt", "pt", [0, 30.0, 45.0, 60.0, 100.0, 200.0]),
        hist.Bin("abseta", "abseta", [0, 0.4, 0.8, 1.1, 1.4, 1.6, 1.9, 2.2, 2.5]),
    ),
    "ptabseta_flip_el0" : hist.Hist(
        "Events",
        hist.Cat("sam", "sam"),
        hist.Bin("pt", "pt", [0, 30.0, 45.0, 60.0, 100.0, 200.0]),
        hist.Bin("abseta", "abseta", [0, 0.4, 0.8, 1.1, 1.4, 1.6, 1.9, 2.2, 2.5]),
    ),
    "ptabseta_noflip_el1" : hist.Hist(
        "Events",
        hist.Cat("sam", "sam"),
        hist.Bin("pt", "pt", [0, 30.0, 45.0, 60.0, 100.0, 200.0]),
        hist.Bin("abseta", "abseta", [0, 0.4, 0.8, 1.1, 1.4, 1.6, 1.9, 2.2, 2.5]),
    ),
    "ptabseta_flip_el1" : hist.Hist(
        "Events",
        hist.Cat("sam", "sam"),
        hist.Bin("pt", "pt", [0, 30.0, 45.0, 60.0, 100.0, 200.0]),
        hist.Bin("abseta", "abseta", [0, 0.4, 0.8, 1.1, 1.4, 1.6, 1.9, 2.2, 2.5]),
    )
}
