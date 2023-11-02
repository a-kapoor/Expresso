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
    selections.add('Primary Vertex Selection', events.has_atleast_onegood_PV)
    selections.add('Lepton Selection', events.has_one_goodlepton)
    selections.add('MET Filters', events.passesMETfilters)
    selections.add('b-Jet Selection', events.bjet_selection)
    selections.add('W-Jet Selection', events.wjet_selection)
    selections.add('Trigger Selection', events.trigger_selection)
 
    #-----------Add your pre selection here----------------------#
    return events, selections


