from modules.ExpressoPlotter import ExpressoPlotter,normalplot

DY='Output/Analysis/HiggsToZG_Run2/output/analysis/DYJetsToLL_passop_Xsecweight.pkl.gz'
DM='Output/Analysis/HiggsToZG_Run2/output/analysis/DoubleMuonRun2016D_passop_Xsecweight.pkl.gz'
ZG='Output/Analysis/HiggsToZG_Run2/output/analysis/ZGToLLG_passop_Xsecweight.pkl.gz'

lumi=(4.4/1.5)*1000.0*(1075761.0 / 33861745.0)
print(f'lumi = {lumi/1000} fb-1')
plotter=ExpressoPlotter("2016")
plotter.histolocation('./')
plotter.savelocation('./Output/plots/')
plotter.settings('modules/plotsettings.yaml')
plotter.print_stat()
#plotter.noyerr()
plotter.plot_ratio()
plotter.addfile('Data',DM,'black','nostack',1,isdata=True) ## -1 normalizes to 1
plotter.addfile('DYJetsToLL',DY,'red','stack',lumi) ## -1 normalizes to 1
plotter.addfile('ZGToLLG',ZG,'blue','stack',lumi) ## -1 normalizes to 1


p=normalplot(plotter,filename="Mmmg",hi='Mmmg',axis='Mmmg',pklfiles='DYJetsToLL,ZGToLLG,Data',rebin=10)
p=normalplot(plotter,filename="Mmm",hi='Mmm',axis='Mmm',pklfiles='DYJetsToLL,ZGToLLG,Data',rebin=10)
p=normalplot(plotter,filename="PhotonpT",hi='PhotonpT',axis='PhotonpT',pklfiles='DYJetsToLL,ZGToLLG,Data',rebin=10)
p=normalplot(plotter,filename="LeadingMuonpT",hi='LeadingMuonpT',axis='LeadingMuonpT',pklfiles='DYJetsToLL,ZGToLLG,Data',rebin=10)
#p=normalplot(plotter,filename="cutflow_individual",hi='cutflow_individual',axis='x',pklfiles='DYJetsToLL,ZGToLLG,Data')
#p=normalplot(plotter,filename="cutflow",hi='cutflow',axis='x',pklfiles='DYJetsToLL,ZGToLLG,Data')
#p=normalplot(plotter,filename="sumw",hi='sumw',axis='sumw',pklfiles='DYJetsToLL,ZGToLLG,Data')
#p=normalplot(plotter,filename="sumw",hi='sumw',axis='sumw',pklfiles='DYJetsToLL,ZGToLLG,Data')
