from modules.ExpressoPlotter import ExpressoPlotter,normalplot

DY1='Output/Analysis/ZTomumu/output/analysis/DYJetsToLL_passop_Xsecweight.pkl.gz'
DM='Output/Analysis/ZTomumu/output/analysis/DoubleMuonRun2016D_passop_Xsecweight.pkl.gz'
ZG='Output/Analysis/ZTomumu/output/analysis/ZGToLLG_passop_Xsecweight.pkl.gz'
DY='Output/Analysis/ZTomumu/output/analysis//DYJetsToLL_M-50_madgraph_passop_Xsecweight.pkl.gz'
ZG='Output/Analysis/ZTomumu/output/analysis//ZGToLLG_passop_Xsecweight.pkl.gz'
#lumi=(4.4/1.5)*1000.0*(1075761.0 / 33861745.0)
lumi=128.130400
#lumi=1.0
print(f'lumi = {lumi/1000} fb-1')
plotter=ExpressoPlotter("2016")
plotter.histolocation('./')
plotter.savelocation('./Output/plots/')
plotter.settings('modules/plotsettings.yaml')
plotter.print_stat()
#plotter.noyerr()
plotter.plot_ratio()
plotter.addfile('Data',DM,'black','nostack',1,isdata=True) ## -1 normalizes to 1
plotter.addfile('DYJetsToLL_amc',DY1,'blue','nostack',lumi) ## -1 normalizes to 1
plotter.addfile('DYJetsToLL',DY,'red','stack',lumi) ## -1 normalizes to 1
plotter.addfile('ZGToLLG',ZG,'green','stack',lumi) ## -1 normalizes to 1


p=normalplot(plotter,filename="Mmmg",hi='Mmmg',axis='Mmmg',pklfiles='DYJetsToLL,ZGToLLG,Data',rebin=1)
p=normalplot(plotter,filename="Mmm",hi='Mmm',axis='Mmm',pklfiles='DYJetsToLL,ZGToLLG,DYJetsToLL_amc,Data',rebin=1)
p=normalplot(plotter,filename="Mmm_zoom",hi='Mmm_zoom',axis='Mmm',pklfiles='DYJetsToLL,ZGToLLG,Data',rebin=1)
p=normalplot(plotter,filename="PhotonpT",hi='PhotonpT',axis='PhotonpT',pklfiles='DYJetsToLL,ZGToLLG,Data',rebin=2)
p=normalplot(plotter,filename="LeadingMuonpT",hi='LeadingMuonpT',axis='LeadingMuonpT',pklfiles='DYJetsToLL,DYJetsToLL_amc,ZGToLLG,Data',rebin=4)
#p=normalplot(plotter,filename="cutflow_individual",hi='cutflow_individual',axis='x',pklfiles='DYJetsToLL,ZGToLLG,Data')
#p=normalplot(plotter,filename="cutflow",hi='cutflow',axis='x',pklfiles='DYJetsToLL,ZGToLLG,Data')
#p=normalplot(plotter,filename="sumw",hi='sumw',axis='sumw',pklfiles='DYJetsToLL,ZGToLLG,Data')
#p=normalplot(plotter,filename="sumw",hi='sumw',axis='sumw',pklfiles='DYJetsToLL,ZGToLLG,Data')
