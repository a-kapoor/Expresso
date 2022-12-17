from modules.ExpressoPlotter import ExpressoPlotter,normalplot
import sys
import coffea, uproot3, numpy
from modules.ExpressoTools import cprint
import os

print(f'Input:{sys.argv[1]}')

plotter=ExpressoPlotter('')
plotter.histolocation('./')
plotter.savelocation('./')
plotter.addfile('myfile',sys.argv[1],'','',1)

filename=sys.argv[1].replace(".pkl.gz",".root")
try:
    fout = uproot3.create(filename)
except:
    commandtodel=f'rm {filename}'
    os.system(commandtodel)
    cprint(f'Root file with same name already present: {filename}, Overwritten now.','OKCYAN')
    fout = uproot3.create(filename)
    exit()

print('Writing Histograms to root!')
for histname in list(plotter._files[0]['coffehists'].keys()):
    h=plotter._files[0]['coffehists'][histname]
    #print(h.fields)
    if 'process' in h.fields:
        h=h.integrate('process')
        #print(histname)
        if len(h.fields)<2:
            fout[histname] = coffea.hist.export1d(h.project(h.fields[0]))
        else:
            print('Only 1d histos supported')
    if 'selection' in h.fields:
        h=h.integrate('selection')
        #print(histname)
        if len(h.fields)<2:
            fout[histname] = coffea.hist.export1d(h.project(h.fields[0]))
        else:
            print('Only 1d histos supported')
    else:
        #print(histname)
        if len(h.fields)<2:
            fout[histname] = coffea.hist.export1d(h.project(h.fields[0]))
        else:
            print('Only 1d histos supported')
fout.close()
print(f'Done! Wrote: {sys.argv[1].replace(".pkl.gz",".root")}')

