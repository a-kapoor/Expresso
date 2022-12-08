from modules.ExpressoPlotter import ExpressoPlotter,normalplot
import sys
import coffea, uproot3, numpy

print(f'Input:{sys.argv[1]}')

plotter=ExpressoPlotter('')
plotter.histolocation('./')
plotter.savelocation('./')
plotter.addfile('myfile',sys.argv[1],'','',1)

fout = uproot3.create(sys.argv[1].replace(".pkl.gz",".root"))
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
print(f'Done! Wrote: {sys.argv[1].replace(".pkl.gz",",root")}')

