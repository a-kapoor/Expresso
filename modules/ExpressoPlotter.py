import numpy as np
import gzip
from pathlib import Path
import yaml
import os
try:
    import pickle5 as pickle
except ImportError:
    import pickle
#import hist
import coffea.hist as hist
import hist as skhist
from tabulate import tabulate
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 200
plt.rcParams['savefig.dpi'] = 200
import mplhep as hep
from hist.intervals import ratio_uncertainty, clopper_pearson_interval, poisson_interval
import numpy as np

class ExpressoPlotter():
    def __init__(self,year):
        print('Expresso Plotter Initialized!')
        import numpy as np
        import gzip
        import yaml
        import os
        try:
            import pickle5 as pickle
        except ImportError:
            import pickle
        #import hist
        import coffea.hist as hist
        from tabulate import tabulate
        import matplotlib.pyplot as plt
        plt.rcParams['figure.dpi'] = 200
        plt.rcParams['savefig.dpi'] = 200
        import mplhep as hep
        from hist.intervals import ratio_uncertainty, clopper_pearson_interval, poisson_interval
        import numpy as np
        self._files=[]
        self._LABELS=[]
        self._year=year
        self._plots=[]
        self.lumi=0.0
        self.plot_live=False
        self.printmeanvar=False
        self.isnotebook=False
        self._yerr=True
        self.plotratio=False
        self._SSaveLocation='./'
        settingspath = Path('modules/plotsettings.yaml')
        with open(settingspath) as stream:
            self._ps=yaml.safe_load(stream)
    def notebook(self):
        self.isnotebook=True
    def noyerr(self):
        self._yerr=None
    def print_stat(self):
        self.printmeanvar=True
        self.printstat=True
    def plot_ratio(self):
        self.plotratio=True
        

    def settings(self,file):
        path = Path(file)
        with open(path) as stream:
            self._ps = yaml.safe_load(stream)

    def get_hist_from_pkl(self, path_to_pkl,allow_empty=True,tohist=False):
        h = pickle.load( gzip.open(path_to_pkl) )
        #print(h['Mmm'].values())

        if not allow_empty:
            h = {k:v for k,v in h.items() if v.values() != {}}
        if tohist:
            return h.to_hist()
        else:
            return h
        
    def geths(self,h,scale=-1):
        # Return scaled histogram
        return h*(1/sum(h.counts())) if scale==-1 else h*scale
    
    def dictprint(self, di):
        for key, value in di.items():
            print(key, ' : ', value)

    def geterrratio(self,hi,typeunc='p'):
        ratio = (hi[0].values()/hi[1].values())
        if typeunc=='p':
            err_down, err_up  = ratio_uncertainty(hi[0].values(), hi[1].values(), 'efficiency')
        else:
            err_down, err_up = clopper_pearson_interval(hi[0].values(), hi[1].values())
        labels=[]
        for ra, u, d in zip(ratio.ravel(), err_up.ravel(), err_down.ravel()):
            ra, u, d = f'{ra:.6f}', f'{u:.6f}', f'{d:.6f}'
            st = '$'+ra+'_{-'+d+'}^{+'+u+'}$'
            labels.append(st)
        return ratio,labels

    def histolocation(self,loc):
        self._loc=loc
    def savelocation(self,loc):
        # Set the save location to the given location
        self._SSaveLocation=loc
        # Check if the save location exists
        isExist = os.path.exists(self._SSaveLocation)
        # If the save location does not exist, create it
        if not isExist:
            os.makedirs(self._SSaveLocation)
            print("created folder : ", self._SSaveLocation)


    def addfile(self,label,thefile,color,stack=False,scale=1,isdata=False):
        if label in self._LABELS:
            print(f"Error: Two files can not have same name: {label}")
            exit()
        self._files.append({'label':label,'file':thefile,'color':color,'stack':stack,
                            'scale':scale,'isdata':isdata,
                            'coffehists':self.get_hist_from_pkl(self._loc+'/'+thefile),
                            #'h':self.get_hist_from_pkl(self._loc+'/'+file,tohist=True)
                           })
        self._LABELS.append(label)
        
    def gethist(self,filename,hi):
        file_a=None
        for plotterfile in self._files:
            if plotterfile['label']==filename:
                file_a=plotterfile['coffehists'][hi].to_hist()
        if file_a != None:
            return file_a
        else:
            print("wrong filename or hi or axis")
            return 0

    def getfigurewithaxis(self):
        if self.plotratio:
            fig, (ax1, ax2) = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios': [5, 1]}, constrained_layout=True)
        else:
            fig, ax1 = plt.subplots(1, constrained_layout=True)
            ax2 = plt.gca()
        return fig, ax1, ax2

    def getfiles(self):
        return self._files

class normalplot():
    def __init__(self,plotter,filename,hi,axis,rebin=1,colors='',pklfiles='',fontsize='xx-small',ncol=1,printplotnames=True,xlabel=None):
        print(f'Plotting {filename}')
        self.data=''
        self.his=hi
        self.axes=axis
        self.colors=colors
        self.plotter=plotter
        self.rebin=rebin
        self.printplotnames=printplotnames
        SSaveLocation=plotter._SSaveLocation
        hep.style.use(hep.style.ROOT) if plotter._ps['hepstyle'] == 'ROOT' else None
        hep.style.use("CMS") if plotter._ps['hepstyle']=='CMS' else None
        fig,ax1,ax2=plotter.getfigurewithaxis()
        hep.cms.label(plotter._ps["PrivateLabel"], data=plotter._ps["withdata"], year=plotter._year,ax=ax1,lumi=round(plotter.lumi/1000,3)) if plotter._ps['hepstyle']=='CMS' else None
        self.his=self.his.split(','); self.axes=self.axes.split(',');

        (nostack, stack, nostacklabels, nostackcolors, stacklabels,stackcolors,stackscales, nostackscales, stackerrors, nostackerrors,
         stack_isdatalist, nostack_isdatalist) = self.get_stacks_and_overlays(plotter,pklfiles)
    
        if len(stack)!=0:
            scaledstacks=[plotter.geths(st,scaleit) for st,scaleit in zip(stack,stackscales)]
            hep.histplot(scaledstacks,lw=1,ax=ax1,stack=True,histtype='fill',label=stacklabels, color=stackcolors)
            self.fullstack=sum(scaledstacks)
            try:
                ax1.fill_between(self.fullstack.axes[0].centers,#,ax=ax1,
                                 self.fullstack.values()-np.sqrt(self.fullstack.variances()),
                                 self.fullstack.values()+np.sqrt(self.fullstack.variances()),
                                 hatch='/////', zorder=2, color= 'none',edgecolor='k',step='mid',alpha=0.6,label='stat_unc')
            except:
                print("Not plotting errors")
                        
        if len(nostack)!=0:
            
            for nst,scaleit,labelit,colorit,isdatait in zip(nostack,nostackscales,nostacklabels,
                                                            nostackcolors,nostack_isdatalist):
                if not isdatait:
                    here_h=plotter.geths(nst,scaleit)
                    hep.histplot(here_h,ax=ax1,
                                 lw=1,stack=False,histtype='step',label=labelit,color=colorit,yerr=False)
                    
                    if plotter._yerr:
                        ax1.fill_between(here_h.axes[0].centers,
                                         here_h.values()-np.sqrt(here_h.variances()),
                                         here_h.values()+np.sqrt(here_h.variances()),
                                         hatch='', zorder=2, fc='grey',step='mid',alpha=0.6)
                    
                else:
                    if scaleit!=1:
                        print("Warning: You are scaling data, I hope you know what you are doing")
                    self.data=plotter.geths(nst,scaleit)
                    hep.histplot(self.data,ax=ax1,
                                 lw=2,stack=False,histtype='errorbar',label=labelit,color=colorit)

            plain_ratio=0        
            if plotter.plotratio and self.data:
                xlbl = ax1.get_xaxis().get_label()
                ratio = (self.data.values()/self.fullstack.values())
                err_down, err_up  = ratio_uncertainty(self.data.values(), self.fullstack.values())##Needs to be checked
                if plotter.printstat:
                    plain_ratio=round(sum(self.data.values())/sum(self.fullstack.values()),2)
                ax2.errorbar(self.data.axes[0].centers,ratio,yerr=(err_down, err_up),fmt='k.',markersize=8)
                if not xlabel:xlabel=xlbl.get_text()
                ax2.set_xlabel(xlabel);
                ax1.set_xlabel(None);
                ax2.set_ylim(0,2);
                ax2.axhline(y=1, color='black', linestyle='--')
            if plotter.printstat:
                ax1.legend(title=f'Data/Bkg= {plain_ratio}',loc='best',fontsize=fontsize,ncol=1,fancybox=True,title_fontsize=fontsize)
            else:
                ax1.legend(loc='best',fontsize=fontsize,ncol=1,fancybox=True)
                    
        if not plotter.plotratio and not self.data:
            if xlabel:
                ax1.set_xlabel(xlabel)
        #plt.tight_layout()
        plt.savefig(f'{SSaveLocation}/normal_{filename}.pdf', dpi=200)
        if plotter.isnotebook:
            plt.show()
            self.stackerrors=stackerrors
            self.hists=stack
            self.scaledhists=[plotter.geths(st,scaleit) for st,scaleit in zip(stack,stackscales)]
        else:
            plt.close()


    def get_stacks_and_overlays(self,plotter,pklfiles):
        # Pick the pkl files to plot and also get the labels
        LABELS=pklfiles.split(',') if pklfiles else [myfile['label'] for myfile in plotter._files]                
        listcolors=self.colors
        nostack, stack, nostacklabels, nostackcolors, stacklabels,stackcolors,stackscales, nostackscales = [], [], [], [], [], [], [], []
        stackerrors, nostackerrors, stack_isdatalist, nostack_isdatalist = [], [], [], []

        for file in plotter._files:
            ## Only plot pkl files that have been specified
            if file['label'] not in LABELS: continue
            
            ### fill list of colors in given as options
            listcolors= [file['color']]*len(self.his) if not self.colors else self.colors.split(',')
            
            for hi,axis,color in zip(self.his,self.axes,listcolors):
                _ch=file['coffehists'][hi]; _color=color; _stack=file['stack']; _scale=file['scale']; _isdata=file['isdata']
                project=_ch.project(axis).to_hist()[::skhist.rebin(self.rebin)]
                _label=file['label']+ f'({hi})' if self.printplotnames else file['label']
                
                #Print stat on plot if needed
                if plotter.printmeanvar:
                    meanvar=getmeanvar(h=plotter.geths(project,_scale))
                    mean,var,num=meanvar['mean'],meanvar['var'],meanvar['num']
                    _label = _label + f'( mean = {mean}, var = {var}, integral= {num} )'
                    
                if(_stack=='stack'):
                    with np.errstate(divide='ignore',invalid='ignore'):
                        stackerrors.append(np.nan_to_num(np.sqrt(project.variances())/project.values()))
                    stack_isdatalist.append(_isdata); stack.append(project); stacklabels.append(_label); stackcolors.append(_color); stackscales.append(_scale)
                if(_stack=='nostack'):
                    with np.errstate(divide='ignore',invalid='ignore'):
                        nostackerrors.append(np.nan_to_num(np.sqrt(project.variances())/project.values()))
                    nostack_isdatalist.append(_isdata);nostack.append(project); nostacklabels.append(_label); nostackcolors.append(_color); nostackscales.append(_scale)
                    
        return (nostack, stack, nostacklabels, nostackcolors, stacklabels,stackcolors,stackscales, nostackscales, stackerrors, nostackerrors, stack_isdatalist, nostack_isdatalist)

    
    def getstack(self):
        return self.hists,self.scaledhists,self.stackerrors

def getmeanvar(h):
    n, bins = h.to_numpy()
    #print(n)
    mids=0.5*(bins[1:] + bins[:-1])
    num=np.sum(n)
    mean=np.average(mids, weights=n).round(2)
    var=np.average((mids - mean)**2, weights=n).round(2)
    return {'mean':mean,'var':var,'num':num}
