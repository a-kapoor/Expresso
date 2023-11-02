'''
This code defines a class hcoll which is the main histogramming class of expresso. Every instance of hcoll contains a collection of hist histograms. The fill method of the class fills the histogram with the given data. The get method returns the histogram. The binning function returns a list of numbers from a to b with a step size of c. The __main__ block prints some information about the class and a message to see the testAnalysis folder to create your own barebones analysis.
'''
import numpy as np
import awkward as ak
# Define the hcoll class
class hcoll:
    # Define the __init__ method
    def __init__(self, h, isData, xsec, sow, doweight, **conf):
        self.h = h
        self.conf = conf
        self.isData = isData
        self.xsec = xsec
        self.sow = sow
        self.doweight=doweight
        
    # Define the fill method
    def fill(self, name, weights, mask, obj, cat={}, flatten=False, **axes):
        fullhist = {}
        # Loop over the axes
        #print('----------########---')
        for ini, axis in enumerate(axes.keys()):
            # Evaluate the object and mask
            arrr = eval(f"obj.{axes[axis]}[mask]")
            # If flatten is True
            if flatten:
                # If this is the first axis
                if ini==0:
                    # Flatten the array and weights
                    fullhist[axis],weights = ak.flatten(ak.zip(arrr,weights))
                    weights=weights[ak.flatten(mask)]
                else:
                    # Flatten only the array
                    fullhist[axis] = ak.flatten(arrr)
            else:
                # If this is the first axis
                if ini==0:
                    weights=weights[mask]
                fullhist[axis] = arrr
        # If doweight is True
        if self.doweight:
            # Fill the histogram with weights
            self.h[name].fill(weight=weights, **cat, **fullhist, **self.conf)
            
        else:
            # Fill the histogram without weights
            self.h[name].fill(**cat, **fullhist, **self.conf)

    def get(self):
        # Define the get method
        return self.h

def binning(a,b,c):
    # Define the binning function
    return list(np.arange(a,b+c,c))

# If this is the main block
if __name__=='__main__':

    print('''
        This the main histogramming class of expresso

    Every instance of hcoll contains a collection of hist histograms

    ''')

    print('''
        See testAnalysis folder to create your own barebones analysis
    ''')
