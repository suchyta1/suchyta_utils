import matplotlib as _mpl
import matplotlib.pyplot as _plt
import os as _os


def Setup():
    dir = _os.path.dirname(_os.path.realpath(__file__))
    style = _ReadStyle(_os.path.join(dir,'custom-sytle.mpl'))
    style = _SetStyle(style)


def _ReadStyle(file):                                                                                     
    d = {}                                                                                               
    with open(file) as f:                                                                                
        for line in f:                                                                                   
            l = line.strip()                                                                             
            if (l=='') or (l[0]== '#') :                                                                 
                continue                                                                                 
                                                                                                         
            key, val = l.split(':')                                                                      
            k = key.strip()                                                                              
            v = val.strip()                                                                              
                                                                                                         
            if v[0] in ["'", '"']:                                                                       
                d[k] = v                                                                                 
            else:                                                                                        
                try:                                                                                     
                    d[k] = float(v)                                                                      
                except:                                                                                  
                    d[k] = v                                                                             
    return d                                                                                             

                                                                                                         
def _SetStyle(style):                                                                                     
    for k,v in style.iteritems():                                                                        
        _mpl.rcParams[k]=v 


def NTicks(ax, nxticks=None, nyticks=None):
    if nyticks is not None:
        ax.yaxis.set_major_locator( _mpl.ticker.MaxNLocator(nbins=(nyticks-1)) )
        ax.yaxis.set_major_locator( _mpl.ticker.MaxNLocator(nbins=(nyticks-1)) )
    if nxticks is not None:
        ax.xaxis.set_major_locator( _mpl.ticker.MaxNLocator(nbins=(nxticks-1)) )
        ax.xaxis.set_major_locator( _mpl.ticker.MaxNLocator(nbins=(nxticks-1)) )


def OffsetX(r, offset=0, log=False):
    if log:
        newr = np.log10(r)
    else:
        newr = copy.copy(r)
    
    newr = newr + offset
    if log:
        newr = np.power(10, newr)

    return newr


def LineSegment(ax=None, left=None, right=None, plotkwargs={}):
    if ax is None:
        fig, ax = _plt.subplots(1,1)

    ax.plot( [left[0], right[0]], [left[1], right[1]], **plotkwargs )
    return ax
