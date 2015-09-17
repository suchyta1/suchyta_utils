import matplotlib as _mpl
import os as _os

def Setup():
    dir = _os.path.dirname(os.path.realpath(__file__))
    style = _ReadStyle(os.path.join(dir,'custom-sytle.mpl')
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
