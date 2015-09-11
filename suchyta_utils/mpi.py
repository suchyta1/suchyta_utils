"""
Each of these functions can be used with as many arguments as you want, e.g.:
    arr1, arr2, arr3 = Scatter(arr1, arr2, arr3)
"""

import _mpifunctions as _mpif

#__all__ = ['Gather']

def Broadcast(*args):
    """
    Send the same array(s) everywhere.
    """
    if len(args)==1:
        return _mpif._broadcast(args[0])
    else:
        return [_mpif._broadcast(arg) for arg in args]

def Gather(*args):
    """
    Gather array(s) together from various nodes/CPUss you're using. 
    Basically way to piece pack to togther results after you've done a Scatter.
    You gather the results along the 0th axis, and each element can be an aribrarily dimensioned array.
    """
    if len(args)==1:
        return _mpif._gather(args[0])
    else:
        return [_mpif._gather(arg) for arg in args]

def Scatter(*args):
    """
    Scatter array(s) along 0th dimension. In other words, you'll split up the array along the 0th axis to distribute the data to the nodes/cpus you're using.
    Each element along the 0th dimension can be an array of as many dimensions as you want.
    """
    if len(args)==1:
        return _mpif._scatter(args[0])
    else:
        return [_mpif._scatter(arg) for arg in args]
