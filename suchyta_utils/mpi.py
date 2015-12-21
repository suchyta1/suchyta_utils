"""
:mod:`suchyta_utils.mpi` is some functions for working with mpi4py.
Each of these functions can be used with as many arguments as you want, e.g.::

    arr1, arr2, arr3 = Scatter(arr1, arr2, arr3)

.. note::
    :mod:`suchyta_utils.mpi` is not imported by default. This is to avoid
    nested calls to mpi4py.MPI during import, if trying to use other :mod:`suchyta_utils` functionality 
    from a program launched by the MPI runner.
"""

import _mpifunctions as _mpif

#__all__ = ['Gather']

def Broadcast(*args):
    """
    Send the same array(s) everywhere.

    Parameters
    ----------
    *args
        Arguments to broadcast everywhere

    Returns
    -------
    sent
        The broadcasted objects(s)

    """
    if len(args)==1:
        return _mpif._broadcast(args[0])
    else:
        return [_mpif._broadcast(arg) for arg in args]

def Gather(*args):
    """
    Gather array(s) together from various CPUs you're using. 
    Basically, a way to piece back togther results after you've done a Scatter.
    Each item can be an aribrarily dimensioned array.

    Paramters
    ---------
    *args
        The array(s) to join back together into single array(s) of all the results

    Returns
    -------
    gathered
        The joined results

    """
    if len(args)==1:
        return _mpif._gather(args[0])
    else:
        return [_mpif._gather(arg) for arg in args]

def Scatter(*args):
    """
    Scatter array(s) along 0th dimension. 
    In other words, you'll split up the array along the 0th axis to distribute the data to the nodes/cpus you're using.
    Each element along the 0th dimension can be an array of as many dimensions as you want.

    Parameters
    ----------
    *args
        The array(s) to distribute

    Returns
    -------
    scattered
        The distributed array(s)

    """
    if len(args)==1:
        return _mpif._scatter(args[0])
    else:
        return [_mpif._scatter(arg) for arg in args]


