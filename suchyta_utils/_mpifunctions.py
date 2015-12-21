import sys
import itertools
import numpy as np
from mpi4py import MPI


def _scatter(arr):
    if arr is not None:
        arr = np.array_split(arr, MPI.COMM_WORLD.size, axis=0)
    arr = MPI.COMM_WORLD.scatter(arr, root=0)
    return arr

def _broadcast(arr):
    return MPI.COMM_WORLD.bcast(arr, root=0)

def _gather(arr):
    arr = MPI.COMM_WORLD.gather(arr, root=0)
    
    if MPI.COMM_WORLD.Get_rank()==0:
        arr = np.array(arr)
        
        """
        isarr = (type(arr).__module__ == np.__name__)
        try:
            len(arr[0])
            islen = True
        except:
            islen = False
        reshape = (isarr and islen)
        """

        #if (arr.ndim > 1) | (reshape):
        if (arr.ndim > 1) | ((type(arr).__module__ == np.__name__) and (type(arr[0]).__module__ == np.__name__)):
            dt = arr[0].dtype
            size = arr[0].ndim
            s = np.array(arr[0].shape)
            rs = None
            if len(s) > 1:
                rs = s[1:]

            for i in range(size):
                arr = itertools.chain.from_iterable(arr)

            #arr = np.fromiter(arr, dtype=np.float)
            arr = np.fromiter(arr, dtype=dt)
            
            if rs is not None:
                z = arr.shape[0] / np.prod(rs)
                ns = np.insert(rs, 0, z)
                arr = np.reshape(arr, ns)
            
    return arr


