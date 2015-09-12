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
        if arr.ndim > 1:
            size = arr[0].ndim
            s = np.array(arr[0].shape)
            rs = None
            if len(s) > 1:
                rs = s[1:]

            for i in range(size):
                arr = itertools.chain.from_iterable(arr)
            arr = np.fromiter(arr, dtype=np.float)
            
            if rs is not None:
                z = arr.shape[0] / np.prod(rs)
                ns = np.insert(rs, 0, z)
                arr = np.reshape(arr, ns)
            
    return arr


