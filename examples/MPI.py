#!/usr/bin/env python

"""
This file contains examples for the three MPI functions available from suchyta_utils:
    Scatter: distribute job across various CPUs
    Gather: collect data from various CPUs into a single entity
    Broadcast: send the same piece of data everywhere
"""

import suchyta_utils as es
from mpi4py import MPI

import numpy as np
import matplotlib.pyplot as plt


if __name__ == "__main__":
   
    # 1D array being distributed to the various available processors, and then returned back
    size = 1e6
    if MPI.COMM_WORLD.Get_rank()==0:
       random = np.random.randn(size)
    else:
        random = None
    
    random = es.mpi.Scatter(random)
    random = random + MPI.COMM_WORLD.Get_rank() * 4
    random = es.mpi.Gather(random)

    if MPI.COMM_WORLD.Get_rank()==0:
        bins = np.arange(-4, 4*MPI.COMM_WORLD.size, 0.1)
        fig, ax = plt.subplots(1,1)
        ax.hist(random, bins=bins)
        fig.savefig('out-hist.png')



    # 2D array being distributed to the various available processors, and then returned back
    if MPI.COMM_WORLD.Get_rank()==0:
        arr = np.arange(0, 9*MPI.COMM_WORLD.size)
        arr = np.reshape(arr, (arr.shape[0]/3, 3))
        print 'scatter this: ', arr
    else:
        arr = None
    
    arr = es.mpi.Scatter(arr)
    arr = 2*arr
    print 'rank = %i, arr ='%(MPI.COMM_WORLD.Get_rank()), arr

    arr = es.mpi.Gather(arr)
    if MPI.COMM_WORLD.Get_rank()==0:
        print 'final array', arr


    # You can also Scatter/Gather multiple things in one command if you want.
    random, arr = es.mpi.Scatter(random, arr)
    random, arr = es.mpi.Gather(random, arr)



    # 1D array being sent everywhere
    if MPI.COMM_WORLD.Get_rank()==0:
        r = np.random.randn(10)
    else:
        r = None
    rr = es.mpi.Broadcast(r)
    print 'rank = %i, common array ='%(MPI.COMM_WORLD.Get_rank()), rr


    # You can use Broadcast with multiple arguments as well
    r1, r2 = es.mpi.Broadcast(r, r)
